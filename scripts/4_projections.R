library(here)
library(sf)
library(terra)
library(ncdf4)
library(ggplot2)
library(patchwork)

bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))

# Dropping empty polygons
bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]
UrbanAreas <- st_read(here("data", "2020_Arizona_Census_Urban_Areas", "2020_Arizona_Census_Urban_Areas.shp"))
UrbanAreas <- st_transform(UrbanAreas, st_crs(bg_sf))
UrbanAreas <- st_make_valid(UrbanAreas)
bg_centroids <- st_centroid(bg_sf)
inside <- st_intersects(bg_centroids, st_union(UrbanAreas), sparse = FALSE)[, 1]
bg_sf <- bg_sf[inside, ]
bg <- data.frame(bg_sf)

colnames(bg)[colnames(bg) == "geoid20"]   <- "GEOID"
colnames(bg)[colnames(bg) == "evp_prp"]   <- "evap_prop"
colnames(bg)[colnames(bg) == "med_inc"]    <- "med_income"
colnames(bg)[colnames(bg) == "pct_mnr"]    <- "pct_minority"
colnames(bg)[colnames(bg) == "ave_age"]    <- "med_year_built"
colnames(bg)[colnames(bg) == "pct_rnt"]    <- "pct_renter"
colnames(bg)[colnames(bg) == "pct_sfr"]    <- "pct_sfh"
colnames(bg)[colnames(bg) == "covennt"]    <- "covenant"

bg_sf <- st_as_sf(bg)
bg_sf <- st_transform(bg_sf, 4326)
coords <- st_coordinates(st_centroid(bg_sf))
bg$lon <- coords[, 1]
bg$lat <- coords[, 2]
bg_sf <- st_make_valid(bg_sf)

# Baseline
wx <- read.csv(here("data", "weather", "tucson_hourly.csv"))
wx$date <- as.Date(wx$date)
wx$temp_c <- (wx$tmpf - 32) * 5 / 9
wx$month <- as.integer(format(wx$date, "%m"))

calc_twet <- function(temp_c, rh_pct) {
  temp_c * atan(0.151977 * (rh_pct + 8.313659)^0.5) +
    atan(temp_c + rh_pct) -
    atan(rh_pct - 1.676331) +
    0.00391838 * rh_pct^1.5 * atan(0.023101 * rh_pct) -
    4.686035
}

calc_supply <- function(temp_c, rh_pct, eta = 0.65) {
  twet <- calc_twet(temp_c, rh_pct)
  temp_c - eta * (temp_c - twet)
}

eta_values <- c(0.5, 0.6, 0.7, 0.8)

wx$failure <- calc_supply(wx$temp_c, wx$relh) > 27
baseline_by_eta <- lapply(eta_values, function(eta_val) {
  fail <- calc_supply(wx$temp_c, wx$relh, eta = eta_val) > 27
  daily_f <- aggregate(fail ~ wx$date, FUN = sum)
  colnames(daily_f) <- c("date", "failure_hours")
  daily_f$year <- as.integer(format(daily_f$date, "%Y"))
  daily_f$month <- as.integer(format(daily_f$date, "%m"))
  daily_f$failure_day <- as.integer(daily_f$failure_hours > 0)
  summer <- daily_f[daily_f$month >= 6 & daily_f$month <= 9, ]
  annual <- aggregate(failure_day ~ year, data = summer, FUN = sum)
  data.frame(eta = eta_val, baseline = mean(annual$failure_day))
})
baseline_by_eta <- do.call(rbind, baseline_by_eta)

daily_obs <- aggregate(failure ~ date, data = wx, FUN = sum)
colnames(daily_obs)[2] <- "failure_hours"
daily_obs$year <- as.integer(format(as.Date(daily_obs$date), "%Y"))
daily_obs$month <- as.integer(format(as.Date(daily_obs$date), "%m"))
daily_obs$failure_day <- as.integer(daily_obs$failure_hours > 0)

obs_summer <- daily_obs[daily_obs$month >= 6 & daily_obs$month <= 9, ]
obs_annual <- aggregate(failure_day ~ year, data = obs_summer, FUN = sum)
baseline_failure_days <- mean(obs_annual$failure_day)
baseline_failure_se <- sd(obs_annual$failure_day) / sqrt(nrow(obs_annual))
baseline_failure_min <- min(obs_annual$failure_day)
baseline_failure_max <- max(obs_annual$failure_day)

cat("\n=== OBSERVED BASELINE (eta = 0.65) ===\n")
cat("Mean annual failure days:", round(baseline_failure_days, 1), "\n")
cat("SE:", round(baseline_failure_se, 1), "\n")
cat("Range:", baseline_failure_min, "-", baseline_failure_max, "\n")

# NEX-GDDP-CMIP6 projections
dir.create(here("data", "cmip6"), recursive = TRUE, showWarnings = FALSE)
projections_path <- here("data", "cmip6", "tucson_projections.csv")

if (!file.exists(projections_path)) {
  cat("=== Downloading NEX-GDDP-CMIP6 projections ===\n")
  
  tuc_lat <- 32.25
  tuc_lon <- -111.0
  models <- c("ACCESS-CM2", "GFDL-ESM4", "IPSL-CM6A-LR",
              "MPI-ESM1-2-HR", "UKESM1-0-LL")
  ssps <- c("ssp245", "ssp585")
  periods <- list(
    historical = 2005:2014,
    near       = 2025:2034,
    mid        = 2040:2049,
    far        = 2060:2069
  )
  
  variables <- c("tasmax", "hurs")
  s3_base <- "https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com/NEX-GDDP-CMIP6"
  
  all_results <- list()
  
  for (mod in models) {
    for (ssp in ssps) {
      for (period_name in names(periods)) {
        yrs <- periods[[period_name]]
        experiment <- ifelse(period_name == "historical", "historical", ssp)
        
        for (yr in yrs) {
          for (v in variables) {
            
            fname <- sprintf("%s_day_%s_%s_r1i1p1f1_gn_%d.nc", v, mod, experiment, yr)
            url <- sprintf("%s/%s/%s/r1i1p1f1/%s/%s", s3_base, mod, experiment, v, fname)
            
            tmp <- tempfile(fileext = ".nc")
            result <- try(download.file(url, tmp, mode = "wb", quiet = TRUE), silent = TRUE)
            
            if (!inherits(result, "try-error") && file.size(tmp) > 10000) {
              nc_try <- try({
                r <- rast(tmp)
                r <- rotate(r)
                cell_id <- cellFromXY(r, cbind(tuc_lon, tuc_lat))
                vals <- r[cell_id]
                vals <- as.numeric(vals)
                
                is_leap <- (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
                jun1 <- ifelse(is_leap, 153, 152)
                sep30 <- ifelse(is_leap, 274, 273)
                
                if (length(vals) >= sep30) {
                  summer_vals <- vals[jun1:sep30]
                  dates <- seq(as.Date(paste0(yr, "-06-01")),
                               as.Date(paste0(yr, "-09-30")), by = "day")
                  if (length(summer_vals) == length(dates)) {
                    d <- data.frame(
                      model = mod, ssp = ssp, period = period_name,
                      year = yr, date = dates, variable = v,
                      value = summer_vals
                    )
                    all_results[[length(all_results) + 1]] <- d
                  }
                }
              }, silent = TRUE)
              unlink(tmp)
            }
          }
          cat("  ", mod, experiment, yr, "\n")
        }
      }
    }
  }
  proj_raw <- do.call(rbind, all_results)
  write.csv(proj_raw, projections_path, row.names = FALSE)
  cat("Projections saved:", nrow(proj_raw), "rows\n")
}

proj <- read.csv(projections_path, stringsAsFactors = FALSE)
proj$date <- as.Date(proj$date)
cat("Projection rows loaded:", nrow(proj), "\n")

# Cooler failure days
proj_tmax <- proj[proj$variable == "tasmax", c("model", "ssp", "period", "year", "date", "value")]
proj_hurs <- proj[proj$variable == "hurs", c("model", "ssp", "period", "year", "date", "value")]
colnames(proj_tmax)[6] <- "tasmax_K"
colnames(proj_hurs)[6] <- "hurs"

proj_wide <- merge(proj_tmax, proj_hurs,
                   by = c("model", "ssp", "period", "year", "date"))

# Compute model historical bias in hurs
hist_model_rh <- mean(proj_wide$hurs[proj_wide$period == "historical"], na.rm = TRUE)
hist_obs_rh <- mean(wx$relh[wx$month %in% 6:9], na.rm = TRUE)
rh_bias <- hist_model_rh - hist_obs_rh
cat("RH bias (model - observed):", round(rh_bias, 1), "\n")

proj_wide$hurs_corrected <- proj_wide$hurs - rh_bias
proj_wide$tmax_f <- proj_wide$tasmax_K * 9/5 - 459.67
proj_wide$tmax_c <- (proj_wide$tmax_f - 32) / 1.8

for (eta_val in eta_values) {
  col <- paste0("failure_eta", eta_val * 10)
  proj_wide[[col]] <- calc_supply(proj_wide$tmax_c, proj_wide$hurs_corrected, eta = eta_val) > 27
}

proj_wide$failure <- calc_supply(proj_wide$tmax_c, proj_wide$hurs_corrected) > 27
annual_fail <- aggregate(failure ~ model + ssp + period + year,
                         data = proj_wide, FUN = sum)
colnames(annual_fail)[5] <- "failure_days"
period_fail <- aggregate(failure_days ~ model + ssp + period,
                         data = annual_fail, FUN = mean)

period_fail_by_eta <- lapply(eta_values, function(eta_val) {
  col <- paste0("failure_eta", eta_val * 10)
  af <- aggregate(proj_wide[[col]] ~ model + ssp + period + year,
                  data = proj_wide, FUN = sum)
  colnames(af) <- c("model", "ssp", "period", "year", "failure_days")
  pf <- aggregate(failure_days ~ model + ssp + period, data = af, FUN = mean)
  pf$eta <- eta_val
  pf
})
period_fail_all_eta <- do.call(rbind, period_fail_by_eta)

# Delta method
hist_fail <- period_fail[period_fail$period == "historical", ]
colnames(hist_fail)[3] <- "hist_period"
colnames(hist_fail)[4] <- "hist_failure"
future_fail <- period_fail[period_fail$period != "historical", ]
future_fail <- merge(future_fail, hist_fail[, c("model", "hist_failure")], by = "model")
future_fail$delta <- future_fail$failure_days - future_fail$hist_failure

delta_summary <- aggregate(delta ~ ssp + period, data = future_fail,
                           FUN = function(x) c(median = median(x),
                                               p10 = quantile(x, 0.1),
                                               p90 = quantile(x, 0.9)))
delta_df <- data.frame(
  ssp = delta_summary$ssp,
  period = delta_summary$period,
  delta_median = delta_summary$delta[, 1],
  delta_p10 = delta_summary$delta[, 2],
  delta_p90 = delta_summary$delta[, 3]
)

cat("\n=== DELTA SUMMARY ===\n")
print(delta_df)

write.csv(delta_df, here("results", "P4_delta_summary.csv"), row.names = FALSE)

# Table 1: Delta-corrected projected failure days
table1_corrected <- delta_df
table1_corrected$proj_median <- round(baseline_failure_days + delta_df$delta_median, 1)
table1_corrected$proj_lo <- round(baseline_failure_days + delta_df$delta_p10, 1)
table1_corrected$proj_hi <- round(baseline_failure_days + delta_df$delta_p90, 1)
table1_corrected$label <- paste0(table1_corrected$proj_median,
                                 " [", table1_corrected$proj_lo, "-", table1_corrected$proj_hi, "]")

table1_corrected <- rbind(
  data.frame(ssp = "observed", period = "baseline",
             delta_median = 0, delta_p10 = 0, delta_p90 = 0,
             proj_median = round(baseline_failure_days, 1),
             proj_lo = round(baseline_failure_days, 1),
             proj_hi = round(baseline_failure_days, 1),
             label = as.character(round(baseline_failure_days, 1))),
  table1_corrected
)

cat("\n=== TABLE 1: Delta-corrected projected failure days ===\n")
print(table1_corrected[, c("ssp", "period", "label")])

write.csv(table1_corrected, here("results", "P4_Table1_projected_failure_days.csv"), row.names = FALSE)

# Sensitivity table: raw model failure days across eta values (supplement)
eta_proj_sensitivity <- lapply(eta_values, function(eta_val) {
  pf <- period_fail_all_eta[period_fail_all_eta$eta == eta_val, ]
  rows <- list()
  for (s in c("ssp245", "ssp585")) {
    for (p in c("historical", "near", "mid", "far")) {
      sub <- pf[pf$ssp == s & pf$period == p, ]
      if (nrow(sub) == 0 && p == "historical") {
        sub <- pf[pf$period == "historical", ]
      }
      if (nrow(sub) > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          eta = eta_val, ssp = s, period = p,
          median = round(median(sub$failure_days), 1),
          p10 = round(quantile(sub$failure_days, 0.1), 1),
          p90 = round(quantile(sub$failure_days, 0.9), 1))
      }
    }
  }
  do.call(rbind, rows)
})
eta_proj_table <- do.call(rbind, eta_proj_sensitivity)
eta_proj_table$label <- paste0(eta_proj_table$median, " [", eta_proj_table$p10, "-", eta_proj_table$p90, "]")
write.csv(eta_proj_table, here("results", "P4_Table_eta_sensitivity.csv"), row.names = FALSE)

# Project to block groups using delta method
for (i in seq_len(nrow(delta_df))) {
  s <- delta_df$ssp[i]
  p <- delta_df$period[i]
  d <- delta_df$delta_median[i]
  col_name <- paste0("exposure_", s, "_", p)
  bg[[col_name]] <- (baseline_failure_days + d) * bg$evap_prop
}

# Threshold crossing
threshold <- 10
bg$exposure_current <- baseline_failure_days * bg$evap_prop
bg$above_now <- bg$exposure_current > threshold

crossing_summary <- data.frame()
for (s in c("ssp245", "ssp585")) {
  for (p in c("near", "mid", "far")) {
    col <- paste0("exposure_", s, "_", p)
    above <- bg[[col]] > threshold
    newly_above <- above & !bg$above_now
    crossing_summary <- rbind(crossing_summary, data.frame(
      ssp = s, period = p,
      n_above = sum(above),
      n_newly_above = sum(newly_above),
      pct_above = round(mean(above) * 100, 1),
      pct_newly_above = round(mean(newly_above) * 100, 1)
    ))
  }
}

crossing_summary <- rbind(
  data.frame(ssp = "observed", period = "baseline",
             n_above = sum(bg$above_now),
             n_newly_above = 0,
             pct_above = round(mean(bg$above_now) * 100, 1),
             pct_newly_above = 0),
  crossing_summary
)

cat("\n=== THRESHOLD CROSSING ===\n")
print(crossing_summary)

write.csv(crossing_summary, here("results", "P4_Table2_threshold_crossing.csv"), row.names = FALSE)

# Demographics of newly vulnerable blocks
col_far <- "exposure_ssp585_far"
newly_vuln <- bg[bg[[col_far]] > threshold & !bg$above_now, ]
already_vuln <- bg[bg$above_now, ]
not_vuln <- bg[!bg$above_now & bg[[col_far]] <= threshold, ]

compare_vars <- c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built")
compare_labels <- c("Evap. prevalence", "Median income ($)", "Minority (%)",
                    "Renter (%)", "Year built")

tableS1 <- data.frame(
  variable = compare_labels,
  newly_mean = sapply(compare_vars, function(v) round(mean(newly_vuln[[v]], na.rm = TRUE), 2)),
  already_mean = sapply(compare_vars, function(v) round(mean(already_vuln[[v]], na.rm = TRUE), 2)),
  not_vuln_mean = sapply(compare_vars, function(v) round(mean(not_vuln[[v]], na.rm = TRUE), 2)),
  p_newly_vs_not = sapply(compare_vars, function(v) {
    if (nrow(newly_vuln) > 2 & nrow(not_vuln) > 2) {
      signif(wilcox.test(newly_vuln[[v]], not_vuln[[v]])$p.value, 3)
    } else { NA }
  })
)

cat("\n=== NEWLY VULNERABLE DEMOGRAPHICS ===\n")
print(tableS1)

write.csv(tableS1, here("results", "P4_TableS1_newly_vulnerable_demographics.csv"), row.names = FALSE)

# ============================================================
# FIGURES
# ============================================================

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

# ---- Figure 1a: Delta-corrected fan chart ----
ribbon_df <- delta_df
ribbon_df$median <- baseline_failure_days + delta_df$delta_median
ribbon_df$lo <- baseline_failure_days + delta_df$delta_p10
ribbon_df$hi <- baseline_failure_days + delta_df$delta_p90
ribbon_df$period_year <- c(2030, 2045, 2065, 2030, 2045, 2065)[
  match(paste(ribbon_df$ssp, ribbon_df$period),
        paste(rep(c("ssp245", "ssp585"), each = 3), rep(c("near", "mid", "far"), 2)))]

fig1 <- ggplot() +
  # SSP2-4.5 ribbon
  geom_ribbon(data = ribbon_df[ribbon_df$ssp == "ssp245", ],
              aes(x = period_year, ymin = lo, ymax = hi),
              fill = "#2c7bb6", alpha = 0.2) +
  geom_line(data = ribbon_df[ribbon_df$ssp == "ssp245", ],
            aes(x = period_year, y = median), color = "#2c7bb6", linewidth = 1) +
  # SSP5-8.5 ribbon
  geom_ribbon(data = ribbon_df[ribbon_df$ssp == "ssp585", ],
              aes(x = period_year, ymin = lo, ymax = hi),
              fill = "#d73027", alpha = 0.2) +
  geom_line(data = ribbon_df[ribbon_df$ssp == "ssp585", ],
            aes(x = period_year, y = median), color = "#d73027", linewidth = 1) +
  # Observed baseline with eta uncertainty
  geom_point(aes(x = 2020, y = baseline_failure_days), size = 3, shape = 18) +
  geom_errorbar(aes(x = 2020,
                    ymin = min(baseline_by_eta$baseline),
                    ymax = max(baseline_by_eta$baseline)),
                width = 0.25, color = "grey30", linewidth = 0.4) +
  # Labels
  annotate("text", x = 2050, y = max(ribbon_df$hi[ribbon_df$ssp == "ssp585"]),
           label = "SSP5-8.5", color = "#d73027", size = 3, hjust = 0) +
  annotate("text", x = 2050, y = max(ribbon_df$hi[ribbon_df$ssp == "ssp245"]),
           label = "SSP2-4.5", color = "#2c7bb6", size = 3, hjust = 0) +
  scale_x_continuous(breaks = c(2020, 2030, 2045, 2065),
                     labels = c("Obs.", "2030s", "2040s", "2060s")) +
  ylim(0, max(ribbon_df$hi) * 1.05) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "Cooler failure days per summer", title = "a")

# ---- Figure 1b: Map of projected compound exposure ----
bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = exposure_ssp585_far), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Compound\nexposure\n(2060s, SSP5-8.5)") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b")

# ---- Figure 1c: Threshold crossing bar chart ----
cross_plot <- crossing_summary[crossing_summary$ssp != "observed", ]
cross_plot$period <- factor(cross_plot$period, levels = c("near", "mid", "far"))
cross_plot$ssp_label <- ifelse(cross_plot$ssp == "ssp245", "SSP2-4.5", "SSP5-8.5")

fig3 <- ggplot(cross_plot, aes(x = period, y = pct_above, fill = ssp_label)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.8) +
  geom_hline(yintercept = crossing_summary$pct_above[crossing_summary$ssp == "observed"],
             linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("SSP2-4.5" = "#2c7bb6", "SSP5-8.5" = "#d73027"),
                    name = "") +
  scale_x_discrete(labels = c("2030s", "2040s", "2060s")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.85),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "% block groups above threshold", title = "c")

# ---- Figure S1: Current vs future exposure maps ----
figS1a <- ggplot(bg_sf2) +
  geom_sf(aes(fill = exposure_current), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Current", limits = range(c(bg$exposure_current, bg$exposure_ssp585_far))) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a  Current")

figS1b <- ggplot(bg_sf2) +
  geom_sf(aes(fill = exposure_ssp585_far), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "2060s SSP5-8.5", limits = range(c(bg$exposure_current, bg$exposure_ssp585_far))) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "b  SSP5-8.5, 2060s")

# ---- Figure S2: Eta sensitivity of projections ----
eta_ssp585_far <- eta_proj_table[eta_proj_table$ssp == "ssp585", ]

figS2 <- ggplot(eta_ssp585_far, aes(x = factor(period, levels = c("historical", "near", "mid", "far")),
                                    group = eta, color = factor(eta))) +
  geom_line(aes(y = median), linewidth = 0.8) +
  geom_point(aes(y = median), size = 2) +
  geom_ribbon(aes(ymin = p10, ymax = p90, fill = factor(eta)), alpha = 0.08, color = NA) +
  scale_color_manual(values = c("0.5" = "#2c7bb6", "0.6" = "#5b8fa8",
                                "0.7" = "#fee090", "0.8" = "#fc8d59"),
                     name = "eta") +
  scale_fill_manual(values = c("0.5" = "#2c7bb6", "0.6" = "#5b8fa8",
                               "0.7" = "#fee090", "0.8" = "#fc8d59"),
                    name = "eta") +
  scale_x_discrete(labels = c("Hist.", "2030s", "2040s", "2060s")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "", y = "Cooler failure days per summer (SSP5-8.5)", title = "")

# ---- Save all figures ----
pdf(here("results", "P4_FigureS2_eta_sensitivity.pdf"), width = 7, height = 5)
print(figS2)
dev.off()

png(here("results", "P4_FigureS2_eta_sensitivity.png"), width = 7, height = 5, units = "in", res = 300)
print(figS2)
dev.off()

fig_main <- (fig1 | fig2) / (fig3 + plot_spacer()) + plot_layout(heights = c(1, 0.7))

pdf(here("results", "P4_Figure1_projections.pdf"), width = 10, height = 8)
print(fig_main)
dev.off()

png(here("results", "P4_Figure1_projections.png"), width = 10, height = 8, units = "in", res = 300)
print(fig_main)
dev.off()

figS1 <- figS1a + figS1b + plot_layout(ncol = 2)

pdf(here("results", "P4_FigureS1_current_vs_future_maps.pdf"), width = 10, height = 5)
print(figS1)
dev.off()

png(here("results", "P4_FigureS1_current_vs_future_maps.png"), width = 10, height = 5, units = "in", res = 300)
print(figS1)
dev.off()

cat("\n=== DONE ===\n")