library(here)
library(sf)
library(ggplot2)
library(spdep)
library(patchwork)
library(splines)

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

# (AZMET + NWS Tucson)
dir.create(here("data", "weather"), recursive = TRUE, showWarnings = FALSE)
wx_path <- here("data", "weather", "tucson_hourly.csv")

if (!file.exists(wx_path)) {
  
  # NWS Tucson International Airport (KTUS)
  years <- 2018:2026
  wx_list <- list()
  
  for (yr in years) {
    url <- sprintf(
      "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=TUS&data=tmpf&data=relh&tz=America/Phoenix&format=onlycomma&latlon=no&elev=no&missing=M&trace=T&direct=no&report_type=3&year1=%d&month1=5&day1=1&year2=%d&month2=10&day2=1",
      yr, yr
    )
    tmp <- tempfile(fileext = ".csv")
    result <- try(download.file(url, tmp, mode = "w", quiet = TRUE), silent = TRUE)
    
    if (!inherits(result, "try-error") && file.size(tmp) > 500) {
      d <- read.csv(tmp, stringsAsFactors = FALSE)
      if (nrow(d) > 0 && "tmpf" %in% colnames(d)) {
        wx_list[[length(wx_list) + 1]] <- d
        cat("  IEM", yr, "downloaded:", nrow(d), "rows\n")
      }
    }
  }
  
  if (length(wx_list) > 0) {
    wx <- do.call(rbind, wx_list)
    wx$tmpf <- as.numeric(wx$tmpf)
    wx$relh <- as.numeric(wx$relh)
    wx$datetime <- as.POSIXct(wx$valid, format = "%Y-%m-%d %H:%M", tz = "America/Phoenix")
    wx$date <- as.Date(wx$datetime, tz = "America/Phoenix")
    wx$hour <- as.integer(format(wx$datetime, "%H"))
    wx$year <- as.integer(format(wx$datetime, "%Y"))
    wx$month <- as.integer(format(wx$datetime, "%m"))
    wx$doy <- as.integer(format(wx$datetime, "%j"))
    wx <- wx[!is.na(wx$tmpf) & !is.na(wx$relh), ]
    wx$temp_c <- (wx$tmpf - 32) * 5 / 9
    write.csv(wx, wx_path, row.names = FALSE)
  } else {
    stop("Weather download failed")
  }
}

wx <- read.csv(wx_path, stringsAsFactors = FALSE)
wx$datetime <- as.POSIXct(wx$datetime, tz = "America/Phoenix")
wx$date <- as.Date(wx$date)

# Cooler failure days
# Stull (2011) wet-bulb approximation
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

# Compute failure for each eta
for (eta_val in eta_values) {
  col <- paste0("failure_eta", eta_val * 10)
  wx[[col]] <- calc_supply(wx$temp_c, wx$relh, eta = eta_val) > 27
}

# Primary analysis uses eta = 0.65
wx$failure <- calc_supply(wx$temp_c, wx$relh, eta = 0.65) > 27  # TRUE = failure

daily <- aggregate(
  cbind(failure_hours = failure, max_temp = tmpf, max_rh = relh) ~ date + year + month + doy,
  data = wx,
  FUN = function(x) c(sum(x), max(x), max(x))[1]
)

daily_fail <- aggregate(failure ~ date + year + month + doy, data = wx, FUN = sum)
colnames(daily_fail)[5] <- "failure_hours"

daily_tmax <- aggregate(tmpf ~ date, data = wx, FUN = max, na.rm = TRUE)
colnames(daily_tmax)[2] <- "tmax_f"

daily_rhmax <- aggregate(relh ~ date, data = wx, FUN = max, na.rm = TRUE)
colnames(daily_rhmax)[2] <- "rh_max"

daily <- merge(daily_fail, daily_tmax, by = "date")
daily <- merge(daily, daily_rhmax, by = "date")

daily$failure_day <- as.integer(daily$failure_hours > 0)

# Assign season: pre-monsoon (May-Jun), monsoon (Jul-Sep), post-monsoon (Oct)
daily$season <- "Pre-monsoon"
daily$season[daily$month >= 7 & daily$month <= 9] <- "Monsoon"
daily$season[daily$month == 10] <- "Post-monsoon"
daily$season <- factor(daily$season, levels = c("Pre-monsoon", "Monsoon", "Post-monsoon"))

annual_season <- aggregate(
  cbind(failure_days = failure_day, total_failure_hours = failure_hours) ~ year + season,
  data = daily, FUN = sum
)

streak_fun <- function(x) {
  if (sum(x) == 0) return(0)
  r <- rle(x)
  max(r$lengths[r$values == 1])
}

streaks <- aggregate(failure_day ~ year + season, data = daily, FUN = streak_fun)
colnames(streaks)[3] <- "longest_streak"

annual_season <- merge(annual_season, streaks, by = c("year", "season"))

clim <- aggregate(
  cbind(failure_days, total_failure_hours, longest_streak) ~ season,
  data = annual_season, FUN = mean
)
clim_se <- aggregate(
  cbind(failure_days, total_failure_hours, longest_streak) ~ season,
  data = annual_season,
  FUN = function(x) sd(x) / sqrt(length(x))
)
colnames(clim_se)[2:4] <- paste0(colnames(clim_se)[2:4], "_se")

table1 <- merge(clim, clim_se, by = "season")
table1[, 2:7] <- round(table1[, 2:7], 1)

write.csv(table1, here("results", "P3_Table1_failure_climatology.csv"), row.names = FALSE)

# Sensitivity table
eta_sensitivity <- lapply(eta_values, function(eta_val) {
  col <- paste0("failure_eta", eta_val * 10)
  daily_eta <- aggregate(wx[[col]] ~ date, data = wx, FUN = sum)
  colnames(daily_eta) <- c("date", "failure_hours")
  daily_eta$year <- as.integer(format(as.Date(daily_eta$date), "%Y"))
  daily_eta$failure_day <- as.integer(daily_eta$failure_hours > 0)
  annual <- aggregate(failure_day ~ year, data = daily_eta, FUN = sum)
  data.frame(eta = eta_val,
             mean_failure_days = round(mean(annual$failure_day), 1),
             sd = round(sd(annual$failure_day), 1),
             min = min(annual$failure_day),
             max = max(annual$failure_day))
})
eta_sensitivity <- do.call(rbind, eta_sensitivity)
write.csv(eta_sensitivity, here("results", "P3_Table_eta_sensitivity.csv"), row.names = FALSE)
write.csv(annual_season, here("results", "P3_TableS1_annual_season_breakdown.csv"), row.names = FALSE)

# Compound exposure
mean_failure_days <- mean(aggregate(failure_day ~ year, data = daily, FUN = sum)$failure_day)
bg$compound_exposure <- bg$evap_prop * mean_failure_days
bg$exposure_q <- cut(bg$compound_exposure,
                     breaks = quantile(bg$compound_exposure, probs = 0:4/4),
                     labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                     include.lowest = TRUE)

high_exp <- bg[bg$exposure_q == "Q4 (highest)", ]
low_exp  <- bg[bg$exposure_q == "Q1 (lowest)", ]

compare_vars <- c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built")
compare_labels <- c("Evap. prevalence", "Median income ($)", "Minority (%)",
                    "Renter (%)", "Year built")

table2 <- data.frame(
  variable = compare_labels,
  high_mean = sapply(compare_vars, function(v) round(mean(high_exp[[v]], na.rm = TRUE), 2)),
  high_se = sapply(compare_vars, function(v) round(sd(high_exp[[v]], na.rm = TRUE) / sqrt(nrow(high_exp)), 2)),
  low_mean = sapply(compare_vars, function(v) round(mean(low_exp[[v]], na.rm = TRUE), 2)),
  low_se = sapply(compare_vars, function(v) round(sd(low_exp[[v]], na.rm = TRUE) / sqrt(nrow(low_exp)), 2)),
  wilcox_p = sapply(compare_vars, function(v) {
    signif(wilcox.test(high_exp[[v]], low_exp[[v]])$p.value, 3)
  })
)
write.csv(table2, here("results", "P3_Table2_demographics_by_exposure.csv"), row.names = FALSE)

# Spatial weights
coords <- cbind(bg$lon, bg$lat)
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")
moran_exp <- moran.test(bg$compound_exposure, lw)

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)
heatmap_data <- aggregate(
  cbind(mean_temp = temp_c, mean_rh = relh, fail_prop = failure) ~ doy + hour,
  data = wx, FUN = mean, na.rm = TRUE
)

fig1a <- ggplot(heatmap_data, aes(x = doy, y = hour, fill = mean_temp)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Temp (°C)") +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Day of year", y = "Hour", title = "a")

fig1b <- ggplot(heatmap_data, aes(x = doy, y = hour, fill = mean_rh)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#f7f7f7", "#74add1", "#313695"),
                       name = "RH (%)") +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Day of year", y = "Hour", title = "b")

fig1c <- ggplot(heatmap_data, aes(x = doy, y = hour, fill = fail_prop)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#f7f7f7", "#fc8d59", "#b30000"),
                       name = "Failure\nprobability",
                       limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 23, 4)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Day of year", y = "Hour", title = "c")

fig1 <- fig1a / fig1b / fig1c
pdf(here("results", "P3_Figure1_heatmap.pdf"), width = 8, height = 10)
print(fig1)
dev.off()

png(here("results", "P3_Figure1_heatmap.png"), width = 8, height = 10, units = "in", res = 300)
print(fig1)
dev.off()

bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = compound_exposure), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Compound\nexposure") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P3_Figure2_compound_exposure_map.pdf"), width = 7, height = 7)
print(fig2)
dev.off()

png(here("results", "P3_Figure2_compound_exposure_map.png"), width = 7, height = 7, units = "in", res = 300)
print(fig2)
dev.off()

wx$failure_065 <- calc_supply(wx$temp_c, wx$relh, eta = 0.65) > 27
daily_065 <- aggregate(failure_065 ~ date, data = wx, FUN = sum)
colnames(daily_065) <- c("date", "failure_hours")
daily_065$year <- as.integer(format(as.Date(daily_065$date), "%Y"))
daily_065$failure_day <- as.integer(daily_065$failure_hours > 0)
annual_065 <- aggregate(failure_day ~ year, data = daily_065, FUN = sum)

figS1 <- ggplot(annual_065, aes(x = year, y = failure_day)) +
  geom_col(fill = "#d73027", alpha = 0.7, width = 0.6) +
  geom_hline(yintercept = mean(annual_065$failure_day), linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Year", y = "Total failure days (May-Sep)", title = "")

pdf(here("results", "P3_FigureS1_annual_failure_days.pdf"), width = 6, height = 4)
print(figS1)
dev.off()

png(here("results", "P3_FigureS1_annual_failure_days.png"), width = 6, height = 4, units = "in", res = 300)
print(figS1)
dev.off()

figS2 <- ggplot(bg, aes(x = evap_prop, y = compound_exposure)) +
  geom_point(size = 1.2, alpha = 0.5, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "#d73027",
              fill = "#fc8d59", alpha = 0.2, linewidth = 0.7) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Evap. cooler prevalence",
       y = "Compound exposure (prevalence x failure days)",
       title = "")

pdf(here("results", "P3_FigureS2_scatter_compound.pdf"), width = 5, height = 5)
print(figS2)
dev.off()

png(here("results", "P3_FigureS2_scatter_compound.png"), width = 5, height = 5, units = "in", res = 300)
print(figS2)
dev.off()

eta_heatmaps <- lapply(eta_values, function(eta_val) {
  col <- paste0("failure_eta", eta_val * 10)
  hd_eta <- aggregate(wx[[col]] ~ doy + hour, data = wx, FUN = mean, na.rm = TRUE)
  colnames(hd_eta) <- c("doy", "hour", "fail_prop")
  ggplot(hd_eta, aes(x = doy, y = hour, fill = fail_prop)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("#f7f7f7", "#fc8d59", "#b30000"),
                         name = "Failure\nprob.", limits = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 23, 4)) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.2, "cm"),
          plot.title = element_text(size = 9, face = "bold")) +
    labs(x = "Day of year", y = "Hour", title = bquote(eta == .(eta_val)))
})

figS3 <- (eta_heatmaps[[1]] + eta_heatmaps[[2]] + eta_heatmaps[[3]]) /
  (eta_heatmaps[[4]] + plot_spacer())

pdf(here("results", "P3_FigureS3_eta_heatmaps.pdf"), width = 12, height = 7)
print(figS3)
dev.off()

png(here("results", "P3_FigureS3_eta_heatmaps.png"), width = 12, height = 7, units = "in", res = 300)
print(figS3)
dev.off()

figS4 <- ggplot(eta_sensitivity, aes(x = factor(eta), y = mean_failure_days)) +
  geom_col(fill = "#c47a4a", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = pmax(mean_failure_days - sd, 0),
                    ymax = mean_failure_days + sd),
                width = 0.2, color = "grey40") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Saturation efficiency (eta)", y = "Mean failure days per summer", title = "")

pdf(here("results", "P3_FigureS4_eta_sensitivity.pdf"), width = 6, height = 4)
print(figS4)
dev.off()

png(here("results", "P3_FigureS4_eta_sensitivity.png"), width = 6, height = 4, units = "in", res = 300)
print(figS4)
dev.off()

