library(here)
library(sf)
library(ggplot2)
library(spdep)
library(sf)
library(spatialreg)
library(patchwork)
library(maps)
library(cowplot)

bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))

# Drop empty polygons
bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]

# Keep only those inside urban areas
UrbanAreas <- st_read(here("data", "2020_Arizona_Census_Urban_Areas", "2020_Arizona_Census_Urban_Areas.shp"))
UrbanAreas <- st_transform(UrbanAreas, st_crs(bg_sf))
UrbanAreas <- st_make_valid(UrbanAreas)
bg_centroids <- st_centroid(bg_sf)
inside <- st_intersects(bg_centroids, st_union(UrbanAreas), sparse = FALSE)[, 1]
bg_sf <- bg_sf[inside, ]
bg <- data.frame(bg_sf)

# Rename columns
colnames(bg)[colnames(bg) == "geoid20"]   <- "GEOID"
colnames(bg)[colnames(bg) == "evp_prp"]   <- "evap_prop"
colnames(bg)[colnames(bg) == "med_inc"]    <- "med_income"
colnames(bg)[colnames(bg) == "pct_mnr"]    <- "pct_minority"
colnames(bg)[colnames(bg) == "ave_age"]    <- "med_year_built"
colnames(bg)[colnames(bg) == "pct_rnt"]    <- "pct_renter"
colnames(bg)[colnames(bg) == "pct_sfr"]    <- "pct_sfh"
colnames(bg)[colnames(bg) == "covennt"]    <- "covenant"

# Compute centroids for lon/lat
bg_sf <- st_as_sf(bg)
bg_sf <- st_transform(bg_sf, 4326)
coords <- st_coordinates(st_centroid(bg_sf))
bg$lon <- coords[, 1]
bg$lat <- coords[, 2]
bg_sf <- st_make_valid(bg_sf)

# Create quartiles
bg$evap_q <- cut(bg$evap_prop,
                 breaks = quantile(bg$evap_prop, probs = 0:4/4),
                 labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                 include.lowest = TRUE)

summarize_by_quartile <- function(x, group) {
  means <- tapply(x, group, mean, na.rm = TRUE)
  sds <- tapply(x, group, sd, na.rm = TRUE)
  ns <- tapply(x, group, length)
  ses <- sds / sqrt(ns)
  data.frame(quartile = names(means), mean = round(means, 2),
             sd = round(sds, 2), se = round(ses, 2), n = as.integer(ns))
}

vars <- c("med_income", "pct_minority", "med_year_built", "pct_renter", "pct_sfh")
var_labels <- c("Median income ($)", "Minority (%)", "Median year built",
                "Renter (%)", "Single-family (%)")

table1_list <- lapply(vars, function(v) {
  out <- summarize_by_quartile(bg[[v]], bg$evap_q)
  out$variable <- v
  out
})
table1 <- do.call(rbind, table1_list)
table1$var_label <- rep(var_labels, each = 4)

# Kruskal-Wallis test across quartiles for each variable
kw_tests <- sapply(vars, function(v) {
  kruskal.test(bg[[v]] ~ bg$evap_q)$p.value
})
names(kw_tests) <- var_labels

write.csv(table1, here("results", "P1_Table1_quartile_summary.csv"), row.names = FALSE)

# Spearman correlations
cor_vars <- c("med_income", "pct_minority", "med_year_built", "pct_renter")
cor_labels <- c("Median income", "% Minority", "Median year built", "% Renter")

cor_results <- data.frame(
  variable = cor_labels,
  rho = sapply(cor_vars, function(v) cor(bg$evap_prop, bg[[v]], method = "spearman")),
  p = sapply(cor_vars, function(v) cor.test(bg$evap_prop, bg[[v]], method = "spearman")$p.value)
)
cor_results$rho <- round(cor_results$rho, 3)
cor_results$p <- signif(cor_results$p, 3)

# Proportion of evap-cooler households in vs out of covenanted areas
cov_in <- bg$evap_prop[bg$covenant == 1]
cov_out <- bg$evap_prop[bg$covenant == 0]

# Wilcoxon test
wt <- wilcox.test(cov_in, cov_out)
cat("Wilcoxon p-value:", signif(wt$p.value, 3), "\n")

# Odds ratio
evap_high <- as.integer(bg$evap_prop > median(bg$evap_prop))
cov_table <- table(evap_high, bg$covenant)
or_fish <- fisher.test(cov_table)

# Predictors of Evap Cooler prevalence
bg$z_income <- scale(bg$med_income)
bg$z_minority <- scale(bg$pct_minority)
bg$z_year_built <- scale(bg$med_year_built)
bg$z_renter <- scale(bg$pct_renter)

m1 <- glm(evap_prop ~ z_income + z_minority + z_year_built + z_renter + covenant,
          family = quasibinomial, data = bg)

coef_table <- data.frame(
  variable = c("Intercept", "Income (z)", "Minority (z)", "Year built (z)",
               "Renter (z)", "Covenant"),
  estimate = round(coef(m1), 3),
  se = round(summary(m1)$coefficients[, 2], 3),
  p = signif(summary(m1)$coefficients[, 4], 3)
)
write.csv(coef_table, here("results", "P1_Table2_GLM_results.csv"), row.names = FALSE)

# Spatial autocorrelation of residuals
coords <- cbind(bg$lon, bg$lat)
nb <- knn2nb(knearneigh(coords, k = 5))
lw <- nb2listw(nb, style = "W")
moran_resid <- moran.test(residuals(m1), lw)

if (moran_resid$p.value < 0.05) {
  m_spatial <- errorsarlm(evap_prop ~ z_income + z_minority + z_year_built + z_renter + covenant,
                          data = bg, listw = lw)
  print(summary(m_spatial))
  
  se_coefs <- summary(m_spatial)$Coef
  coef_table_spatial <- data.frame(
    variable = rownames(se_coefs),
    estimate = round(se_coefs[, 1], 3),
    se = round(se_coefs[, 2], 3),
    p = signif(se_coefs[, 4], 3)
  )
  write.csv(coef_table_spatial, here("results", "P1_Table2b_spatial_error_model.csv"), row.names = FALSE)
}

# LISA clustering
lisa <- localmoran(bg$evap_prop, lw)
bg$lisa_I <- lisa[, 1]
bg$lisa_p <- lisa[, 5]
bg$lag_evap <- lag.listw(lw, bg$evap_prop)
mean_evap <- mean(bg$evap_prop)
mean_lag <- mean(bg$lag_evap)

bg$lisa_cluster <- "Not significant"
sig <- bg$lisa_p < 0.05
bg$lisa_cluster[sig & bg$evap_prop > mean_evap & bg$lag_evap > mean_lag] <- "High-High"
bg$lisa_cluster[sig & bg$evap_prop < mean_evap & bg$lag_evap < mean_lag] <- "Low-Low"
bg$lisa_cluster[sig & bg$evap_prop > mean_evap & bg$lag_evap < mean_lag] <- "High-Low"
bg$lisa_cluster[sig & bg$evap_prop < mean_evap & bg$lag_evap > mean_lag] <- "Low-High"

# Plots
pal_lisa <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
              "High-Low" = "#fdae61", "Low-High" = "#abd9e9",
              "Not significant" = "grey90")
pal_cov <- c("0" = "grey70", "1" = "#7b3294")

us <- map_data("state")
az <- us[us$region == "arizona", ]
inset <- ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               fill = "grey90", color = "grey70", size = 0.15) +
  geom_polygon(data = az, aes(x = long, y = lat, group = group),
               fill = "#b5d4be", color = "grey40", size = 0.3) +
  geom_point(aes(x = -110.9747, y = 32.2226), size = 1.5, color = "#d73027") +
  annotate("text", x = -110.9747, y = 30.5, label = "Tucson",
           size = 2, color = "#2c2418", fontface = "bold") +
  coord_fixed(1.3) +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  labs(title = "b") +
  theme(plot.title = element_text(size = 10, face = "bold"))
city_boundary <- st_union(bg_sf)

fig1_main <- ggplot() +
  geom_sf(data = bg_sf, aes(fill = evap_prop), color = "white", size = 0.15) +
  geom_sf(data = city_boundary, fill = NA, color = "grey30", size = 0.4, linetype = "solid") +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Evap. cooler\nprevalence") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a")

fig1 <- ggdraw() +
  draw_plot(fig1_main, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(inset, x = 0.62, y = 0.72, width = 0.35, height = 0.3)


bg_sf$lisa_cluster <- bg$lisa_cluster
fig2 <- ggplot(bg_sf) +
  geom_sf(aes(fill = lisa_cluster), color = "white", size = 0.15) +
  geom_sf(data = city_boundary, fill = NA, color = "grey30", size = 0.4, linetype = "solid") +
  scale_fill_manual(values = pal_lisa, name = "LISA cluster") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "a")

fig2 <- ggdraw() +
  draw_plot(fig2, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(inset, x = 0.62, y = 0.72, width = 0.35, height = 0.3)

make_scatter <- function(xvar, xlab, panel_label) {
  rho <- cor(bg$evap_prop, bg[[xvar]], method = "spearman")
  pval <- cor.test(bg$evap_prop, bg[[xvar]], method = "spearman")$p.value
  ptext <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", signif(pval, 2)))
  label <- paste0("rho == ", round(rho, 2))
  
  ggplot(bg, aes_string(x = xvar, y = "evap_prop")) +
    geom_point(size = 1.2, alpha = 0.5, color = "grey30") +
    geom_smooth(method = "loess", se = TRUE, color = "#cb181d",
                fill = "#fb6a4a", alpha = 0.2, size = 0.7) +
    annotate("text", x = Inf, y = Inf, label = label, parse = TRUE,
             hjust = 1.1, vjust = 1.5, size = 3) +
    annotate("text", x = Inf, y = Inf, label = ptext,
             hjust = 1.1, vjust = 3, size = 2.5, color = "grey40") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(x = xlab, y = "Evap. cooler prevalence", title = panel_label)
}

fig3a <- make_scatter("med_income", "Median household income ($)", "c")
fig3b <- make_scatter("pct_minority", "Minority proportion", "d")
fig3c <- make_scatter("med_year_built", "Median year built", "e")
fig3d <- make_scatter("pct_renter", "Renter proportion", "f")

or_dat <- data.frame(
  label = "High evap.\nin covenanted area",
  or = or_fish$estimate,
  lo = or_fish$conf.int[1],
  hi = or_fish$conf.int[2]
)

fig4 <- ggplot(or_dat, aes(x = or, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, size = 0.6) +
  geom_point(size = 3, color = "#7b3294") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(x = "Odds ratio (95% CI)", y = "")


pdf(here("results", "P1_Figure1.pdf"), width = 4, height = 5)
print(fig1)
dev.off()

pdf(here("results", "P1_Figure2.pdf"), width = 4, height = 5)
print(fig2)
dev.off()

fig1_combined <- fig1 + fig2 + plot_layout(ncol = 2)

pdf(here("results", "P1_Figure1_maps.pdf"), width = 10, height = 5)
print(fig1_combined)
dev.off()

png(here("results", "P1_Figure1_maps.png"), width = 10, height = 5, units = "in", res = 300)
print(fig1_combined)
dev.off()

fig3a <- make_scatter("med_income", "Median household income ($)", "a")
fig3b <- make_scatter("pct_minority", "Minority proportion", "b")
fig3c <- make_scatter("med_year_built", "Median year built", "c")
fig3d <- make_scatter("pct_renter", "Renter proportion", "d")

fig2_combined <- (fig3a + fig3b) / (fig3c + fig3d)

pdf(here("results", "P1_Figure2_scatterplots.pdf"), width = 8, height = 7)
print(fig2_combined)
dev.off()

png(here("results", "P1_Figure2_scatterplots.png"), width = 8, height = 7, units = "in", res = 300)
print(fig2_combined)
dev.off()

pdf(here("results", "P1_Figure3_covenant_OR.pdf"), width = 5, height = 3)
print(fig4 + labs(title = ""))
dev.off()

png(here("results", "P1_Figure3_covenant_OR.png"), width = 5, height = 3, units = "in", res = 300)
print(fig4 + labs(title = ""))
dev.off()

