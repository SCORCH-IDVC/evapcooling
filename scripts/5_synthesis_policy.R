library(here)
library(sf)
library(ggplot2)
library(spdep)
library(patchwork)

bg_sf <- st_read(here("data", "Q1 Data Shapefile", "pima_Q1_data.shp"))
bg_sf <- bg_sf[bg_sf$med_inc != 0 & !is.na(bg_sf$med_inc), ]
bg_sf <- bg_sf[!is.na(bg_sf$evp_prp), ]

UrbanAreas <- st_read(here("data", "2020_Arizona_Census_Urban_Areas",
                            "2020_Arizona_Census_Urban_Areas.shp"))
UrbanAreas <- st_transform(UrbanAreas, st_crs(bg_sf))
UrbanAreas <- st_make_valid(UrbanAreas)
bg_centroids <- st_centroid(bg_sf)
inside <- st_intersects(bg_centroids, st_union(UrbanAreas), sparse = FALSE)[, 1]
bg_sf <- bg_sf[inside, ]

bg <- data.frame(bg_sf)
colnames(bg)[colnames(bg) == "geoid20"]  <- "GEOID"
colnames(bg)[colnames(bg) == "evp_prp"]  <- "evap_prop"
colnames(bg)[colnames(bg) == "med_inc"]  <- "med_income"
colnames(bg)[colnames(bg) == "pct_mnr"]  <- "pct_minority"
colnames(bg)[colnames(bg) == "ave_age"]  <- "med_year_built"
colnames(bg)[colnames(bg) == "pct_rnt"]  <- "pct_renter"
colnames(bg)[colnames(bg) == "pct_sfr"]  <- "pct_sfh"
colnames(bg)[colnames(bg) == "covennt"]  <- "covenant"

bg_sf <- st_as_sf(bg)
bg_sf <- st_transform(bg_sf, 4326)
coords <- st_coordinates(st_centroid(bg_sf))
bg$lon <- coords[, 1]
bg$lat <- coords[, 2]
bg_sf <- st_make_valid(bg_sf)

cat("Block groups loaded:", nrow(bg), "\n")

# Cooling access deserts
# A cooling access desert is a block group where:
#   - evap prevalence > 75th percentile
#   - median income < 25th percentile

q75_evap <- quantile(bg$evap_prop, 0.75)
q25_income <- quantile(bg$med_income, 0.25)
bg$cooling_desert <- bg$evap_prop >= q75_evap & bg$med_income <= q25_income
desert <- bg[bg$cooling_desert, ]
non_desert <- bg[!bg$cooling_desert, ]

desert_compare <- data.frame(
  variable = c("Evap. prevalence", "Median income ($)", "Minority (%)",
               "Renter (%)", "Year built"),
  desert_mean = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built"),
                        function(v) round(mean(desert[[v]], na.rm = TRUE), 2)),
  city_mean = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built"),
                      function(v) round(mean(bg[[v]], na.rm = TRUE), 2)),
  wilcox_p = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter", "med_year_built"),
                     function(v) signif(wilcox.test(desert[[v]], non_desert[[v]])$p.value, 3))
)

print(desert_compare)
write.csv(desert_compare, here("results", "P7_Table1_cooling_deserts.csv"), row.names = FALSE)

# Composite vulnerability
bg$z_evap <- as.numeric(scale(bg$evap_prop))
bg$z_income_inv <- as.numeric(scale(-bg$med_income))  # invert: lower = worse
bg$z_minority <- as.numeric(scale(bg$pct_minority))
bg$z_renter <- as.numeric(scale(bg$pct_renter))
bg$z_age_inv <- as.numeric(scale(-bg$med_year_built))  # invert: older = worse

## Equal-weight composite
bg$vulnerability_index <- (bg$z_evap + bg$z_income_inv + bg$z_minority +
                            bg$z_renter + bg$z_age_inv) / 5

## Classify
bg$vuln_q <- cut(bg$vulnerability_index,
                  breaks = quantile(bg$vulnerability_index, probs = 0:4/4),
                  labels = c("Low", "Moderate", "High", "Very high"),
                  include.lowest = TRUE)

cat("\n=== Q3: VULNERABILITY INDEX ===\n")
cat("Range:", round(range(bg$vulnerability_index), 2), "\n")
print(table(bg$vuln_q))

## Demographics by vulnerability quartile
vuln_summary <- data.frame(
  quartile = levels(bg$vuln_q),
  evap = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$evap_prop[bg$vuln_q == q]) * 100, 1)),
  income = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$med_income[bg$vuln_q == q]))),
  minority = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$pct_minority[bg$vuln_q == q]) * 100, 0)),
  renter = sapply(levels(bg$vuln_q), function(q)
    round(mean(bg$pct_renter[bg$vuln_q == q]) * 100, 0))
)
print(vuln_summary)
write.csv(vuln_summary, here("results", "P7_Table3_vulnerability_index.csv"), row.names = FALSE)

# Renter-owner
renter_q <- cut(bg$pct_renter,
                 breaks = quantile(bg$pct_renter, probs = 0:4/4),
                 labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
                 include.lowest = TRUE)

renter_evap <- tapply(bg$evap_prop, renter_q, mean, na.rm = TRUE)
renter_income <- tapply(bg$med_income, renter_q, mean, na.rm = TRUE)

cat("\n=== Q4: RENTER-OWNER DIVIDE ===\n")
cat("Mean evap prevalence by renter quartile:\n")
print(round(renter_evap * 100, 1))
cat("Kruskal-Wallis p:", signif(kruskal.test(bg$evap_prop ~ renter_q)$p.value, 3), "\n")

# Covenant persistence
cov <- bg[bg$covenant == 1, ]
non_cov <- bg[bg$covenant == 0, ]

persistence <- data.frame(
  dimension = c("Evap. prevalence", "Income", "Minority %", "Renter %",
                "Year built", "Vulnerability index"),
  covenanted = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter",
                         "med_year_built", "vulnerability_index"),
                       function(v) round(mean(cov[[v]], na.rm = TRUE), 3)),
  non_covenanted = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter",
                              "med_year_built", "vulnerability_index"),
                            function(v) round(mean(non_cov[[v]], na.rm = TRUE), 3)),
  wilcox_p = sapply(c("evap_prop", "med_income", "pct_minority", "pct_renter",
                       "med_year_built", "vulnerability_index"),
                     function(v) signif(wilcox.test(cov[[v]], non_cov[[v]])$p.value, 3))
)


# Figures
bg_sf2 <- st_as_sf(bg)
if (is.na(st_crs(bg_sf2))) bg_sf2 <- st_set_crs(bg_sf2, 4326)
bg_sf2 <- st_make_valid(bg_sf2)

bg_sf2$cooling_desert <- bg$cooling_desert
fig1 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = factor(cooling_desert, levels = c(FALSE, TRUE))),
          color = "white", size = 0.15) +
  scale_fill_manual(values = c("FALSE" = "#e8e8e8", "TRUE" = "#d73027"),
                    labels = c("No", "Yes"), name = "Cooling\ndesert") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P7_Figure1_cooling_deserts.pdf"), width = 7, height = 7)
print(fig1)
dev.off()

png(here("results", "P7_Figure1_cooling_deserts.png"), width = 7, height = 7, units = "in", res = 300)
print(fig1)
dev.off()

bg_sf2$vulnerability_index <- bg$vulnerability_index

fig2 <- ggplot(bg_sf2) +
  geom_sf(aes(fill = vulnerability_index), color = "white", size = 0.15) +
  scale_fill_gradientn(colors = c("#2c7bb6", "#abd9e9", "#fee090", "#d73027"),
                       name = "Vulnerability\nindex") +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.25),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.3, "cm"),
        plot.title = element_text(size = 10, face = "bold")) +
  labs(title = "")

pdf(here("results", "P7_Figure2_vulnerability_index.pdf"), width = 7, height = 7)
print(fig2)
dev.off()

png(here("results", "P7_Figure2_vulnerability_index.png"), width = 7, height = 7, units = "in", res = 300)
print(fig2)
dev.off()