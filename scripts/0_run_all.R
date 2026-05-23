library(here)

dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("docs"), recursive = TRUE, showWarnings = FALSE)

source(here("scripts", "1_equity_environmental_justice.R"), local = (e1 <- new.env()))
save(list = ls(e1), envir = e1, file = here("docs", "ws_1_equity.RData"))

source(here("scripts", "2_heat_island.R"), local = (e2 <- new.env()))
save(list = ls(e2), envir = e2, file = here("docs", "ws_2_heat.RData"))
source(here("scripts", "3_monsoon_heat.R"), local = (e3 <- new.env()))
save(list = ls(e3), envir = e3, file = here("docs", "ws_3_monsoon.RData"))
source(here("scripts", "4_projections.R"), local = (e4 <- new.env()))
save(list = ls(e4), envir = e4, file = here("docs", "ws_4_projections.RData"))
source(here("scripts", "5_health_outcomes.R"), local = (e5 <- new.env()))
save(list = ls(e5), envir = e5, file = here("docs", "ws_5_health.RData"))
source(here("scripts", "6_projected_health.R"), local = (e6 <- new.env()))
save(list = ls(e6), envir = e6, file = here("docs", "ws_6_projected_health.RData"))
source(here("scripts", "7_synthesis_policy.R"), local = (e7 <- new.env()))
save(list = ls(e7), envir = e7, file = here("docs", "ws_7_synthesis_policy.RData"))

quarto::quarto_render(here("docs", "results.qmd"), output_format = "html", output_file = "index.html")
