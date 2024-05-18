# Mike Vlah
# vlahm13@gmail.com
# last data retrieval dates given in comments below:
# last edit: 2025-05-13

library(dataRetrieval)
library(data.table)
library(lubridate)
library(glue)
library(dygraphs)
library(xts)
library(hydroGOF)
library(ggplot2)
library(gridExtra)
# library(neonUtilities)
library(htmlwidgets)
# library(DAAG)
library(glmnet)
# library(segmented)
library(tidyverse)
# library(macrosheds)
library(parallel)
library(doParallel)
library(yaml)

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

source('src/00_helpers.R')
donor_gauges <- read_yaml('cfg/donor_gauges.yml')
sites <- read_csv('cfg/sites.csv')

## 1. setup ####

build_dir_structure()

results <- tibble(
    site_code = names(donor_gauges), method = NA, nse = NA, nse_cv = NA, kge = NA,
    kge_cv = NA, pbias = NA, pbias_cv = NA, bestmod = NA, adj_r_squared = NA
)

## 2. ridge (glmnet) regression on specific discharge ####

#WARNING: if you're on a unix-like machine, fork clusters will be created.
#these can be unstable when Rstudio is involved, so watch for memory leaks.
#i encountered a few by running glmnet models piecemeal. if you run this whole
#script at once, rather than going through line by line, you should be okay.
#you can also grep 00_helpers.R for clst_type and hard-code it to "PSOCK"

donor_gauges <- read_yaml('cfg/donor_gauges.yml')

regress(site_code = 'erwin', framework = 'glmnet',
        # custom_formula = 'discharge_log ~ `{paste(paste0(gage_ids, "_log"), collapse = "`*`")}`',
        seasons = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1),
        # seasons = c(1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 1, 1),
        bootstrap_ci = FALSE,
        ncores = 4, target_daterange = c('2016-01-01', '2024-04-30'))

## other scenarios (not yet adapted for this workflow): segmented regression ####

#SYCA: erroneous outlier (sites separated by < 10 km on same river with no major intervening
#   tribs. very unlikely for one to report >100L/s while the other reports ~15L/s).
#   The questionable measurement also followed a long sampling hiatus.
syca_d <- assemble_q_df(neon_site = 'SYCA', nearby_usgs_gages = donor_gauges[['SYCA']]) %>%
    filter(! as.Date(datetime) == as.Date('2021-07-26'))
regress(neon_site = 'SYCA', framework = 'segmented', precomputed_df = syca_d)

## other scenarios (not yet adapted for this workflow): intercept shift ####

#LECO: a wildfire in 2016 changed the intercept!
leco_d <- assemble_q_df(neon_site = 'LECO', nearby_usgs_gages = donor_gauges[['LECO']]) %>%
    mutate(dummy = ifelse(datetime > as.POSIXct('2016-11-27'), 1, 0),
           dummy = as.factor(dummy))
# plot(leco_d$`03497300_log`, leco_d$discharge_log, col = leco_d$dummy)
# plot(leco_d$datetime, leco_d$discharge_log, type = 'b', col = leco_d$dummy, pch = 20)
regress(neon_site = 'LECO', framework = 'lm', precomputed_df = leco_d,
        custom_formula = 'discharge_log ~ `{paste0(donor_gauges[["LECO"]], "_log")}` * season + dummy',
        dummy_break = as.Date('2016-11-27'))

## other scenarios (not yet adapted for this workflow): non-USGS Q source ####

#MCRA: uses Andrews Experimental Forest data
mcra_d_ <- read_csv('in/hjandrews_q.txt') %>%
    filter(SITECODE %in% donor_gauges[['MCRA']]) %>%
    select(datetime = DATE_TIME, SITECODE, INST_Q) %>%
    mutate(INST_Q = INST_Q * 28.317) %>% #cfs to L/s
    pivot_wider(names_from = SITECODE, values_from = INST_Q) %>%
    arrange(datetime)
mcra_d <- assemble_q_df(neon_site = 'MCRA', ms_Q_data = mcra_d_)
regress(neon_site = 'MCRA', framework = 'glmnet', precomputed_df = mcra_d, ncores = 9)

## other scenarios (not yet adapted for this workflow): composite regression ####

gauge_ids <- donor_gauges[[neon_site]][1:2]
como_d_ <- ms_q %>%
    filter(site_code %in% gauge_ids) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
como_d1 <- assemble_q_df_daily(neon_site = neon_site, ms_Q_data = como_d_) %>%
    select(-contains('MARTINELLI'))
como_best1 <- regress(
    neon_site = neon_site, framework = 'glmnet', precomputed_df = como_d1,
    custom_formula = "discharge_log ~ `{paste(paste0(gauge_ids[1:2], \"_log\"), collapse = \"`*`\")}` * season",
    custom_gauges = gauge_ids, no_write = TRUE
)

gauge_ids <- donor_gauges[[neon_site]]
como_d_ <- ms_q %>%
    filter(site_code %in% gauge_ids) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
como_d2 <- assemble_q_df_daily(neon_site = neon_site, ms_Q_data = como_d_)
como_best2 <- regress(
    neon_site = neon_site, framework = 'glmnet',
    precomputed_df = como_d2, no_write = TRUE
)

results <- plots_and_results_daily_composite(
    neon_site = neon_site,
    best1 = como_best1,
    best2 = como_best2, #best2 must be the case with more gauges included
    results = results,
    in_df1 = como_d1, in_df2 = como_d2, bootstrap_ci = TRUE, ncores = 16
)

results$method[results$site_code == 'COMO'] <- 'glmnet composite'
write_csv(results, 'out/lm_out/results_specificq.csv')


## 3. ridge regression on absolute discharge ####

rename_dir_structure()
build_dir_structure()

results_specificq <- results
results <- tibble(
    site_code = names(donor_gauges), method = NA, nse = NA, nse_cv = NA, kge = NA,
    kge_cv = NA, pbias = NA, pbias_cv = NA, bestmod = NA, adj_r_squared = NA
)

regress(site_code = 'erwin', framework = 'glmnet', scale_q_by_area = FALSE,
        target_daterange = c('2023-01-01', '2023-12-31'))

qq = read_csv('erwin_q_for_comparison.csv')
qq$NHC_discharge_Ls = qq$NHC_discharge_m3s * 1000
plot(qq$DateTime_UTC, qq$NHC_discharge_Ls, ylim = c(0, 5000))
qp = read_csv('out/lm_out/predictions/erwin.csv')
points(qp$datetime, qp$Q_predicted, col = 'red', pch = '.')
qm = read_csv('in/field_Q/erwin.csv')
points(qm$datetime, qm$discharge, col = 'blue', pch = 1)

plot(qp$datetime, qp$Q_predicted, col = 'red', pch = '.')

write_csv(results, 'out/lm_out/results.csv')
