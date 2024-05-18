#helper functions for the other R scripts. No need to source these directly.
#start with 01_linear_regression.R

# general (regression) ####

Mode <- function(x, na.rm = TRUE){

    if(na.rm){
        x <- na.omit(x)
    }

    ux <- unique(x)
    mode_out <- ux[which.max(tabulate(match(x, ux)))]
    return(mode_out)

}

neglog <- function(x){
    sgn <- sign(x)
    neglog_ <- log(abs(x) + 1) * sgn
    return(neglog_)
}

inv_neglog <- function(x){
    sgn <- sign(x)
    inv_neglog_ <- (exp(x) - 1) * sgn
    return(inv_neglog_)
}

# setup (regression) ####

build_dir_structure <- function(){

    dir.create('out', showWarnings = FALSE)
    dir.create('figs', showWarnings = FALSE)

    dir.create('in/usgs_Q', showWarnings = FALSE)

    dir.create('figs/lm_plots', showWarnings = FALSE)
    dir.create('figs/lm_plots/diag', showWarnings = FALSE)
    dir.create('figs/lm_plots/pred', showWarnings = FALSE)
    dir.create('figs/lm_plots/fit', showWarnings = FALSE)
    dir.create('figs/lm_plots/val', showWarnings = FALSE)

    dir.create('out/lm_out', showWarnings = FALSE)
    dir.create('out/lm_out/predictions', showWarnings = FALSE)
    dir.create('out/lm_out/fit', showWarnings = FALSE)
    dir.create('out/lm_out/summary', showWarnings = FALSE)
}

rename_dir_structure <- function(){

    file.rename('figs/lm_plots', 'figs/lm_plots_specQ')
    file.rename('out/lm_out', 'out/lm_out_specQ')
}

# data prep ####

assemble_q_df <- function(site_code, nearby_usgs_gages = NULL, ms_Q_data = NULL,
                          datetime_snapdist_hrs = 12,
                          seasons = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4),
                          target_daterange = NULL, overwrite = FALSE,
                          scale_q_by_area = TRUE){

    ms_Q_data <- NULL #

    #ms_Q_data: [not ready yet. hard-coded NULL] data.frame with posixct datetime column, and Q columns. Each Q column must
    #   have its MacroSheds site_code as header. transformed columns will be created
    #datetime_snapdist_hrs: numeric. when joining donor gauge data to field
    #   discharge by datetime, allow joins of up to this number of hrs earlier/later
    #seasons: a numeric vector of length 12, indicating how months should be grouped into seasons.
    #target_daterange: vector of two dates for the beginning and end of the
    #   period to be estimated. Donor gauges with insufficient data may be dropped before
    #   modeling, but those gauges might actually have data for your daterange
    #   of interest, so use this parameter to get the most out of the available donor gauges.

    #can't remember if i built this to handle nearby_usgs_gages and ms_Q_data ,
    #but in any case that wasn't relevant for this study.

    if(! paste0(site_code, '.csv') %in% list.files('in/field_Q')){
        stop('File in/field_Q/<site_code>.csv is missing.')
    }

    #load field discharge
    field_q = read_csv(glue('in/field_Q/{site_code}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE) %>%
        mutate(site_code := !!site_code)
    earliest_date = as.character(date(min(field_q$datetime)))

    #checks
    aa = lubridate::interval(target_daterange[1], target_daterange[2])
    bb = lubridate::interval(min(field_q$datetime), max(field_q$datetime))
    if(! int_overlaps(aa, bb)){
        stop('No overlap between target_daterange and field_Q samples')
    }

    if(! is.null(target_daterange)){

        if(length(target_daterange) != 2){
            stop('target_daterange must be a vector of 2')
        }

        tdr <- as.Date(target_daterange)
        if(tdr[2] <= tdr[1]){
            stop('the second date in target_daterange must come after the first')
        }
    } else {
        target_daterange <- c(earliest_date, '2099-01-01')
    }

    if(scale_q_by_area){
        wsa = filter(sites, site_code == !!site_code) %>% pull(ws_area_ha)
        field_q$discharge_manual = field_q$discharge_manual / wsa * 1000 #scale by 1000 so that neglog skew is negligible
    }

    # if(! file.exists(glue('in/usgs_Q/{site_code}.csv')) || overwrite){

    site_nearbyA = NULL
    # if(! is.null(nearby_usgs_gages)){

    usgsq = dataRetrieval::readNWISdata(
        sites = nearby_usgs_gages,
        service = 'iv', #instantaneous values
        parameterCd = '00060', #discharge (cfs)
        startDate = target_daterange[1],
        endDate = target_daterange[2]
    )

    if(! nrow(usgsq)) stop('no instantaneous Q?')

    # if(any(! usgsq$X_00060_00000_cd %in% c('A', 'P', 'A e', 'A R'))) stop()
    if(any(grepl('>', usgsq$X_00060_00000_cd))) warning(glue('">" detected in {site_code} donor Q'))
    if(any(grepl('<', usgsq$X_00060_00000_cd))) warning(glue('"<" detected in {site_code} donor Q'))
    if(any(usgsq$tz_cd != 'UTC')) stop('non-UTC datetime encountered')

    site_nearbyA = usgsq %>%
        as_tibble() %>%
        select(site_no, datetime = dateTime, discharge = X_00060_00000) %>%
        pivot_wider(names_from = site_no, values_from = discharge)

    site_nearbyA <- cull_gauges(site_nearbyA)

    if(scale_q_by_area){

        for(g in nearby_usgs_gages){

            if(g == '06190540'){
                #for NEON runthrough, this was the only site without
                #a watershed area included in dataRetrieval. unlikely
                #to be encountered again?
                site_nearbyA[[g]] = site_nearbyA[[g]] / 51089.73 * 1000
                next
            }

            nwissite = dataRetrieval::readNWISsite(g)

            if(is.null(nwissite)) stop('NWIS watershed area not listed for ', g)

            if(is.na(nwissite$contrib_drain_area_va)){
                wsa = nwissite$drain_area_va * 258.999 #mi^2 -> ha
            } else {
                wsa = nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
            }

            if(is.na(wsa)) stop('NWIS watershed area not listed for ', g)

            site_nearbyA[[g]] = site_nearbyA[[g]] / wsa * 1000
        }
    }

    site_nearbyA = site_nearbyA %>%
        mutate(across(matches('^[0-9]+$'), ~ . * 28.3168)) %>% #cfs -> L/s
        mutate(across(matches('^[0-9]+$'),
                      neglog,
                      # ~boxcox_write(., !!site_code, write = FALSE),
                      .names = '{.col}_log')) %>%
        arrange(datetime)
    # }

    if(! is.null(ms_Q_data)){

        ms_sites = grep('datetime', colnames(ms_Q_data), invert = TRUE, value = TRUE)

        if(scale_q_by_area){

            for(g in ms_sites){

                wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
                ms_Q_data[[g]] = ms_Q_data[[g]] / wsa * 1000
            }
        }

        ms_Q_data = ms_Q_data %>%
            mutate(across(-datetime, neglog, .names = '{.col}_log'))
    }

    if(! is.null(site_nearbyA) && ! is.null(ms_Q_data)){
        site_nearby = full_join(site_nearbyA, ms_Q_data, by = c('datetime'))
    } else if(! is.null(site_nearbyA)){
        site_nearby = site_nearbyA
    } else  {
        site_nearby = ms_Q_data
    }

    # site_nearby <- read_csv(glue('in/usgs_Q/{site_code}.csv'))
    # warning('OI!!')
    write_csv(site_nearby, glue('in/usgs_Q/{site_code}.csv'))

    #
    # } else {
    #    site_nearby = read_csv(glue('in/usgs_Q/{site_code}.csv'))
    # # }
    #
    #     warning('uncomment hard-code file read in assemble_q_df!!')

    #rolling join usgs data to field measurements by datetime
    x = rename(field_q, datetime_x = datetime) %>% as.data.table()
    y = rename(site_nearby, datetime_y = datetime) %>% as.data.table()

    rollmax = 60 * 60 * datetime_snapdist_hrs #join up to this number of hrs off
    x[, `:=` (datetime_min = datetime_x - rollmax,
              datetime_max = datetime_x + rollmax)]
    y[, `:=` (datetime_y_orig = datetime_y)] #this datetime col will be dropped

    #join x rows to y if y's datetime falls within the x range
    joined = y[x, on = .(datetime_y <= datetime_max,
                         datetime_y >= datetime_min)]
    joined = na.omit(joined, cols = 'datetime_y_orig') #drop rows without matches

    # joined = data.table(x=as.Date(c('2019-12-31', '2020-01-01', '2020-01-02',
    #                                 '2019-12-31', '2020-01-01', '2020-01-02',
    #                                 '2019-12-31', '2020-01-01', '2020-01-02')),
    #                     y=as.Date(c('2020-01-01', '2020-02-01', '2020-03-01',
    #                                 '2020-01-01', '2020-01-02', '2020-01-03',
    #                                 '2020-03-01', '2020-02-01', '2020-01-02')),
    #            a = c(NA, 1, 2, NA,NA,NA, NA,NA,3), b=c(1,1,1,2,2,2,3,3,3))
    # joined[, `:=` (datetime_match_diff = abs(x - y))]
    # joined[order(datetime_match_diff), lapply(.SD, function(xx) first(na.omit(xx))), by = x]

    #for any datetimes in x or y that were matched more than once, keep only
    #the nearest non-NA match by time
    joined[, `:=` (datetime_match_diff = abs(datetime_x - datetime_y_orig))]
    joined = joined[order(datetime_match_diff),
                    lapply(.SD, function(z) first(na.omit(z))),
                    by = datetime_x]

    joined[, c('datetime_y', 'datetime_y.1', 'datetime_y_orig', 'datetime_match_diff') := NULL]
    setnames(joined, 'datetime_x', 'datetime')

    joined = joined %>%
        as_tibble() %>%
        arrange(datetime) %>%
        rename(discharge = discharge_manual) %>%
        # mutate(season = factor(lubridate::quarter(datetime)),
        mutate(season = factor(!!seasons[month(joined$datetime)]),
               discharge_log = neglog(discharge)) %>%
        select(site_code, datetime, discharge, discharge_log, everything())

    return(joined)
}

assemble_q_df_daily <- function(site_code, ms_Q_data = NULL, scale_q_by_area = TRUE,
                                overwrite = FALSE){

    #ms_Q_data: data.frame with date column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    field_q = read_csv(glue('in/field_Q/{site_code}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE) %>%
        mutate(date = as.Date(datetime)) %>%
        select(-datetime)

    if(scale_q_by_area){
        wsa = filter(ms_areas, site_code == site_code) %>% pull(ws_area_ha)
        field_q$discharge_manual = field_q$discharge_manual / wsa * 1000
    }

    if(! file.exists(glue('in/usgs_Q/{site_code}.csv')) || overwrite){

        ms_sites = grep('date', colnames(ms_Q_data), invert = TRUE, value = TRUE)

        if(scale_q_by_area){

            for(g in ms_sites){

                wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
                ms_Q_data[[g]] = ms_Q_data[[g]] / wsa * 1000
            }
        }

        ms_Q_data = ms_Q_data %>%
            mutate(across(-date, neglog, .names = '{.col}_log'))

        write_csv(ms_Q_data, glue('in/usgs_Q/{site_code}.csv'))

    } else {
        ms_Q_data = read_csv(glue('in/usgs_Q/{site_code}.csv'))
    }

    joined = left_join(field_q, ms_Q_data, by = 'date') %>%
        arrange(date) %>%
        rename(discharge = discharge_manual) %>%
        # mutate(season = factor(lubridate::quarter(date)),
        mutate(season = factor(seasons[month(date)]),
               discharge_log = neglog(discharge)) %>%
        select(site_code, date, discharge, discharge_log, everything())

    return(joined)
}

cull_gauges <- function(d){

    d <- arrange(d, datetime)

    recording_intervals_m <- c()
    for(col in colnames(d)[-1]){

        dtcol <- d %>%
            select(datetime, !!col) %>%
            filter(! is.na(!!sym(col))) %>%
            pull(datetime)

        recording_intervals_m <- c(recording_intervals_m, Mode(diff(as.numeric(dtcol)) / 60))
    }

    weird_ints <- ! recording_intervals_m %in% c(5, 15, 30, 60)
    if(any(weird_ints)){
        warning('gauge recording intervals of 5, 15, 30, and 60 minutes are ',
                'expected, but interval(s) of ',
                paste(unique(recording_intervals_m[weird_ints]),
                      collapse = ', '),
                ' was/were encountered. Probably nbd.')
    }

    max_intvl <- max(recording_intervals_m)

    highres_col_ind <- which(recording_intervals_m < max_intvl) + 1
    highres_cols <- colnames(d)[highres_col_ind]

    cat(
        'Recording intervals (minutes) by gauge:\n\t',
        paste0(paste0(str_pad(paste0(colnames(d)[-1], ': '),
                              width = '12',
                              side = 'right')),
               paste0(recording_intervals_m, '\n\t'))
    )

    nrow0 <- nrow(d)
    d <- filter(d, ! if_all(-all_of(c('datetime', highres_cols)), is.na))

    if(any(recording_intervals_m != max_intvl)){
        message('Coercing all donor gauge series to the coarsest recording interval (',
                max_intvl, ' mins). Dropping ', nrow0 - nrow(d), ' rows.')
    }

    nas_by_gauge <- apply(d[, -1], 2, function(x) sum(is.na(x)))
    na_pcts_by_gauge <- apply(d[, -1], 2, function(x) round(sum(is.na(x)) / length(x), 2))
    na_pcts_by_gauge <- sub('^0$', '~0', na_pcts_by_gauge)

    cat(
        'Missing values by gauge (any timestep with a missing value cannot be predicted):\n\t',
        paste0(str_pad(paste0(names(nas_by_gauge), ': '),
                       width = '12',
                       side = 'right'),
               str_pad(paste0(unname(nas_by_gauge)),
                       width = '6',
                       side = 'right'),
               paste0(' (', unname(na_pcts_by_gauge), '%)\n\t'),
               collapse = ' ')
    )

    return(d)
}

# regression ####

glmnet_wrap <- function(data, full_spec, unscale_q_by_area = TRUE,
                        k_folds = 10, cv_ensemble = 10, plot_marg = TRUE,
                        alpha = 0, ...){

    #data: a data frame with a date column, and all columns referenced in model_list
    #full_spec: formula; the fully specified model to regularize via L2 (ridge regression)
    #cv_ensemble: the number of cross-validation procedures to carry out. optimal
    #   lambda varies according to how observations are randomly allocated to CV folds.
    #   so, select `cv_ensemble` optimal lambdas and average them.
    #plot_marg: logical. plot marginal relationships
    #alpha: 0 = ridge; 1 = lasso; see ?glmnet
    #...: additional arguments passed to glmnet AND cv.glmnet

    # prepare data
    x <- try(model.matrix(full_spec, data = data), silent = TRUE)
    if(inherits(x, 'try-error')){
        if(grepl('or more levels', attr(x, 'condition')$message)){
            message('Only one season represented in the field data. Removing seasonal term.')
            full_spec <- update(full_spec, . ~ . - season - .:season)
            x <- model.matrix(full_spec, data = data)
        } else stop(x)
    }
    mod_inds <- complete.cases(data)
    y <- data$discharge_log[mod_inds]
    data_filt <- data[mod_inds, ]

    ## get cross-validated predictions and metrics

    nobs <- length(y)
    if(nobs < k_folds){
        message('k_folds is ', k_folds, ', but there are only ', nobs,
                ' observations. Setting k_folds to ', nobs, ' (LOOCV).')
        k_folds <- min(nobs, k_folds)
    }

    is_grpd <- ifelse(nobs / k_folds < 3, FALSE, TRUE)

    best_mse <- Inf
    site_code <- data$site_code[1]
    wsa <- filter(sites, site_code == !!site_code) %>% pull(ws_area_ha)
    for(i in 0:1){

        n_samp <- length(y)
        # fold_ids <- sample(rep(1:10, length.out = n_samp), size = n_samp)

        sim_mat <- obs_mat <- matrix(NA_real_, nrow = n_samp, ncol = cv_ensemble)
        lambda_vec <- rep(NA_real_, cv_ensemble)
        for(j in 1:cv_ensemble){

            cv_model_ <- cv.glmnet(x, y, alpha = alpha, nfolds = k_folds, grouped = is_grpd,
                                   intercept = as.logical(i), type.measure = 'mse', ...)
            lambda_vec[j] <- cv_model_$lambda.min
            cvpred <- c(predict(cv_model_, s = 'lambda.min', newx = x))

            if(unscale_q_by_area){
                sim_mat[, j] <- inv_neglog(cvpred) * wsa / 1000 #un-scale by 1000. see assemble_q_df
                obs_mat[, j] <- data_filt$discharge * wsa / 1000
            } else {
                sim_mat[, j] <- inv_neglog(cvpred)
                obs_mat[, j] <- data_filt$discharge
            }
        }

        best_lambda_ <- mean(lambda_vec)
        sim_ <- apply(sim_mat, 1, mean, na.rm = TRUE)
        obs_ <- apply(obs_mat, 1, mean, na.rm = TRUE)

        # print(tibble(nse = hydroGOF::NSE(sim_, obs_),
        #              kge = hydroGOF::KGE(sim_, obs_)))

        new_mse <- hydroGOF::mse(sim_, obs_)
        if(new_mse < best_mse){
            best_mse <- new_mse
            # cv_model <- cv_model_
            sim <- sim_
            obs <- obs_
            best_lambda <- best_lambda_
            intcpt <- as.logical(i)
        }
    }

    metr_cv <- tibble(
        nse = hydroGOF::NSE(sim, obs),
        kge = hydroGOF::KGE(sim, obs),
        mse = hydroGOF::mse(sim, obs),
        mae = hydroGOF::mae(sim, obs),
        pbias = hydroGOF::pbias(sim, obs)
    )

    ## plot marginal associations

    trms <- rownames(attributes(terms(full_spec))$factors)
    dep <- gsub('`', '', trms[1])
    indeps <- gsub('`', '', trms[-1])
    site_indeps <- grep('season', indeps, invert = TRUE, value = TRUE)

    ggps <- list()
    for(i in seq_along(site_indeps)){
        d2 <- filter(data, ! is.na(discharge_log) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] <- ggplot(d2, aes(x = !!sym(site_indeps[i]), y = discharge_log)) +
            geom_point()
    }
    gd <- do.call("grid.arrange", c(ggps))
    if(plot_marg) plot(gd)

    ## handle novel seasonal factor values

    if('season' %in% indeps){
        modeled_seasons <- sort(as.character(unique(data_filt$season)))
        modseas_inds <- as.character(data_filt$season) %in% modeled_seasons
    } else {
        modeled_seasons <- as.character(1:4)
        modseas_inds <- rep(TRUE, nrow(data_filt))
    }

    ## run regression with optimal lambda and generate predictions + metrics

    best_model <- glmnet(x, y, alpha = alpha, lambda = best_lambda, intercept = intcpt, ...)

    data$lm <- NA_real_
    bestpred <- predict(best_model, s = best_lambda, newx = x[modseas_inds, ])

    if(unscale_q_by_area){
        data$lm[mod_inds][modseas_inds] <- inv_neglog(bestpred) * wsa / 1000
        data$discharge <- data$discharge * wsa / 1000
    } else {
        data$lm[mod_inds][modseas_inds] <- inv_neglog(bestpred)
    }

    data$lm[data$lm < 0] <- 0

    metr <- tibble(
        nse = hydroGOF::NSE(data$lm, data$discharge),
        kge = hydroGOF::KGE(data$lm, data$discharge),
        mse = hydroGOF::mse(data$lm, data$discharge),
        mae = hydroGOF::mae(data$lm, data$discharge),
        pbias = hydroGOF::pbias(data$lm, data$discharge)
    )

    ## trim leading/trailing NAs, clean up data for plotting
    first_non_na <- Position(function(x) ! is.na(x), data$lm)
    last_non_na <- nrow(data) - Position(function(x) ! is.na(x), rev(data$lm)) + 1
    plot_data <- data[first_non_na:last_non_na, ]

    site_indeps <- c(site_indeps, sub('_log', '', site_indeps))
    drop_cols <- grep('^[0-9]', colnames(plot_data), value = TRUE)
    drop_cols <- drop_cols[! drop_cols %in% site_indeps]
    plot_data <- select(plot_data, -any_of(drop_cols), -ends_with('_log')) %>%
        select(site_code, any_of(c('date', 'datetime')),
               Q_field = discharge, Q_predicted = lm,
               everything())

    ## unscale usgs Q

    if(unscale_q_by_area){

        nearby_usgs_gages <- grep('^[0-9]+$', colnames(plot_data), value = TRUE)
        for(g in nearby_usgs_gages){

            if(g == '06190540'){
                plot_data[[g]] <- plot_data[[g]] * 51089.73 #ha
                next
            }

            nwissite <- dataRetrieval::readNWISsite(g)

            if(is.null(nwissite)) stop('cannot find nwis watershed area')

            if(is.na(nwissite$contrib_drain_area_va)){
                wsa <- nwissite$drain_area_va * 258.999 #mi^2 -> ha
            } else {
                wsa <- nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
            }

            plot_data[[g]] <- plot_data[[g]] * wsa / 1000
        }

        nearby_ms_gages <- select(plot_data,
                                  -site_code,
                                  -any_of(c('date', 'datetime', 'season')),
                                  -starts_with('Q_'),
                                  -any_of(nearby_usgs_gages)) %>%
            colnames()

        for(g in nearby_ms_gages){
            wsa <- filter(sites, site_code == !!g) %>% pull(ws_area_ha)
            plot_data[[g]] <- plot_data[[g]] * wsa / 1000
        }
    }

    full_spec_ignoreint <- full_spec
    if(! intcpt){
        full_spec <- update(full_spec, ~ . - 1)
    }

    out <- list(best_model = full_spec,
                best_model_ignoreint = full_spec_ignoreint,
                best_model_object = best_model,
                alpha = alpha,
                prediction = unname(data$lm),
                lm_data = plot_data,
                metrics = metr,
                metrics_crossval = metr_cv,
                fits = gd)

    return(out)
}

bootstrap_ci_glmnet <- function(ncores, in_df, frm, best, has_intcpt, newx_){

    clst_type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
    clst <- makeCluster(spec = ncores, type = clst_type)
    registerDoParallel(clst)
    error_indicator <- tempfile()

    newx <- model.matrix(update(frm, NULL ~ .), newx_)

    on.exit(stopCluster(clst))

    bootstrap_samps <- parallel::parSapply(clst, 1:1000, function(x){

        if(file.exists(error_indicator)) return()

        set.seed(as.numeric(Sys.time()) + sample(1:1000, 1))

        if('season' %in% rownames(attributes(terms(frm))$factors)){
            resamp <- stratified_resample(in_df, 'season')
        } else {
            resamp <- in_df[sample(seq_len(nrow(in_df)), replace = TRUE), ]
        }

        x <- model.matrix(frm, data = resamp)
        mod_inds <- complete.cases(resamp)
        y <- resamp$discharge_log[mod_inds]

        s <- try(cv.glmnet(x, y, alpha = best$alpha, intercept = has_intcpt)$lambda.min,
                 silent = TRUE)

        if(inherits(s, 'try-error')){
            if(grepl('is constant', attr(s, 'condition')$message)){
                write('STOP', error_indicator)
                return()
            } else {
                stop(s)
            }
        }

        m <- glmnet(x, y, alpha = best$alpha, intercept = has_intcpt, lambda = s)

        predict(m, s = s, newx = newx)
    }, simplify = TRUE)

    if(sum(sapply(bootstrap_samps, is.null)) > 50){
        warning('Not enough y-variation to perform cross-validation')
        return(list(ci_lwr = rep(NA_real_, nrow(newx)),
                    ci_upr = rep(NA_real_, nrow(newx))))
    }

    ci <- parallel::parApply(clst, bootstrap_samps, 1, quantile, probs = c(0.025, 0.975))

    return(list(ci_lwr = ci[1, ], ci_upr = ci[2, ]))
}

generate_nested_formulae <- function(full_spec, d,
                                     interactions = TRUE, through_origin = TRUE,
                                     min_points_per_param = 15,
                                     max_interaction = Inf){

    #full_spec: a formula representing the full model specification (all possible predictors)
    #d: a data frame. must include dependent variable from full_spec as a column name
    #interactions: are interaction terms allowed? see max_interaction
    #through_origin: should 0-intercept models be generated?
    #min_points_per_param: at least this many data points must be present for each
    #   parameter included in resulting formulae.
    #max_interaction: the maximum number of terms that may interact, e.g. if 3,
    #   then formulae with x:y:z:w will be filtered from the results.

    #if one of the terms is "dummy" it will not be considered for interaction

    #returns a list of formulae

    trms = rownames(attributes(terms(full_spec))$factors)

    dummy_present <- 'dummy' %in% trms
    if(dummy_present){
        trms <- grep('dummy', trms, value = TRUE, invert = TRUE)
    }

    y = trms[1]
    y_nobacktick = gsub('`', '', y)
    trms = trms[! trms == y]
    ntrms = length(trms)

    indeps = list()
    for(i in seq_len(ntrms)){
        indeps = append(indeps, combn(trms, i, simplify = FALSE))
    }

    if(interactions){
        int_trms = Filter(function(x) length(x) > 1, indeps) %>%
            map(~paste(., collapse = ':')) %>%
            unlist()
    }

    max_params = sum(! is.na(d[[y_nobacktick]])) %/% min_points_per_param #https://statisticsbyjim.com/regression/overfitting-regression-models/
    indeps = Filter(function(x) length(x) <= max_params, indeps)

    if(! length(indeps)) stop('not enough data')

    mods = list()
    modind = 0
    for(i in seq_along(indeps)){

        modind = modind + 1
        mods[[modind]] = as.formula(glue('{y} ~ {paste(indeps[[i]], collapse = "+")}'))

        if(through_origin){
            modind = modind + 1
            mods[[modind]] = as.formula(glue('{y} ~ 0 + {paste(indeps[[i]], collapse = "+")}'))
        }

        if(interactions){

            valid_ints_for_indep = map(int_trms, ~str_split(., ':')[[1]]) %>%
                map(~all(. %in% indeps[[i]])) %>%
                unlist()
            ints_to_tack_on = int_trms[valid_ints_for_indep]
            ints_to_tack_on = Filter(function(x) length(strsplit(x, ':')[[1]]) < max_params, ints_to_tack_on)
            ints_to_tack_on = Filter(function(x) length(strsplit(x, ':')[[1]]) <= max_interaction, ints_to_tack_on)

            jj = min(length(ints_to_tack_on), (max_params - 2))
            for(j in seq_len(jj)){

                int_trm_grps = combn(ints_to_tack_on, j, simplify = FALSE)
                int_trm_grps = Filter(function(x) length(x) < (max_params - 2), int_trm_grps)

                for(grp in int_trm_grps){

                    modind = modind + 1
                    mods[[modind]] = as.formula(glue('{y} ~ {xs} + {xs_int}',
                                                     xs = paste(indeps[[i]], collapse = '+'),
                                                     xs_int = paste(grp, collapse = '+')))

                    if(through_origin){
                        modind = modind + 1
                        mods[[modind]] = as.formula(glue('{y} ~ 0 + {xs} + {xs_int}',
                                                         xs = paste(indeps[[i]], collapse = '+'),
                                                         xs_int = paste(grp, collapse = '+')))
                    }
                }
            }
        }
    }

    if(dummy_present){
        mods <- c(mods, lapply(mods, function(x) update(x, . ~ . + dummy)))
    }

    param_counts = sapply(mods, function(x) sum(grepl('season', attr(terms(x), 'term.labels')) * 2 + 1)) + 1
    is_through_origin = sapply(mods, function(x) grepl('^0 +', as.character(x)[3]))
    param_counts = param_counts - is_through_origin
    mods = mods[param_counts <= max_params]

    mods <- Filter(function(x){
        trms <- attr(terms(x), 'term.labels')
        ! all(trms %in% c('season', 'dummy'))
    }, mods)

    return(mods)
}

lm_wrap <- function(data, model_list,
                    unscale_q_by_area = TRUE, k_folds = 10){

    #data: a data frame with a date column, and all columns referenced in model_list
    #model_list: a list of lm formulae, as generated by generate_nested_formulae

    #lots of vestigial code in this function. sincere apologies.

    log = 'xy'
    best_mod = NA
    best_score = Inf
    for(i in seq_along(model_list)){

        try({

            trms = rownames(attributes(terms(model_list[[i]]))$factors)
            indeps = gsub('`', '', trms[-1])
            dep = gsub('`', '', trms[1])
            dep_transformed = sub('_log', '', dep)

            dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))

            tryCatch({
                dd = CVlm(data = dd, model_list[[i]], m = k_folds,
                          plotit = FALSE, printit = FALSE)
            }, warning = function(w) stop('rank-deficiency detected. moving on.'))

            if(log == 'xy'){

                site_code = data$site_code[1]
                if(unscale_q_by_area){
                    wsa = filter(ms_areas, site_code == site_code) %>% pull(ws_area_ha)
                    sim = inv_neglog(dd$cvpred) * wsa / 1000
                    # nsenum = sum((inv_neglog(dd$cvpred) * wsa / 1000 - inv_neglog(dd[[dep]]) * wsa / 1000)^2) #just checking
                } else {
                    sim = inv_neglog(dd$cvpred)
                    # nsenum = sum((inv_neglog(dd$cvpred) - inv_neglog(dd[[dep]]))^2)
                }
            }

            if(unscale_q_by_area){
                obs = pull(dd[dep_transformed]) * wsa / 1000
            } else {
                obs = pull(dd[dep_transformed])
            }

            # nsedenom = sum((obs - mean(obs))^2)
            # nse_ = 1 - (nsenum / nsedenom)

            # nse = hydroGOF::NSE(sim, obs)
            mserr = hydroGOF::mse(sim, obs)

            # if(! all.equal(nse, nse_)) stop('oi')

            if(mserr < best_score){
                best_score = mserr
                best_mod = i

                metr_cv = tibble(
                    nse = hydroGOF::NSE(sim, obs),
                    kge = hydroGOF::KGE(sim, obs),
                    mse = mserr,
                    mae = hydroGOF::mae(sim, obs),
                    pbias = hydroGOF::pbias(sim, obs)
                )
            }
        })
    }

    trms = rownames(attributes(terms(model_list[[best_mod]]))$factors)
    dep = gsub('`', '', trms[1])
    indeps = gsub('`', '', trms[-1])
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps = grep('season|dummy', indeps, invert = TRUE, value = TRUE)

    dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))
    m = lm(model_list[[best_mod]], data = dd)

    ggps = list()
    yvar = ifelse(log %in% c('y', 'xy'), 'discharge_log', 'discharge')
    for(i in seq_along(site_indeps)){
        d2 = filter(data, ! is.na(!!sym(yvar)) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] = ggplot(d2, aes(x = !!sym(site_indeps[i]), y = !!sym(yvar))) +
            geom_point()

        if(! 'dummy' %in% indeps){
            ggps[[i]] <- ggps[[i]] + stat_smooth(method = "lm", col = "red")
        }
    }
    gd = do.call("grid.arrange", c(ggps))
    print(gd)

    newdata = select(data, all_of(indeps))

    if('season' %in% indeps){
        modeled_seasons = m$xlevels$season
        modseas_inds = as.character(newdata$season) %in% modeled_seasons
    } else {
        modeled_seasons = as.character(1:4)
        modseas_inds = rep(TRUE, nrow(newdata))
    }

    data$lm = NA_real_

    if(log %in% c('y', 'xy')){

        if(unscale_q_by_area){
            data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ])) * wsa / 1000
            data$discharge = data$discharge * wsa / 1000
        } else {
            data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ]))
        }

    } else {
        stop('unused')
        data$lm[modseas_inds] = predict(m, newdata = newdata[modseas_inds, ])
    }

    data$lm[data$lm < 0] = 0

    nse_out = hydroGOF::NSE(data$lm, data$discharge)

    metr = tibble(
        nse = nse_out,
        kge = hydroGOF::KGE(data$lm, data$discharge),
        mse = hydroGOF::mse(data$lm, data$discharge),
        mae = hydroGOF::mae(data$lm, data$discharge),
        pbias = hydroGOF::pbias(data$lm, data$discharge)
    )

    first_non_na = Position(function(x) ! is.na(x), data$lm)
    last_non_na = nrow(data) - Position(function(x) ! is.na(x), rev(data$lm)) + 1
    plot_data = data[first_non_na:last_non_na, ]

    site_indeps = c(site_indeps, sub('_log', '', site_indeps))
    drop_cols = grep('^[0-9]', colnames(plot_data), value = TRUE)
    drop_cols = drop_cols[! drop_cols %in% site_indeps]
    plot_data = select(plot_data, -any_of(drop_cols), -ends_with('_log')) %>%
        select(site_code, any_of(c('date', 'datetime')),
               Q_field = discharge, Q_predicted = lm,
               everything())

    if(unscale_q_by_area){

        nearby_usgs_gages = grep('^[0-9]+$', colnames(plot_data), value = TRUE)
        for(g in nearby_usgs_gages){

            if(g == '06190540'){
                plot_data[[g]] = plot_data[[g]] * 51089.73 #ha
                next
            }

            nwissite = dataRetrieval::readNWISsite(g)

            if(is.null(nwissite)) stop('cannot find nwis watershed area')

            if(is.na(nwissite$contrib_drain_area_va)){
                wsa = nwissite$drain_area_va * 258.999 #mi^2 -> ha
            } else {
                wsa = nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
            }

            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }

        nearby_ms_gages = select(plot_data, -site_code, -any_of(c('date', 'datetime', 'dummy')),
                                 -starts_with('Q_'), -season, -any_of(nearby_usgs_gages)) %>%
            colnames()

        for(g in nearby_ms_gages){
            wsa = filter(sites, site_code == !!g) %>% pull(ws_area_ha)
            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }
    }

    out = list(best_model = model_list[[best_mod]],
               best_model_object = m,
               prediction = unname(data$lm),
               lm_data = plot_data,
               # score = nse_out,
               # score_crossval = best_score,
               metrics = metr,
               metrics_crossval = metr_cv,
               # plot = dg,
               fits = gd)

    return(out)
}

segmented_wrap <- function(data, model_list,
                           unscale_q_by_area = TRUE, k_folds = 10){

    #data: a data frame with a date column, and all columns referenced in model_list
    #model_list: a list of lm formulae, as generated by generate_nested_formulae

    best_mod <- NA
    best_score <- Inf
    for(i in seq_along(model_list)){

        try({

            trms <- rownames(attributes(terms(model_list[[i]]))$factors)
            indeps <- gsub('`', '', trms[-1])
            dep <- gsub('`', '', trms[1])
            dep_transformed <- sub('_log', '', dep)

            dd <- filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))

            tryCatch({
                dd <- CVlm(data = dd, model_list[[i]], m = k_folds,
                          plotit = FALSE, printit = FALSE)
            }, warning = function(w) stop('rank-deficiency detected. moving on.'))

            site_code <- data$site_code[1]
            if(unscale_q_by_area){
                wsa <- filter(sites, site_code == !!site_code) %>% pull(ws_area_ha)
                sim <- inv_neglog(dd$cvpred) * wsa / 1000
            } else {
                sim <- inv_neglog(dd$cvpred)
            }

            if(unscale_q_by_area){
                obs <- pull(dd[dep_transformed]) * wsa / 1000
            } else {
                obs <- pull(dd[dep_transformed])
            }

            mserr <- hydroGOF::mse(sim, obs)

            if(mserr < best_score){
                best_score <- mserr
                best_mod <- i
                metr_cv <- tibble(nse = NA_real_, kge = NA_real_, mse = mserr,
                                  mae = NA_real_, pbias = NA_real_)
            }
        })
    }

    trms <- rownames(attributes(terms(model_list[[best_mod]]))$factors)
    dep <- gsub('`', '', trms[1])
    indeps <- gsub('`', '', trms[-1])
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps <- grep('season', indeps, invert = TRUE, value = TRUE)

    #segmented package can't handle variable names with backticks
    model_list[[best_mod]] <- model_list[[best_mod]] %>%
        as.character() %>%
        str_replace_all('`?([0-9]{2,}(?:_log)?)`?', 'x\\1') %>%
        {paste(.[c(2, 1, 3)], collapse = ' ')} %>%
        as.formula()

    data_filt <- filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.))) %>%
        rename_with(~sub('`?([0-9]{2,}_log)`?', 'x\\1', .))

    seg_vars <- update(model_list[[best_mod]], NULL ~ . - season)
    if(! attr(terms(model_list[[best_mod]]), 'intercept')){
        seg_vars <- update(seg_vars, NULL ~ . + 1)
    }

    m <- lm(model_list[[best_mod]], data = data_filt) %>%
        segmented(seg.Z = seg_vars, npsi = 1)

    ggps = list()
    for(i in seq_along(site_indeps)){
        d2 = filter(data, ! is.na(discharge_log) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] = ggplot(d2, aes(x = !!sym(site_indeps[i]), y = discharge_log)) +
            geom_point()
    }
    gd = do.call("grid.arrange", c(ggps))
    print(gd)

    newdata = select(data, all_of(indeps)) %>%
        rename_with(~sub('`?([0-9]{2,}_log)`?', 'x\\1', .))

    if('season' %in% indeps){
        modeled_seasons = m$xlevels$season
        modseas_inds = as.character(newdata$season) %in% modeled_seasons
    } else {
        modeled_seasons = as.character(1:4)
        modseas_inds = rep(TRUE, nrow(newdata))
    }

    data$lm = NA_real_

    if(unscale_q_by_area){
        data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ])) * wsa / 1000
        data$discharge = data$discharge * wsa / 1000
    } else {
        data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ]))
    }

    data$lm[data$lm < 0] = 0

    metr = tibble(
        nse = hydroGOF::NSE(data$lm, data$discharge),
        kge = hydroGOF::KGE(data$lm, data$discharge),
        mse = hydroGOF::mse(data$lm, data$discharge),
        mae = hydroGOF::mae(data$lm, data$discharge),
        pbias = hydroGOF::pbias(data$lm, data$discharge)
    )

    first_non_na = Position(function(x) ! is.na(x), data$lm)
    last_non_na = nrow(data) - Position(function(x) ! is.na(x), rev(data$lm)) + 1
    plot_data = data[first_non_na:last_non_na, ]

    site_indeps = c(site_indeps, sub('_log', '', site_indeps))
    drop_cols = grep('^[0-9]', colnames(plot_data), value = TRUE)
    drop_cols = drop_cols[! drop_cols %in% site_indeps]
    plot_data = select(plot_data, -any_of(drop_cols), -ends_with('_log')) %>%
        select(site_code, any_of(c('date', 'datetime')),
               Q_field = discharge, Q_predicted = lm,
               everything())

    if(unscale_q_by_area){

        nearby_usgs_gages = grep('^[0-9]+$', colnames(plot_data), value = TRUE)
        for(g in nearby_usgs_gages){

            if(g == '06190540'){
                plot_data[[g]] = plot_data[[g]] * 51089.73 #ha
                next
            }

            nwissite = dataRetrieval::readNWISsite(g)

            if(is.null(nwissite)) stop('cannot find nwis watershed area')

            if(is.na(nwissite$contrib_drain_area_va)){
                wsa = nwissite$drain_area_va * 258.999 #mi^2 -> ha
            } else {
                wsa = nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
            }

            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }

        nearby_ms_gages = select(plot_data, -site_code, -any_of(c('date', 'datetime')),
                                 -starts_with('Q_'), -season, -any_of(nearby_usgs_gages)) %>%
            colnames()

        for(g in nearby_ms_gages){
            wsa = filter(sites, site_code == !!g) %>% pull(ws_area_ha)
            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }
    }

    out = list(best_model = model_list[[best_mod]],
               best_model_object = m,
               prediction = unname(data$lm),
               lm_data = plot_data,
               metrics = metr,
               metrics_crossval = metr_cv,
               fits = gd)

    return(out)
}

regress <- function(site_code, framework, ..., scale_q_by_area = TRUE,
                    precomputed_df = NULL, custom_formula = NULL,
                    dummy_break = NULL, custom_gauges = NULL,
                    seasons = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4),
                    no_write = FALSE, bootstrap_ci = TRUE,
                    target_daterange = NULL,
                    ncores = max(parallel::detectCores() %/% 4, 2)){

    #site_code: a site code included in cfg/sites.csv
    #framework: string. 'lm', 'glmnet', or 'segmented'
    #...: other arguments to lm_wrap, glmnet_wrap, or segmented_wrap
    #     depending on the value of `framework`
    #scale_q_by_area: logical. should discharge be scaled by watershed area?
    #precomputed_df: data.frame as produced by assemble_q_df or
    #   assemble_q_df_daily. only supply this if you need to make modifications.
    #   otherwise it will be generated automatically by this function. may
    #   include a factor column called "dummy" that contains a dummy variable
    #custom_formula: a `glue` specification . if not provided, the fully specified model
    #   is assumed. See function definition for examples. only set up to work
    #   with special cases of framework == 'lm'.
    #dummy_break: if a dummy variable is used, this specifies the date on which to
    #   shift from dummy = 0 to dummy = 1. ignored otherwise
    #custom_gauges: character vector. use if you need to override the defaults
    #   provided by cft/donor_gauges.yml, e.g. if for building composite models.
    #seasons: a numeric vector of length 12, indicating how months should be
    #   grouped into seasons. The first element corresponds to January, the last
    #   to December.
    #no_write: logical. if TRUE, do not write to the `results` data.frame in the
    #   global environment, and do not write data or plots via plots_and_results.
    #bootstrap_ci: logical. should 95% confidence intervals for glmnet models
    #   be computed? Ignored if framework != 'glmnet'. see ncores.
    #target_daterange: vector of two dates for the beginning and end of the
    #   period to be estimated. Donor gauges with insufficient data may be dropped before
    #   modeling, but those gauges might actually have data for your daterange
    #   of interest, so use this parameter to get the most out of the available donor gauges.
    #ncores: integer. number of cores to use for parallel computation of 95%
    #   confidence intervals

    #updates `results` in memory and out/lm_out/results_specificq.csv (unless no_write is TRUE)
    #   if scale_q_by_area is FALSE, writes instead to out/lm_out/results.csv
    #returns list including:
    #   best_model: chosen model, as a formula
    #   best_model_ignoreint: glmnet only. same thing as best_model, but if
    #       regression through the origin was successful, this won't indicate it.
    #       used internally.
    #   best_model_object: the glmnet, lm, or segmented fit object
    #   prediction: vector of best-model predictions for each observation
    #   lm_data: a tibble of inputs and outputs. misnomer in case of glmnet or segmented
    #   metrics: nse, kge, rmse, mae, pbias for best model
    #   metrics_crossval: the same, but computed on out-of-sample folds during
    #       cross-validation
    #   fits: misnomer. just a plot of the marginal associations between each
    #       predictor gauge and the response. regression line included if model
    #       framework == 'lm'

    if(! framework %in% c('lm', 'glmnet', 'segmented')){
        stop('framework must be "lm", "glmnet", or "segmented"')
    }

    if(is.null(custom_gauges)){
        if(! site_code %in% names(donor_gauges)){
            stop(site_code, ' entry missing from cfg/donor_gauges.yml')
        }
        gage_ids <- donor_gauges[[site_code]]
    } else {
        gage_ids <- custom_gauges
    }

    if(is.null(precomputed_df)){
        in_df <- assemble_q_df(site_code = site_code,
                               nearby_usgs_gages = gage_ids,
                               scale_q_by_area = scale_q_by_area,
                               target_daterange = target_daterange,
                               seasons = seasons)
    } else {
        in_df <- precomputed_df
    }

    formulaA <- 'discharge_log ~ `{paste(paste0(gage_ids, "_log"), collapse = "`+`")}` + season'
    formulaB <- 'discharge_log ~ `{paste(paste0(gage_ids, "_log"), collapse = "`*`")}` * season'
    if(! is.null(custom_formula)){
        formulaA <- formulaB <- custom_formula
    }

    if(framework == 'lm'){

        mods <- generate_nested_formulae(
            full_spec = as.formula(glue(formulaA)),
            d = in_df,
            max_interaction = 3
        )

        best <- lm_wrap(
            data = in_df,
            model_list = mods,
            unscale_q_by_area = scale_q_by_area,
            ...
        )

    } else if(framework == 'segmented'){

        mods <- generate_nested_formulae(
            full_spec = as.formula(glue(formulaA)),
            d = in_df,
            max_interaction = 3
        )

        best <- segmented_wrap(
            data = in_df,
            model_list = mods,
            unscale_q_by_area = scale_q_by_area,
            ...
        )

    } else { #glmnet

        best <- glmnet_wrap(
            data = in_df,
            full_spec = as.formula(glue(formulaB)),
            unscale_q_by_area = scale_q_by_area,
            ...
        )
    }

    if(! no_write){

        results <- plots_and_results(
            site_code, best, results, in_df,
            unscale_q_by_area = scale_q_by_area,
            dummy_break = dummy_break,
            seasons = seasons,
            bootstrap_ci = bootstrap_ci,
            ncores = ncores
        )
        results[results$site_code == site_code, 'method'] <- framework
        results <<- results

        if(scale_q_by_area){
            write_csv(results, 'out/lm_out/results_specificq.csv')
        } else {
            write_csv(results, 'out/lm_out/results.csv')
        }
    }

    return(best)
}

plots_and_results <- function(site_code, best, results, in_df,
                              return_plot = FALSE, unscale_q_by_area = TRUE,
                              seasons,
                              dummy_break = NULL, bootstrap_ci, ncores){

    #load corroborating usgs/ms site data
    sites_nearby = read_csv(glue('in/usgs_Q/{site_code}.csv'))

    if(inherits(best$best_model_object, 'segmented')){
        sites_nearby <- rename_with(sites_nearby, ~sub('`?([0-9]{2,}(?:_log)?)`?', 'x\\1', .))
    }

    trms = rownames(attributes(terms(best$best_model))$factors)

    dummy_present <- 'dummy' %in% trms
    if(dummy_present){
        trms <- grep('dummy', trms, value = TRUE, invert = TRUE)
        if(is.null(dummy_break)) stop('dummy_break must be supplied if dummy variable present')
    }
    dep = gsub('`', '', trms[1])
    indeps = gsub('`', '', trms[-1])
    indeps = gsub('_log', '', indeps)
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps = grep('season', indeps, invert = TRUE, value = TRUE)
    site_indeps_log = paste0(site_indeps, '_log')

    sites_nearby = sites_nearby %>%
        select(datetime, all_of(site_indeps), all_of(site_indeps_log))

    if('season' %in% indeps){
        # sites_nearby$season = factor(lubridate::quarter(sites_nearby$datetime))
        sites_nearby$season = factor(seasons[month(sites_nearby$datetime)])
    }

    # #assemble neon sensor data, filtered neon sensor data, neon field data
    # #into one frame and plot it
    #
    # neon_q_auto = read_csv(glue('in/neon_continuous_Q/{site_code}.csv')) %>%
    #     filter(! is.na(discharge)) %>%
    #     rename(discharge_auto = discharge)
    #
    # #predict Q for all datetimes with predictor data
    #
    # qall = left_join(sites_nearby, neon_q_auto, by = 'datetime')
    qall = sites_nearby

    if(dummy_present){
        qall$dummy <- 0
        qall$dummy[qall$datetime > dummy_break] <- 1
        qall$dummy <- as.factor(qall$dummy)
    }

    if(inherits(best$best_model_object, 'segmented')){

        qall <- predict(best$best_model_object,
                        newdata = select(qall, all_of(site_indeps_log), any_of('season')),
                        interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog)) %>%
            bind_cols(qall)

    } else if(inherits(best$best_model_object, 'glmnet')){

        betas <- best$best_model_object$beta@Dimnames[[1]]

        newx_ <- select(qall, all_of(site_indeps_log), any_of('season'))

        if('season' %in% colnames(qall)){

            seasons_present <- str_extract(betas, 'season[1-4]') %>%
                unique() %>%
                na.omit() %>%
                str_split_i('season', 2)

            seasons_present <- c(seasons_present,
                                 min(as.numeric(seasons_present)) - 1)

            newx_ <- newx_ %>%
                mutate(season = ifelse(! season %in% seasons_present, NA, season),
                       season = as.factor(season))
        }

        valid_obs <- complete.cases(newx_)
        frm <- best$best_model_ignoreint
        has_intcpt <- '(Intercept)' %in% betas

        pred <- predict(
            best$best_model_object,
            s = best$best_model_object$lambda,
            newx = model.matrix(update(frm, NULL ~ .), newx_),
            type = 'response'
        )

        if(bootstrap_ci){
            cat('Bootstrapping 95% confidence intervals on', ncores, 'cores. Watch your RAM!\n')
            ci <- bootstrap_ci_glmnet(ncores, in_df, frm, best, has_intcpt, newx_)
        } else {
            ci <- list(ci_lwr = rep(NA_real_, length(pred)),
                       ci_upr = rep(NA_real_, length(pred)))
        }

        pred[pred < 0] <- 0
        ci$ci_lwr[ci$ci_lwr < 0] <- 0
        ci$ci_upr[ci$ci_upr < 0] <- 0

        qall$upr <- qall$lwr <- qall$fit <- NA_real_
        qall[valid_obs, 'fit'] <- inv_neglog(pred)
        qall[valid_obs, 'lwr'] <- inv_neglog(ci$ci_lwr)
        qall[valid_obs, 'upr'] <- inv_neglog(ci$ci_upr)

    } else { #lm

        qall <- predict(best$best_model_object,
                        newdata = qall,
                        interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog)) %>%
            bind_cols(qall)
    }

    if(unscale_q_by_area){
        wsa = filter(sites, site_code == !!site_code) %>% pull(ws_area_ha)
        qall = mutate(qall, fit = fit * wsa / 1000, lwr = lwr * wsa / 1000, upr = upr * wsa / 1000)
    }

    out_data = qall %>%
        select(-ends_with('_log'), -any_of('season')) %>%
        full_join(select(best$lm_data, datetime, Q_field), by = 'datetime') %>%
        select(datetime, Q_predicted = fit, #Q_used_in_regression = discharge,
               Q_pred_int_2.5 = lwr, Q_pred_int_97.5 = upr, Q_field) %>%
               # Q_neon_continuous_filtered = discharge_auto_qc,
               #Q_neon_continuous = discharge_auto) %>%
        filter(if_any(-datetime, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(datetime)

    dg = dygraphs::dygraph(xts(x = select(out_data, -datetime, -Q_pred_int_2.5, -Q_pred_int_97.5) %>% tail(5e5),
                               order.by = tail(out_data$datetime, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('figs/lm_plots/pred/{site_code}_log.html'))

    if(inherits(best$fits, 'grob')){

        #make diagnostic and "fit" plots

        png(glue('figs/lm_plots/diag/{site_code}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        defpar = par(mfrow=c(2, 2))
        if(inherits(best$best_model_object, 'segmented')){
            m <- best$best_model_object
            plot(m, conf.level=.95, is=TRUE, isV=FALSE, col=2, shade = TRUE, res=TRUE, res.col=4, pch=3)
            plot(fitted(m), m$residuals, xlab = 'Fitted values', ylab = 'Residuals')
            qqnorm(m$residuals, ylab = 'non-standardized residuals')
            qqline(m$residuals, lty = 3, col = 'gray40', lwd = 2)
        } else if(inherits(best$best_model_object, 'glmnet')){
            plot(best$best_model_object, xvar = 'norm', label = TRUE)
            plot(best$best_model_object, xvar = 'lambda', label = TRUE)
            plot(best$best_model_object, xvar = 'dev', label = TRUE)
        } else {
            plot(best$best_model_object)
        }
        par(defpar)
        dev.off()

        png(glue('figs/lm_plots/fit/{site_code}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        plot(best$fits)
        dev.off()
    }

    #plot predictions versus field measurements. need to round field meas to Q interval
    time_diffs = diff(sites_nearby$datetime)
    units(time_diffs) = 'mins'
    time_interval = Mode(time_diffs)

    zz = out_data %>%
        mutate(datetime = round_date(datetime, paste(time_interval, 'min')))

    field_dts = filter(zz, ! is.na(Q_field)) %>% pull(datetime)

    zz = zz %>%
        filter(datetime %in% field_dts) %>%
        group_by(datetime) %>%
        summarize(across(everything(), ~mean(., na.rm = TRUE))) %>%
        ungroup()

    axlim = c(0, max(c(zz$Q_predicted, zz$Q_field), na.rm = TRUE))

    asterisk <- if(inherits(best$best_model_object, 'segmented')) '*' else ''
    png(glue('figs/lm_plots/val/{site_code}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(zz$Q_field, zz$Q_predicted, xlab = 'Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)',
         main = glue('Site: {site_code}; KGE: {kge1}; KGE crossval{a}: {kge2}',
                     kge1 = round(best$metrics$kge, 2),
                     kge2 = round(best$metrics_crossval$kge, 2),
                     a = asterisk),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()

    out_data = filter(out_data, ! is.na(Q_predicted)) %>% select(-Q_field)

    #save predictions as CSV
    out_data = left_join(out_data,
                         select(sites_nearby, datetime, all_of(site_indeps), any_of('season')),
                         by = 'datetime')
    if(dummy_present){
        out_data$dummy <- 0
        out_data$dummy[out_data$datetime > dummy_break] <- 1
        out_data$dummy <- as.factor(out_data$dummy)
    }

    write_csv(out_data, glue('out/lm_out/predictions/{site_code}.csv'))

    #save fit data as CSV
    write_csv(best$lm_data, glue('out/lm_out/fit/{site_code}.csv'))

    #save model summary
    sink(glue('out/lm_out/summary/{site_code}.txt'))
    if(inherits(best$best_model_object, 'glmnet')){
        print(best$best_model_object)
        cat('\n---\n')
        print(summary(best$best_model_object))
        cat('\n---\n')
        print(best$best_model_object$beta)
    } else {
        print(summary(best$best_model_object))
    }
    sink()

    #return results frame, updated with this site's results
    results$bestmod[results$site_code == site_code] = as.character(best$best_model)[3]
    results$nse[results$site_code == site_code] = best$metrics$nse
    results$nse_cv[results$site_code == site_code] = best$metrics_crossval$nse
    results$kge[results$site_code == site_code] = best$metrics$kge
    results$kge_cv[results$site_code == site_code] = best$metrics_crossval$kge
    results$pbias[results$site_code == site_code] = best$metrics$pbias
    results$pbias_cv[results$site_code == site_code] = best$metrics_crossval$pbias
    if(! inherits(best$best_model_object, 'glmnet')){
        results$adj_r_squared[results$site_code == site_code] = summary(best$best_model_object)$adj.r.squared
    }

    if(unscale_q_by_area){
        write_csv(results, 'out/lm_out/results_specificq.csv')
    } else {
        write_csv(results, 'out/lm_out/results.csv')
    }

    if(return_plot){
        return(list(plot = dg, results = results))
    } else {
        return(results)
    }
}

plots_and_results_daily_composite <- function(site_code, best1, best2, results,
                                              unscale_q_by_area = TRUE, in_df1,
                                              in_df2, bootstrap_ci, ncores){

    #best2 and lm_df2 should represent the model with more terms included

    if(! setequal(class(best1), class(best2))) stop('best1 and best2 must have the same class(es)')

    #load corroborating usgs/ms site data
    sites_nearby = read_csv(glue('in/usgs_Q/{site_code}.csv'))

    #filter usgs/ms site data that didn't end up in the model
    trms1 = rownames(attributes(terms(best1$best_model))$factors)
    dep1 = gsub('`', '', trms1[1])
    indeps1 = gsub('`', '', trms1[-1])
    indeps1 = gsub('_log', '', indeps1)
    if(length(indeps1) == 1 && indeps1 == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps1 = grep('season', indeps1, invert = TRUE, value = TRUE)
    site_indeps_log1 = paste0(site_indeps1, '_log')

    trms2 = rownames(attributes(terms(best2$best_model))$factors)
    dep2 = gsub('`', '', trms2[1])
    indeps2 = gsub('`', '', trms2[-1])
    indeps2 = gsub('_log', '', indeps2)
    if(length(indeps2) == 1 && indeps2 == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps2 = grep('season', indeps2, invert = TRUE, value = TRUE)
    site_indeps_log2 = paste0(site_indeps2, '_log')

    if(length(site_indeps2) <= length(site_indeps1)) stop('model 2 should have more sites than 1')

    sites_nearby = sites_nearby %>%
        select(date, all_of(site_indeps2), all_of(site_indeps_log2))

    if('season' %in% indeps1 | 'season' %in% indeps2){
        # sites_nearby$season = factor(lubridate::quarter(sites_nearby$date))
        sites_nearby$season = factor(seasons[month(sites_nearby$datetime)])
    }

    #assemble neon sensor data, filtered neon sensor data, neon field data
    #into one frame and plot it

    neon_q_daily = read_csv(glue('in/neon_continuous_Q/{site_code}.csv')) %>%
        group_by(date = date(datetime)) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_daily = discharge)

    #predict Q for all datetimes with predictor data

    qall = left_join(sites_nearby, neon_q_daily, by = 'date')
        # left_join(neon_q_daily_qc, by = 'date')

    if(inherits(best1$best_model_object, 'lm')){

        pred1 = predict(best1$best_model_object,
                        newdata = select(qall, all_of(site_indeps_log1), any_of('season')),
                        interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog))

        pred = predict(best2$best_model_object,
                       newdata = select(qall, all_of(site_indeps_log2), any_of('season')),
                       interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog))

        stop('nothing below this point has been tested to work with lm since glmnet took over at COMO')

    } else if(inherits(best1$best_model_object, 'glmnet')){

        qall1 <- qall2 <- qall

        #mod1
        betas <- best1$best_model_object$beta@Dimnames[[1]]
        seasons_present <- str_extract(betas, 'season[1-4]') %>%
            unique() %>%
            na.omit() %>%
            str_split_i('season', 2)

        seasons_present <- c(seasons_present,
                             min(as.numeric(seasons_present)) - 1)

        newx_ <- select(qall, all_of(site_indeps_log1), any_of('season')) %>%
            mutate(season = ifelse(! season %in% seasons_present, NA, season),
                   season = as.factor(season))
        valid_obs <- complete.cases(newx_)
        frm <- best1$best_model_ignoreint
        has_intcpt <- '(Intercept)' %in% betas

        pred1 <- predict(
            best1$best_model_object,
            s = best1$best_model_object$lambda,
            newx = model.matrix(update(frm, NULL ~ .), newx_),
            type = 'response'
        )

        if(bootstrap_ci){
            cat('Bootstrapping 95% confidence intervals on', ncores, 'cores. Watch your RAM!')
            ci <- bootstrap_ci_glmnet(ncores, in_df1, frm, best1, has_intcpt, newx_)
        } else {
            ci <- list(ci_lwr = rep(NA_real_, length(pred1)),
                       ci_upr = rep(NA_real_, length(pred1)))
        }

        pred1[pred1 < 0] <- 0
        ci$ci_lwr[ci$ci_lwr < 0] <- 0
        ci$ci_upr[ci$ci_upr < 0] <- 0

        qall1$upr <- qall1$lwr <- qall1$fit <- NA_real_
        qall1[valid_obs, 'fit'] <- inv_neglog(pred1)
        qall1[valid_obs, 'lwr'] <- inv_neglog(ci$ci_lwr)
        qall1[valid_obs, 'upr'] <- inv_neglog(ci$ci_upr)

        #mod2
        betas <- best2$best_model_object$beta@Dimnames[[1]]
        seasons_present <- str_extract(betas, 'season[1-4]') %>%
            unique() %>%
            na.omit() %>%
            str_split_i('season', 2)

        seasons_present <- c(seasons_present,
                             min(as.numeric(seasons_present)) - 1)

        newx_ <- select(qall, all_of(site_indeps_log2), any_of('season')) %>%
            mutate(season = ifelse(! season %in% seasons_present, NA, season),
                   season = as.factor(season))
        valid_obs <- complete.cases(newx_)
        frm <- best2$best_model_ignoreint
        has_intcpt <- '(Intercept)' %in% betas

        pred2 <- predict(
            best2$best_model_object,
            s = best2$best_model_object$lambda,
            newx = model.matrix(update(frm, NULL ~ .), newx_),
            type = 'response'
        )

        if(bootstrap_ci){
            cat('Bootstrapping 95% confidence intervals on', ncores, 'cores. Watch your RAM!')
            ci <- bootstrap_ci_glmnet(ncores, in_df2, frm, best2, has_intcpt, newx_)
        } else {
            ci <- list(ci_lwr = rep(NA_real_, length(pred2)),
                       ci_upr = rep(NA_real_, length(pred2)))
        }

        pred2[pred2 < 0] <- 0
        ci$ci_lwr[ci$ci_lwr < 0] <- 0
        ci$ci_upr[ci$ci_upr < 0] <- 0

        qall2$upr <- qall2$lwr <- qall2$fit <- NA_real_
        qall2[valid_obs, 'fit'] <- inv_neglog(pred2)
        qall2[valid_obs, 'lwr'] <- inv_neglog(ci$ci_lwr)
        qall2[valid_obs, 'upr'] <- inv_neglog(ci$ci_upr)

    } else stop('this function only specified for lm and glmnet models')

    #qall[complete.cases(newx_), c('fit', 'lwr', 'upr')]
    #qall$upr <- qall$lwr <- qall$fit <- NA_real_

    pred1 <- select(qall1, fit, lwr, upr)
    pred <- select(qall2, fit, lwr, upr)

    missing_preds = is.na(pred$fit)
    pred[missing_preds, ] = pred1[missing_preds, ]
    qall = bind_cols(qall, pred)

    if(unscale_q_by_area){
        wsa = filter(sites, site_code == !!site_code) %>% pull(ws_area_ha)
        qall = mutate(qall, fit = fit * wsa / 1000, lwr = lwr * wsa / 1000, upr = upr * wsa / 1000)
    }

    out_data = qall %>%
        select(-ends_with('_log'), -any_of('season')) %>%
        full_join(select(best2$lm_data, date, Q_field), by = 'date') %>%
        select(date, Q_predicted = fit,
               Q_pred_int_2.5 = lwr, Q_pred_int_97.5 = upr, Q_field,
               # Q_neon_continuous_filtered = discharge_daily_qc,
               Q_neon_continuous = discharge_daily) %>%
        filter(if_any(-date, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(date)

    dg = dygraphs::dygraph(xts(x = select(out_data, -date, Q_pred_int_2.5, -Q_pred_int_97.5) %>% tail(5e5),
                               order.by = tail(out_data$date, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('figs/lm_plots/pred/{site_code}_log.html'))

    png(glue('figs/lm_plots/diag/{site_code}_diag_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    par(defpar)
    dev.off()

    #make diagnostic and fit plots
    png(glue('figs/lm_plots/diag/{site_code}_diag_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    if(inherits(best1$best_model_object, 'glmnet')){
        plot(best1$best_model_object, xvar = 'norm', label = TRUE)
        plot(best1$best_model_object, xvar = 'lambda', label = TRUE)
        plot(best1$best_model_object, xvar = 'dev', label = TRUE)
    } else {
        plot(best1$best_model_object)
    }
    par(defpar)
    dev.off()

    png(glue('figs/lm_plots/diag/{site_code}_diag_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    if(inherits(best2$best_model_object, 'glmnet')){
        plot(best2$best_model_object, xvar = 'norm', label = TRUE)
        plot(best2$best_model_object, xvar = 'lambda', label = TRUE)
        plot(best2$best_model_object, xvar = 'dev', label = TRUE)
    } else {
        plot(best2$best_model_object)
    }
    par(defpar)
    dev.off()

    png(glue('figs/lm_plots/fit/{site_code}_fit_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(best1$fits)
    dev.off()

    png(glue('figs/lm_plots/fit/{site_code}_fit_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(best2$fits)
    dev.off()

    #plot predictions versus field measurements. need to round field meas to Q interval
    nse = hydroGOF::NSE(out_data$Q_predicted, out_data$Q_field)
    kge = hydroGOF::KGE(out_data$Q_predicted, out_data$Q_field)
    pbias = hydroGOF::pbias(out_data$Q_predicted, out_data$Q_field)
    axlim = c(0, max(c(out_data$Q_predicted, out_data$Q_field), na.rm = TRUE))

    png(glue('figs/lm_plots/val/{site_code}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(out_data$Q_field, out_data$Q_predicted, xlab = 'Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)',
         main = glue('Site: {site_code}; KGE: {kge1}; KGE crossval*: {kge2}',
                     kge1 = round(kge, 2),
                     kge2 = NA),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()

    #save predictions as CSV
    out_data = filter(out_data, is.na(Q_field)) %>% select(-Q_field)
    out_data = left_join(out_data,
                         select(sites_nearby, date, all_of(site_indeps2), any_of('season')),
                         by = 'date')
    write_csv(out_data, glue('out/lm_out/predictions/{site_code}.csv'))

    #save fit data as CSV
    write_csv(best2$lm_data, glue('out/lm_out/fit/{site_code}.csv'))

    #save model summary
    sink(glue('out/lm_out/summary/{site_code}.txt'))
    if(inherits(best1$best_model_object, 'glmnet')){
        print(best1$best_model_object)
        cat('\n---\n')
        print(summary(best1$best_model_object))
        cat('\n---\n')
        print(best1$best_model_object$beta)

        cat('\n\n--- MODEL B (2) ---\n\n')

        print(best2$best_model_object)
        cat('\n---\n')
        print(summary(best2$best_model_object))
        cat('\n---\n')
        print(best2$best_model_object$beta)
    } else {
        print(summary(best1$best_model_object))
        cat('\n---\n')
        print(summary(best2$best_model_object))
    }
    sink()

    #return results frame, updated with this site's results
    results$bestmod[results$site_code == site_code] = paste(as.character(best1$best_model)[3],
                                                            as.character(best1$best_model)[3],
                                                            sep = '  &  ')
    results$nse[results$site_code == site_code] = nse
    results$nse_cv[results$site_code == site_code] = mean(best1$metrics_crossval$nse, best2$metrics_crossval$nse)
    results$adj_r_squared[results$site_code == site_code] = NA
    # results$adj_r_squared[results$site_code == site_code] = mean(summary(best1$best_model_object)$adj.r.squared,
    #                                                              summary(best2$best_model_object)$adj.r.squared)

    results$kge[results$site_code == site_code] = kge
    results$kge_cv[results$site_code == site_code] = mean(best1$metrics$kge, best2$metrics$kge)
    results$pbias[results$site_code == site_code] = pbias
    results$pbias_cv[results$site_code == site_code] = mean(best1$metrics$pbias, best2$metrics$pbias)

    if(unscale_q_by_area){
        write_csv(results, 'out/lm_out/results_specificq.csv')
    } else {
        write_csv(results, 'out/lm_out/results.csv')
    }

    return(results)
}

stratified_resample <- function(df, factor){

    indices <- split(seq_len(nrow(df)), df[[factor]])
    resampled_indices <- lapply(indices, sample, replace = TRUE)
    resampled_indices <- unlist(resampled_indices)

    return(df[resampled_indices,])
}

# misc helpers required to generate figures/datasets ####

dms_to_decdeg <- function(x){

    #x: an integer or character vector of latlongs in dms format, e.g. "123456" or 1234567
    #                                                                   DDMMSS     DDDMMSS

    decdeg = c()
    for(i in seq_along(x)){

        xx = x[i]

        if(is.na(xx)){
            decdeg[i] = NA
            next
        }

        if(! nchar(xx) %in% 6:7){
            warning(paste(nchar(xx), 'characters in x. need 6 (or 7 for some longitudes)'))
            decdeg[i] = NA
        }

        deginc = if(nchar(xx) == 7) 1 else 0

        degs = as.numeric(substr(xx, 1, 2 + deginc))
        mins = as.numeric(substr(xx, 3 + deginc, 4 + deginc))
        secs = as.numeric(substr(xx, 5 + deginc, 6 + deginc))

        decdeg[i] = degs + mins / 60 + secs / 3600
    }

    return(decdeg)
}

gaugeid_to_location <- function(id){

    out <- dataRetrieval::readNWISsite(id)
    lat <- dms_to_decdeg(as.integer(out$lat_va))
    long <- -dms_to_decdeg(as.integer(out$long_va))

    sfobj <- tibble(lat = lat, long = long) %>%
        st_as_sf(coords = c('long', 'lat'), crs = 4269) %>%
        st_transform(crs = 4326)

    return(sfobj)
}

get_osm_streams <- function(extent_raster, outfile = NULL){

    #extent_raster: either a terra spatRaster or a rasterLayer. The output
    #   streams will have the same crs, and roughly the same extent, as this raster.
    #outfile: string. If supplied, output shapefile will be written to this
    #   location. If not supplied, the output will be returned.

    message('Downloading streams layer from OpenStreetMap')

    extent_raster <- terra::rast(extent_raster)
    # rast_crs <- as.character(extent_raster@crs)
    rast_crs <- terra::crs(extent_raster,
                           proj = TRUE)

    extent_raster_wgs84 <- terra::project(extent_raster,
                                          y = 'epsg:4326')

    dem_bounds <- terra::ext(extent_raster_wgs84)[c(1, 3, 2, 4)]

    streams_query <- osmdata::opq(dem_bounds) %>%
        osmdata::add_osm_feature(key = 'waterway',
                                 value = c('river', 'stream'))

    streams_query$prefix <- sub('timeout:25', 'timeout:180', streams_query$prefix)

    streams <- osmdata::osmdata_sf(streams_query)
    streams <- streams$osm_lines$geometry

    streams_proj <- streams %>%
        sf::st_transform(crs = rast_crs) %>%
        sf::st_union() %>%
        # sf::st_transform(crs = WGS84) %>%
        sf::st_as_sf() %>%
        rename(geometry = x) %>%
        mutate(FID = 0:(n() - 1)) %>%
        dplyr::select(FID, geometry)

    if(! is.null(outfile)){

        sf::st_write(streams_proj,
                     dsn = outfile,
                     layer = 'streams',
                     driver = 'ESRI Shapefile',
                     delete_layer = TRUE,
                     quiet = TRUE)

        message(paste('OSM streams layer written to', outfile))

    } else {
        return(streams_proj)
    }
}

locate_test_results <- function(runid, strtg, loc){

    phase = ifelse(strtg == 'specialist', 'finetune', 'run')

    run_dirs <- list.files(loc)

    rundir <- grep(paste0('^', phase, runid, '_'), run_dirs, value = TRUE)
    if(length(rundir) > 1){

        warning('multiple ', phase, ' dirs found')
        return()

    } else if(length(rundir) == 0){

        warning('no ', phase, ' dirs found')
        return()
    }

    testdir <- file.path('test', list.files(file.path(loc, rundir, 'test')))
    if(length(testdir) != 1){
        warning('0 or multiple test dirs')
        return()
    }

    res <- file.path(loc, rundir, testdir, 'test_results.p')

    return(res)
}

retrieve_test_results <- function(runids, strategy, loc = nh_dir){

    nse_out <- kge_out <-
        matrix(NA_real_,
               nrow = length(site_codes),
               ncol = length(runids),
               dimnames = list(site_codes, NULL))

    for(i in seq_along(runids)){

        td <- try(locate_test_results(runids[i], strtg = strategy, loc = loc), silent = TRUE)
        if(inherits(td, 'try-error') || is.null(td)) next
        # write_lines(td, '/tmp/rundirs.txt', append = T)

        xx = reticulate::py_load_object(td)

        for(s in site_codes){
            try({
                nse_out[rownames(nse_out) == s, i] <- xx[[paste0(s, '_MANUALQ')]]$`1D`$discharge_NSE
            }, silent = TRUE)
            try({
                kge_out[rownames(kge_out) == s, i] <- xx[[paste0(s, '_MANUALQ')]]$`1D`$discharge_KGE
            }, silent = TRUE)
        }
        # if(all(is.na(nse_out[, i]))) stop('no nse', i)
        # if(all(is.na(kge_out[, i]))) stop('no kge', i)
    }

    return(list(nse = nse_out, kge = kge_out))
}

polygon_with_gaps <- function(df, border = NA){

    rl <- rle(is.na(df$Q_pred_int_2.5))
    vv <- ! rl$values
    chunkfac <- rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] <- 1
    chunks <- split(df, chunkfac)
    noNAchunks <- lapply(chunks, function(x) x[! is.na(x$Q_pred_int_2.5), ])

    for(i in 1:length(noNAchunks)){
        polygon(x = c(noNAchunks[[i]]$datetime, rev(noNAchunks[[i]]$datetime)),
                y = c(noNAchunks[[i]]$Q_pred_int_2.5, rev(noNAchunks[[i]]$Q_pred_int_97.5)),
                col = adjustcolor('red', alpha.f = 0.2), border = border)
    }
}

polygon_with_gaps2 <- function(df, gapcol, lowval, highval, col = 'blue',
                               invert = FALSE, border = col, hashed = FALSE){

    if(invert){
        rl = rle(is.na(df[[gapcol]]))
    } else {
        rl = rle(! is.na(df[[gapcol]]))
    }

    vv = ! rl$values
    chunkfac = rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] = 1
    chunks = split(df, chunkfac)

    if(invert){
        NAchunks = lapply(chunks, function(x) x[! is.na(x[[gapcol]]), ])
    } else {
        NAchunks = lapply(chunks, function(x) x[is.na(x[[gapcol]]), ])
    }

    NAchunks = Filter(function(x) difftime(x$datetime[nrow(x)], x$datetime[1], units = 'hours') >= 6,
                      NAchunks)

    if(length(NAchunks)){

        for(i in 1:length(NAchunks)){

            d <- NAchunks[[i]]

            if(hashed){
                polygon(x = c(d$datetime, rev(d$datetime)),
                        y = c(rep(lowval, nrow(d)), rep(highval, nrow(d))),
                        col = col, border = NA, angle = 45, density = 30)
            } else {
                polygon(x = c(d$datetime, rev(d$datetime)),
                        y = c(rep(lowval, nrow(d)), rep(highval, nrow(d))),
                        col = col, border = border)
            }
        }
    }
}

plot_gap_points <- function(d, dfill, dmask = NULL, mingap, col = 'black'){

    rl <- rle(! is.na(d$discharge))
    vv <- ! rl$values
    chunkfac <- rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] <- 1
    chunks <- split(d, chunkfac)
    NAchunks <- lapply(chunks, function(x) x[is.na(x$discharge), ])
    NAchunks <- Filter(function(x) difftime(x$datetime[nrow(x)], x$datetime[1], units = 'hours') >= mingap,
                       NAchunks)

    if(! is.null(dmask)){

        rl <- rle(! is.na(dmask$discharge))
        vv <- ! rl$values
        chunkfac <- rep(cumsum(vv), rl$lengths)
        chunkfac[chunkfac == 0] <- 1
        chunks <- split(dmask, chunkfac)
        maskchunks <- lapply(chunks, function(x) x[is.na(x$discharge), ])
        maskchunks <- Filter(function(x) difftime(x$datetime[nrow(x)], x$datetime[1], units = 'hours') >= mingap,
                             maskchunks)

        NAchunks <- Filter(function(x) any(sapply(maskchunks, function(y){
            snapdist <- if(s == 'CARI') 60 else 10
            abs(difftime(min(x$datetime), min(y$datetime))) < snapdist
        })), NAchunks)
    }

    if(length(NAchunks)){

        for(j in 1:length(NAchunks)){

            missing_dts <- NAchunks[[j]]$datetime[is.na(NAchunks[[j]]$discharge)]
            points(dfill$datetime[dfill$datetime %in% missing_dts],
                   dfill$Q_predicted[dfill$datetime %in% missing_dts],
                   col = col, pch = '.')
        }
    }
}

ts_plot <- function(site, yr, boldgray = FALSE, ymax = Inf, scaled = TRUE, border = NA){

    if(scaled){
        pred <- read_csv(paste0('out/lm_out_specQ/predictions/', site, '.csv')) %>%
            select(datetime, starts_with('Q_pred'))
        fit <- read_csv(paste0('out/lm_out_specQ/fit/', site, '.csv')) %>%
            select(datetime, Q_field)
    } else {
        pred <- read_csv(paste0('out/lm_out/predictions/', site, '.csv')) %>%
            select(datetime, starts_with('Q_pred'))
        fit <- read_csv(paste0('out/lm_out/fit/', site, '.csv')) %>%
            select(datetime, Q_field)
    }

    withf <- read_csv(glue('in/neon_continuous_Q_withflags/{site}.csv'))
    withoutf <- read_csv(glue('in/neon_continuous_Q/{site}.csv'))

    plotd <- full_join(pred, fit, by = 'datetime') %>%
        left_join(select(withf, datetime, qwith = discharge), by = 'datetime') %>%
        left_join(select(withoutf, datetime, qwo = discharge), by = 'datetime') %>%
        filter(datetime >= as.POSIXct(paste0(yr, '-01-01')),
               datetime <= as.POSIXct(paste0(yr, '-12-31'))) %>%
        arrange(datetime)

    flagdts <- plotd %>%
        filter(is.na(qwo) & ! is.na(qwith)) %>%
        pull(datetime)

    ymax_ <- max(c(plotd$qwith, plotd$Q_pred_int_97.5, plotd$Q_field), na.rm = TRUE)
    rle_ <- rle2(as.numeric(month(plotd$datetime)))

    if(nrow(rle_) < 12){
        rle_ <- rle_[-1, ]
    }

    xaxis_labs <- substr(month.abb, 1, 1)[rle_$values]

    ylm <- c(0, min(ymax_, ymax))
    plot(plotd$datetime, plotd$qwith, type = 'n',
         ylim = ylm, xlab = '', ylab = '', xaxt = 'n')
    axis(1, plotd$datetime[rle_$starts], xaxis_labs)
    axis(1, plotd$datetime[rle_$stops[nrow(rle_)]], yr, tick = F, font = 2)
    polygon_with_gaps(plotd, border = border)
    if(boldgray){
        lines(plotd$datetime, plotd$qwith, col = 'gray50', lwd = 3)
    } else {
        lines(plotd$datetime, plotd$qwith, col = 'gray50')
    }
    lines(plotd$datetime, plotd$Q_predicted, col = 'red')
    points(flagdts, rep(ylm[2] * -0.03, length(flagdts)), pch = 39, col = 'black')
    points(plotd$datetime, plotd$Q_field, col = 'black', pch = 1)
    mtext(site, 3, -1, adj = 0.01, padj = 0.1)
}

rle2 <- function(x){

    r <- rle(x)
    ends <- cumsum(r$lengths)

    r <- tibble(values = r$values,
                starts = as.integer(c(1, ends[-length(ends)] + 1)),
                stops = ends,
                lengths = r$lengths)

    return(r)
}

reduce_results <- function(res, name, f, ...){

    r_ <- suppressWarnings(apply(res, 1, f, ...))

    r_[is.infinite(r_)] <- NA

    out <- tibble(site = names(r_),
                  x = unname(r_)) %>%
        rename(!!sym(name) := x)

    return(out)
}

approxjoin_datetime <- function(x,
                                y,
                                rollmax = '7:30',
                                keep_datetimes_from = 'x',
                                indices_only = FALSE){

    #x and y: data.frames with datetime and val columns. may have others, but they will be dropped.
    #rollmax: the maximum snap duration for matching elements of x and y. Must
    #   be a string of the form "HH:MM:SS"
    #keep_datetimes_from: string. either 'x' or 'y'. the datetime column from
    #   the corresponding tibble will be kept, and the other will be dropped
    #indices_only: logical. if TRUE, a join is not performed. rather,
    #   the matching indices from each tibble are returned as a named list of vectors..

    #tests
    if(! keep_datetimes_from %in% c('x', 'y')) stop('keep_datetimes_from must be "x" or "y"')
    if(! 'datetime' %in% colnames(x) || ! 'datetime' %in% colnames(y)){
        stop('both x and y must have "datetime" columns containing POSIXct values')
    }
    if(! is.logical(indices_only)) stop('indices_only must be a logical')
    if(! grepl('[0-9]{2}:[0-9]{2}:[0-9]{2}', rollmax)){
        stop('rollmax must be a string of the form "HH:MM:SS"')
    }

    x <- x %>%
        rename_with(.fn = ~paste0(., '_x'),
                    .cols = everything()) %>%
        data.table::as.data.table()

    y <- y %>%
        rename_with(.fn = ~paste0(., '_y'),
                    .cols = everything()) %>%
        data.table::as.data.table()

    #convert the desired maximum roll distance from string to integer seconds
    h_m_s <- as.numeric(str_match(rollmax, '([0-9]{2}):([0-9]{2}):([0-9]{2})')[,2:4])
    h_ <- h_m_s[1] * 60 * 60
    m_ <- h_m_s[2] * 60
    s_ <- h_m_s[3]
    rollmax <- h_ + m_ + s_

    #create columns in x that represent the snapping window around each datetime
    x[, `:=` (datetime_min = datetime_x - rollmax,
              datetime_max = datetime_x + rollmax)]
    y[, `:=` (datetime_y_orig = datetime_y)] #datetime col will be dropped from y

    #join x rows to y if y's datetime falls within the x range
    joined <- y[x, on = .(datetime_y <= datetime_max,
                          datetime_y >= datetime_min)]
    joined <- na.omit(joined, cols = 'datetime_y_orig') #drop rows without matches

    #for any datetimes in x or y that were matched more than once, keep only
    #the nearest match
    joined[, `:=` (datetime_match_diff = abs(datetime_x - datetime_y_orig))]
    joined <- joined[, .SD[which.min(datetime_match_diff)], by = datetime_x]
    joined <- joined[, .SD[which.min(datetime_match_diff)], by = datetime_y_orig]

    if(indices_only){
        y_indices <- which(y$datetime_y %in% joined$datetime_y_orig)
        x_indices <- which(x$datetime_x %in% joined$datetime_x)
        return(list(x = x_indices, y = y_indices))
    }

    #drop and rename columns (data.table makes weird name modifications)
    if(keep_datetimes_from == 'x'){
        joined[, c('datetime_y', 'datetime_y.1', 'datetime_y_orig', 'datetime_match_diff') := NULL]
        data.table::setnames(joined, 'datetime_x', 'datetime')
    } else {
        joined[, c('datetime_x', 'datetime_y.1', 'datetime_y', 'datetime_match_diff') := NULL]
        data.table::setnames(joined, 'datetime_y_orig', 'datetime')
    }

    joined <- dplyr::select(joined, datetime, val_x, val_y)

    return(joined)
}

drop_var_prefix <- function(x){

    unprefixed <- substr(x, 4, nchar(x))

    return(unprefixed)
}

mode_interval_dt <- function(x){
    #mode interval in minutes for datetime vector
    Mode(diff(as.numeric(x)) / 60)
}

ms_write_netcdf <- function(df_list, path){

    vars_units0 <- c(discharge = 'cms', dayl = 's', prcp = 'mm/day',
                     srad = 'W/m2', swe = 'mm', tmax = 'C', tmin = 'C', vp = 'Pa',
                     pet = 'mm')

    for(i in seq_along(df_list)){

        d_chunk <- df_list[[i]]

        if(! 'pet' %in% colnames(d_chunk)){
            vars_units <- vars_units0[names(vars_units0) != 'pet']
        } else {
            vars_units <- vars_units0
        }

        ncdf_loc <- glue('{pth}/{s}.nc',
                         pth = path,
                         s = d_chunk$site_code[1])

        dim_time <- ncdim_def('date',
                              units = 'days since 1970-01-01 00:00',
                              calendar = 'standard',
                              vals = d_chunk$date)

        ncdf_varlist <- list()
        for(j in seq_along(vars_units)){

            ncdf_varlist[[j]] <- ncvar_def(name = names(vars_units)[j],
                                           units = unname(vars_units)[j],
                                           dim = list(dim_time),
                                           missval = NULL)
        }

        con <- nc_create(filename = ncdf_loc,
                         vars = ncdf_varlist)

        for(j in seq_along(vars_units)){

            varname <- ncdf_varlist[[j]]$name
            ncvar_put(con, varname, d_chunk[[varname]])
        }

        nc_close(con)
    }
}

restore_transient <- function(x, orig_head, orig_tail, trans){

    #restore leading
    x[1:trans, ] <- left_join(select(x[1:trans, ], datetime),
                              orig_head,
                              by = 'datetime') %>%
        relocate(datetime, .after = 'discharge')

    #restore trailing
    trail_trans <- (nrow(x) - (trans - 1)):nrow(x)
    x[trail_trans, ] <- left_join(select(x[trail_trans, ], datetime),
                                  orig_tail,
                                  by = 'datetime') %>%
        relocate(datetime, .after = 'discharge')

    return(x)
}

get_populatable_indices <- function(src_column, mingap){

    popul_inds <- rle2(is.na(src_column)) %>%
        filter(values, lengths >= mingap) %>%
        {mapply(function(x, y) x:y, .$starts, .$stops)} %>%
        unlist()

    popul_inds <- seq_along(src_column) %in% popul_inds #as logical

    return(popul_inds)
}

prepare_q_neon <- function(site, smooth_plot = FALSE){

    window_size <- 15 #keep it odd
    q <- read_csv(glue('in/neon_continuous_Q/{site}.csv'))

    if(site != 'TOMB'){

        transient <- window_size %/% 2
        q_bounds <- select(q, -site_code, -discharge)
        q_orig <- select(q, datetime, discharge_orig = discharge)
        q <- select(q, datetime, discharge)

        #fill out missing data (make explicit) to one minute interval
        q <- complete(q, datetime = seq(min(datetime), max(datetime), by = '1 min'))

        #pad head and tail for moving average window
        buffer1 <- tibble(datetime = seq(q$datetime[1] - (60 * transient),
                                         q$datetime[1] - 60,
                                         by = '1 min'))
        buffer2 <- tibble(datetime = seq(q$datetime[nrow(q)] + 60,
                                         q$datetime[nrow(q)] + (60 * transient),
                                         by = '1 min'))
        q <- bind_rows(buffer1, q, buffer2)

        q_dt <- q$datetime
        q_na <- is.na(q$discharge)

        #traingular moving average to smooth over bayesian obs error
        q <- xts(q$discharge, q$datetime, tzone = 'UTC')
        q <- rollmean(q, window_size, fill = NA, align = 'center', na.rm = TRUE)
        q <- rollmean(q, window_size, fill = NA, align = 'center', na.rm = TRUE)
        q <- suppressWarnings(as.vector(q))

        #restore missing values
        q[q_na] <- NA_real_

        #need to interpolate uncertainty too, so it isn't lost during downsampling
        q <- left_join(tibble(discharge = q, datetime = q_dt), q_bounds, by = 'datetime') %>%
            slice((transient + 1):(length(q) - transient)) %>%
            mutate(across(-datetime, ~na_interpolation(., , maxgap = 14)))
        # mutate(across(-datetime, ~na_seadec(ts(., deltat = 1/1440), maxgap = 14)))

        #downsample to even 5 minute interval. priority breakdown:
        #   1. preserve observations originally on minutes divisible by 5
        #   2. accept newly interpolated values where the original interval was
        #      2-15 minutes and/or offset from a natural 5-minute sequence
        #   3. shift original values by up to 2 minutes if interpolation was
        #      not performed, e.g. where the sampling interval is 16-60 minutes
        q$m5 <- FALSE
        q$m5[minute(q$datetime) %% 5 == 0] <- TRUE
        q$datetime <- round_date(q$datetime, unit = '5 min')
        q <- group_by(q, datetime) %>%
            slice(if(any(m5)) which(m5)[1] else if(any(! is.na(discharge))) which(! is.na(discharge))[1] else 1) %>%
            ungroup() %>%
            select(-m5)

        if(any(duplicated(q$datetime))) warning('dupes in ', site)

        if(smooth_plot){

            require(dygraphs)
            require(htmlwidgets)

            zz = full_join(q, q_orig, by = 'datetime')
            dygraphs::dygraph(xts(x = select(zz, discharge_orig, discharge_lower, discharge_upper, discharge),
                                  order.by = zz$datetime)) %>%
                dyRangeSelector() %>%
                saveWidget(glue('figs/smooth_plots/{site}.html'))
        }
    }

    q <- q %>%
        mutate(source = 'NEON') %>%
        select(datetime, discharge_Ls = discharge, discharge_lower95_Ls = discharge_lower,
               discharge_upper95_Ls = discharge_upper, source)

    return(q)
}

prepare_q_lm <- function(site){

    lm_res <- read_csv('out/lm_out/results.csv') %>%
        filter(site_code == site) %>%
        pull(kge)

    lm_sq_res <- read_csv('out/lm_out_specQ/results_specificq.csv') %>%
        filter(site_code == site) %>%
        pull(kge)

    if(is.na(lm_res)) stop('no lm result for this site')

    strategy <- ifelse(lm_res > lm_sq_res, 'lm_out', 'lm_out_specQ')

    q <- read_csv(glue('out/{strategy}/predictions/{site}.csv'))

    if('date' %in% colnames(q)){
        q$datetime <- as_datetime(paste(q$date, '12:00:00'))
    }

    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5)

    # plot(q$datetime, q$discharge_Ls, type = 'l', ylim = c(0, 1000))
    # points(q$datetime[!q$m5], q$discharge_Ls[!q$m5], col = 'red')

    q <- complete(q, datetime = seq(min(datetime), max(datetime), by = '5 min'))
    q$m5 <- FALSE
    q$m5[minute(q$datetime) %% 5 == 0] <- TRUE
    q$datetime <- round_date(q$datetime, unit = '5 min')
    q <- group_by(q, datetime) %>%
        slice(if(any(m5)) which(m5)[1] else which(! is.na(discharge_Ls))[1]) %>%
        ungroup() %>%
        select(-m5)

    if(any(duplicated(q$datetime))) warning('dupes in ', site)

    # q$discharge_Ls <- na_seadec(ts(q$discharge_Ls, deltat = 1/288), maxgap = 2)
    q$discharge_Ls <- na_interpolation(q$discharge_Ls, maxgap = 2)

    q$source <- ifelse(strategy == 'lm_out', 'Linreg', 'Linreg_scaled')

    return(q)
}

prepare_q_lstm <- function(site){

    q <- read_csv(glue('out/lstm_out/predictions/{site}.csv'))

    q$source <- filter(results, site_code == !!site)$strategy
    q$datetime <- as_datetime(paste(q$date, '12:00:00'))

    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5, source)

    return(q)
}
