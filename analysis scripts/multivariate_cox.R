# Univariate cox modeling. Numeric variables are splinned with rms::rcs

  insert_head()
  
# container list -----
  
  multi_cox <- list()
  
# globals: median-centered analysis table -----
  
  insert_msg('Analysis table and folds')
  
  ## first order variables
  
  multi_cox$variables <- pah_study$mod_variables$variable
  
  ## analysis table: first and second order terms, median centering
  
  multi_cox$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, -surv_months, -death_study) %>% 
    map(~map_dfc(.x, function(var) if(is.numeric(var)) var^2 else NULL))
  
  multi_cox$analysis_tbl <- multi_cox$analysis_tbl %>% 
    map(set_names, paste0(names(multi_cox$analysis_tbl[[1]]), '_sec')) %>% 
    map2(., pah_study[c('IBK_0', 'LZ_0')], cbind) %>% 
    map(as_tibble) %>% 
    map(~map_dfc(.x, function(var) if(is.numeric(var)) scale(var, center = median(var))[, 1] else var)) %>% 
    map(select, -surv_months, -death_study, -event3, -event5, -death_study_fct) %>% 
    map(column_to_rownames, 'ID')

  ## folds
  
  set.seed(1234)
  
  multi_cox$folds <- 1:200 %>% 
    map(~createFolds(y = multi_cox$analysis_tbl$IBK_0[, 1], k = 10, list = FALSE, returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:200))
  
  ## n numbers
  
  multi_cox$n_numbers <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(count, death_study) %>% 
    map(~c(total = sum(.x$n), 
           events = .x$n[2]))
  
# finding the optimal lambda values by 100-repetition, 10-fold CV, IBK cohort ------
  
  insert_msg('Lambda tuning')
  
  plan('multisession')
  
  multi_cox$lambda_tune <- multi_cox$folds %>% 
    future_map(~cv.glmnet(x = model.matrix(~., multi_cox$analysis_tbl$IBK_0), 
                          y = pah_study$surv_obj$IBK_0, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = 0.5))
  
  plan('sequential')
  
  multi_cox$lambda_tbl <- multi_cox$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')])) %>% 
    map2_dfr(., multi_cox$lambda_tune, ~filter(.x, lambda == .y[['lambda.min']]))
  
  multi_cox$opt_lambda <- multi_cox$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# fitting the training model, linear predictor scores for the training and the test cohort -----
  
  insert_msg('Fitting the training LASSO model')
  
  multi_cox$model_IBK_0 <- glmnet(x = model.matrix(~., multi_cox$analysis_tbl$IBK_0), 
                                  y = pah_study$surv_obj$IBK_0, 
                                  family = 'cox', 
                                  alpha = 0.5, 
                                  lambda = multi_cox$opt_lambda$lambda)
  
  ## linear predictor scores
  
  multi_cox[c('score_IBK_0', 'score_LZ_0')] <- multi_cox$analysis_tbl %>% 
    map(~model.matrix(~., .x)) %>% 
    map(~predict(multi_cox$model_IBK_0, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'ID') %>% 
    map2(., map(pah_study[c('IBK_0', 'LZ_0')], ~.x[c('ID', 'surv_months', 'death_study')]), 
         left_join, by = 'ID') %>% 
    map(set_names, c('ID', 'score', 'surv_months', 'death_study')) %>% 
    map(mutate, ## score tertiles for plotting
        score_quart = cut(score, 
                          c(-Inf, quantile(score, c(0.33, 0.66)), Inf), 
                          c('T1', 'T2', 'T3')))
  
  ## univariate coxph and cph models
  
  multi_cox[c('coxph_IBK_0', 'coxph_LZ_0')] <- multi_cox[c('score_IBK_0', 'score_LZ_0')] %>% 
    map(~model_spline_cox(data = .x, 
                          event_variable = 'death_study', 
                          time_variable = 'surv_months', 
                          indep_variable = 'score', 
                          sec_order = FALSE))
    
  multi_cox[c('cph_IBK_0', 'cph_LZ_0')] <- multi_cox[c('score_IBK_0', 'score_LZ_0')] %>% 
    map(~cph(Surv(surv_months, death_study) ~ score, data = .x, x = TRUE, y = TRUE, time.inc = 12, surv = TRUE))
  
# Fit statistics of the Cox models based on the linear predictor scores -----
  
  insert_msg('Model fit stats')
  
  multi_cox$coxph_stats <- multi_cox[c('coxph_IBK_0', 'coxph_LZ_0')] %>% 
    map2_dfr(., c('IBK_0', 'LZ_0'), ~mutate(.x$summary, cohort = .y))
  
  multi_cox$cph_stats <- multi_cox[c('cph_IBK_0', 'cph_LZ_0')] %>% 
    map(validate) %>% 
    map(unclass) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'statistic') %>% 
    as_tibble %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
# Model calibration -----
  
  insert_msg('Model calibration check')
  
  multi_cox$calibration <- multi_cox[c('cph_IBK_0', 'cph_LZ_0')] %>% 
    map(function(cohort) c(12, 24, 36, 48) %>% 
          map(~calibrate(fit = cohort, u = .x)) %>% 
          set_names(c('mo12', 'mo24', 'mo36', 'mo48'))) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
# END ----
  
  insert_tail()