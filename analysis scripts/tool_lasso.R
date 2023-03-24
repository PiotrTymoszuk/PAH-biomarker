# LASSO ensemble of the established risk assessment tools
# and of the newly developed Elastic Net score

  insert_head()

# container ------

  lasso_tools <- list()

# parallel backend -------

  insert_msg('Parallel backend')
  
  plan('multisession')

# analysis globals ------

  insert_msg('Analysis globals')

  ## analysis variables
  
  lasso_tools$var_lexicon <- 
    rbind(tibble(variable = 'score', 
                 label = 'ElasticNet'), 
          pah_study$comparators[c('variable', 'label')]) %>% 
    filter(!stri_detect(variable, fixed = 'Reveal'))

# Analysis tables -------
  
  insert_msg('Analysis tables')
  
  ## with established risk tools as factors
  
  lasso_tools$analysis_tbl <- pah_study$data_master %>% 
    select(ID, timepoint, any_of(lasso_tools$var_lexicon$variable)) %>% 
    map_dfc(function(x) if(is.numeric(x)) factor(x) else x) %>% 
    dlply('timepoint', select, -timepoint)
  
  ## adding the linear predictor scores of the Elastic Net
  
  lasso_tools$analysis_tbl <- 
    map2(multi_cox$lp_scores %>% 
           map(~.x[c('ID', 'score')]), 
         lasso_tools$analysis_tbl, 
         left_join, by = 'ID') %>% 
    map(column_to_rownames, 'ID')

# CV folds ------
  
  insert_msg('CV folds')
  
  set.seed(1234)
  
  lasso_tools$folds <- 1:200 %>% 
    map(~createFolds(y = lasso_tools$analysis_tbl$IBK_0[, 1], 
                     k = 10, 
                     list = FALSE, 
                     returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:200))
  
# finding the optimal lambda values by 100-repetition, 10-fold CV, IBK cohort ------
  
  insert_msg('Lambda tuning')
  
  lasso_tools$lambda_tune <- lasso_tools$folds %>% 
    future_map(~cv.glmnet(x = model.matrix(~., lasso_tools$analysis_tbl$IBK_0), 
                          y = pah_study$surv_obj$IBK_0, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = 1), 
               .options = furrr_options(seed = TRUE))
  
  lasso_tools$lambda_tbl <- lasso_tools$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo')])) %>% 
    map2_dfr(., lasso_tools$lambda_tune, 
             ~filter(.x, lambda == .y[['lambda.min']]))
  
  lasso_tools$opt_lambda <- lasso_tools$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# fitting the training model, linear predictor scores for the training and the test cohort -----
  
  insert_msg('Fitting the training Elastic Net model')
  
  lasso_tools$glmnet_IBK_0 <- 
    glmnet(x = model.matrix(~., lasso_tools$analysis_tbl$IBK_0), 
           y = pah_study$surv_obj$IBK_0, 
           family = 'cox', 
           alpha = 1, 
           lambda = lasso_tools$opt_lambda$lambda)
  
  ## linear predictor scores
  
  lasso_tools$lp_scores <- lasso_tools$analysis_tbl %>% 
    map(~model.matrix(~., .x)) %>% 
    map(~predict(lasso_tools$glmnet_IBK_0, newx = .x)) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'ID') %>% 
    map2(., 
         map(pah_study[c('IBK_0', 'LZ_0')], 
             ~.x[c('ID', 'surv_months', 'death_study')]), 
         left_join, by = 'ID') %>% 
    map(set_names, c('ID', 'ensemble_score', 'surv_months', 'death_study'))
  
  ## univariate coxph and cph models
  
  lasso_tools$coxph_models <- lasso_tools$lp_scores %>% 
    map(~model_spline_cox(data = .x, 
                          event_variable = 'death_study', 
                          time_variable = 'surv_months', 
                          indep_variable = 'ensemble_score', 
                          sec_order = FALSE))
  
  lasso_tools$cph_models  <- lasso_tools$lp_scores %>% 
    map(~cph(Surv(surv_months, death_study) ~ ensemble_score, 
             data = .x, 
             x = TRUE, 
             y = TRUE, 
             time.inc = 12, 
             surv = TRUE))
  
# Fit statistics of the Cox models based on the linear predictor scores -----
  
  insert_msg('Model fit stats')
  
  lasso_tools$coxph_stats <- lasso_tools$coxph_models %>% 
    map(~.x$summary) %>% 
    compress(names_to = 'cohort')
  
  lasso_tools$cph_stats <- lasso_tools$cph_models %>% 
    map(validate) %>% 
    map(unclass) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'statistic') %>% 
    compress(names_to = 'cohort') %>% 
    as_tibble
  
# Model calibration with RMS -----
  
  insert_msg('Model calibration check with RMS')
  
  lasso_tools$rms_calibration <- lasso_tools$cph_models %>% 
    map(function(cohort) c(12, 24, 36, 48) %>% 
          map(~calibrate(fit = cohort, u = .x)) %>% 
          set_names(c('mo12', 'mo24', 'mo36', 'mo48')))
  
# Model calibration, D'Agustino - Nam ---------
  
  insert_msg('Calibration with the DAgustino - Nam method')
  
  lasso_tools$dn_calibration <- lasso_tools$coxph_models %>% 
    map(~.x$model) %>% 
    map(get_cox_calibration, 
        n = 3, 
        labels = c('T1', 'T2', 'T3'))
  
# Model estimates -------
  
  insert_msg('Model estimates')
  
  ## variable extraction regex
  
  lasso_tools$var_regex <- lasso_tools$var_lexicon$variable %>% 
    sort(decreasing = TRUE) %>% 
    paste(collapse = '|')
  
  ## GLMNET model coefficients
  
  lasso_tools$coefs <- lasso_tools$glmnet_IBK_0 %>% 
    coef %>% 
    as.matrix %>% 
    as.data.frame %>% 
    filter(s0 != 0) %>% 
    rownames_to_column('parameter') %>% 
    mutate(variable = stri_extract(parameter, 
                                   regex = lasso_tools$var_regex), 
           level = stri_replace(parameter, 
                                regex =  lasso_tools$var_regex, 
                                replacement = ''), 
           regulation = ifelse(s0 < 0, 
                               'negative', 'positive'), 
           hr = exp(s0), 
           y_txt = ifelse(!is.na(level), 
                          paste(exchange(variable, 
                                         dict = lasso_tools$var_lexicon), 
                                level, sep = ': '), 
                          exchange(variable,
                                   dict = lasso_tools$var_lexicon)))
  
  ## Forest plot
  
  lasso_tools$coef_plot <- lasso_tools$coefs %>% 
    mutate(variable = reorder(variable, -s0)) %>% 
    ggplot(aes(x = hr, 
               y = reorder(y_txt, abs(s0)), 
               color = regulation)) + 
    geom_vline(xintercept = 1, 
               linetype = 'dashed') + 
    geom_segment(aes(x = 1, 
                     xend = hr, 
                     y = reorder(y_txt, hr), 
                     yend = reorder(y_txt, hr))) + 
    geom_point(size = 2, 
               shape = 16) + 
    geom_text(aes(label = signif(hr, 3)), 
              size = 2.75, 
              hjust = 0.5, 
              vjust = -0.9) + 
    scale_color_manual(values = c('negative' = 'steelblue', 
                                  'positive' = 'coral3'), 
                       labels = c('negative' = 'favorable', 
                                  'positive' = 'unfavorable'), 
                       name = '') + 
    scale_x_continuous(trans = 'log', 
                       labels = function(x) signif(x, 2)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    labs(title = 'LASSO risk assessment tool ensemble', 
         plot_tag = paste0('Total: n = ', 
                           multi_cox$n_numbers$IBK_0['total'], 
                           ', Events: n = ', 
                           multi_cox$n_numbers$LZ_0['events']), 
         x = expression('HR'[LASSO]))
  
# END ------
  
  plan('sequential')
  
  insert_tail()
  
