# Comparison of performance of risk stratification tools at predicting
# overall survival

  insert_head()
  
# container ------
  
  surv_tools <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## analysis variables
  
  surv_tools$var_lexicon <- 
    rbind(tibble(variable = c('score', 'ensemble_score', 'RF ensemble'), 
                 label = c('ElasticNet', 'LASSO ensemble', 'RF ensemble')), 
          pah_study$comparators[c('variable', 'label')])
  
  ## analysis tables
  
  surv_tools$analysis_tbl <- 
    map2(multi_cox$lp_scores, 
         pah_study$data_master %>% 
           blast(timepoint) %>% 
           map(select, 
               ID, any_of(surv_tools$var_lexicon$variable)), 
         left_join, by = 'ID') %>% 
    map2(lasso_tools$lp_scores %>% 
           map(~.x[c('ID', 'ensemble_score')]), 
         left_join, by = 'ID') %>% 
    map(column_to_rownames, 'ID')
  
  surv_tools$analysis_tbl$LZ_0 <-  surv_tools$analysis_tbl$LZ_0 %>% 
    select(-starts_with('Reveal'))
  
  ## cohort specific variables (no REVEAL for the LZ/W cohort!)

  surv_tools$variables <- surv_tools$analysis_tbl %>% 
    map(select, any_of(surv_tools$var_lexicon$variable)) %>% 
    map(names)
  
  ## model formulas
  
  surv_tools$formulas <- surv_tools$variables %>% 
    map(~paste('Surv(surv_months, death_study) ~', .x) %>% 
          map(as.formula)) %>% 
    map2(., surv_tools$variables, 
         set_names)
    
# Cox models ------
  
  insert_msg('Cox models')
  
  ## working with a call to keep the entire formula
  ## with a direct reference to the dataset variables
  ## in the Cox model
  
  surv_tools$models <- 
    map2(surv_tools$formulas, 
         surv_tools$analysis_tbl, 
         function(form, data) form %>% 
           map(~call2(.fn = 'coxph', 
                      formula = .x, 
                      data = data, 
                      x = TRUE))) %>% 
    map(~map(.x, eval)) %>% 
    map2(., surv_tools$analysis_tbl, 
         function(mod, data) mod %>% 
           map(~as_coxex(cox_model = .x, 
                         data = data)))
  
# Model assumptions ------
  
  insert_msg('Model assumptions')
  
  ## proportional hazard assumption is the key
  
  surv_tools$assumptions <- surv_tools$models %>% 
    map(~map(.x, summary, 'assumptions')) %>% 
    map(compress, names_to = 'variable')
  
# Performance statistics ------
  
  insert_msg('Model performance')

  surv_tools$fit_stats <- surv_tools$models %>% 
    map(~map(.x, summary, 'fit')) %>% 
    map(compress, names_to = 'variable') %>% 
    compress(names_to = 'cohort') %>% 
    full_rbind(rf_tools$fit_stats %>% 
                 mutate(variable = 'RF ensemble')) %>% 
    blast(cohort)

# Plotting C-index against IBS ------
  
  insert_msg('R-squared and IBS scatter plots')
  
  ## IBS: integrated Brier score
  
  surv_tools$fit_plots <- 
    list(x = surv_tools$fit_stats, 
         y = c('Training: IBK', 'Test: LZ/W'), 
         z = globals$center_colors[names(surv_tools$analysis_tbl)]) %>% 
    pmap(function(x, y, z) x %>% 
           ggplot(aes(x = 1 - ibs_model, 
                      y = c_index)) + 
           geom_point(shape = 21, 
                      size = 2, 
                      fill = z) + 
           geom_text_repel(aes(label = exchange(variable, 
                                                dict = surv_tools$var_lexicon)), 
                           size = 2.75) + 
           globals$common_theme + 
           labs(title = y, 
                x = '1 - IBS', 
                y = 'C-index'))
  
# Assessing calibration of the score and established tools -------
  
  insert_msg('Assessing calibration of the score and the tools')
  
  ## calibrator objects
  
  surv_tools$calibrator_obj$IBK_0 <-
    list(fit = surv_tools$models$IBK_0, 
         use_unique = c(FALSE, FALSE, 
                        rep(TRUE, 7)), 
         labels = list(score = paste0('T', 1:3), 
                       ensemble_score = paste0('T', 1:3), 
                       mRasp = c('low', 'int', 'high'), 
                       COMPERA = c('low', 'int', 'high'), 
                       SPAHR = c('low', 'int', 'high'), 
                       FRENCH3p = as.character(0:3), 
                       FRENCH4p = as.character(0:4), 
                       Reveal_lite2_3_cat = c('low', 'int', 'high'), 
                       Reveal2_risk_3_cat = c('low', 'int', 'high'))) %>% 
    pmap(calibrate)
  
  surv_tools$calibrator_obj$LZ_0 <-
    list(fit = surv_tools$models$LZ_0, 
         use_unique = c(FALSE, FALSE, 
                        rep(TRUE, 5)), 
         labels = list(score = paste0('T', 1:3), 
                       ensemble_score = paste0('T', 1:3), 
                       mRasp = c('low', 'int', 'high'), 
                       COMPERA = c('low', 'int', 'high'), 
                       SPAHR = c('low', 'int', 'high'), 
                       FRENCH3p = as.character(0:3), 
                       FRENCH4p = as.character(0:4))) %>% 
    pmap(calibrate)
  
  ## global calibration stats
  
  surv_tools$calibration_stats <- surv_tools$calibrator_obj %>% 
    map(~map(.x, summary)) %>% 
    map(compress, names_to = 'variable')
  
  ## differences in survival between the strata
  
  surv_tools$calibration_test <- surv_tools$calibrator_obj %>% 
    map(~map(.x, ~.x$surv_fit) %>% 
          map(surv_pvalue, 
              method = 'survdiff')) %>% 
    map(compress, names_to = 'variable') %>% 
    map(mutate, 
        p_adjusted = p.adjust(pval, 'BH'), 
        significance = ifelse(p_adjusted >= 0.05, 
                              paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                              ifelse(p_adjusted < 0.01, 
                                     'p < 0.001', 
                                     paste('p =', signif(p_adjusted, 2)))))
    
  ## plots of survival within the Elastic Net score strata
  ## and risk groups defined by the remaining scores
  
  surv_tools$calibration_plots$IBK_0 <-
    list(x = surv_tools$calibrator_obj$IBK_0, 
         title = names(surv_tools$calibrator_obj$IBK_0) %>% 
           exchange(dict = surv_tools$var_lexicon) %>% 
           paste('Training: IBK,', .), 
         legend.title = c('Tertile', 
                          'Tertile', 
                          rep('Risk strata', 3), 
                          rep('# risk factors', 2), 
                          rep('Risk strata', 2)), 
         palette = list(c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue4', 
                          'steelblue2', 
                          'coral2'),
                        c('darkolivegreen4', 
                          'steelblue4', 
                          'steelblue2', 
                          'coral2', 
                          'coral4'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'))) %>% 
    pmap(plot,
         xlab = 'Overall survival, months', 
         show_cox = FALSE) %>% 
    map(~.x + globals$common_theme)
  
  surv_tools$calibration_plots$LZ_0 <-
    list(x = surv_tools$calibrator_obj$LZ_0, 
         title = names(surv_tools$calibrator_obj$LZ_0) %>% 
           exchange(dict = surv_tools$var_lexicon) %>% 
           paste('Test: LZ/W,', .), 
         legend.title = c('Tertile', 
                          'Tertile', 
                          rep('Risk strata', 3), 
                          rep('# risk factors', 2)), 
         palette = list(c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue', 
                          'coral3'), 
                        c('darkolivegreen4', 
                          'steelblue4', 
                          'steelblue2', 
                          'coral2'),
                        c('darkolivegreen4', 
                          'steelblue4', 
                          'steelblue2', 
                          'coral2', 
                          'coral4'))) %>% 
    pmap(plot,
         xlab = 'Overall survival, months', 
         show_cox = FALSE) %>% 
    map(~.x + globals$common_theme)
  
  ## annotating the KM plots with the log-rank test results
  
  surv_tools$calibration_plots <- 
    map2(surv_tools$calibration_plots, 
         surv_tools$calibration_test, 
         function(plot, stat) map2(plot, stat$significance, 
                                   ~.x + 
                                     annotate('text', 
                                              label = .y, 
                                              size = 2.75, 
                                              x = 10, 
                                              y = 0.1,
                                              hjust = 0, 
                                              vjust = 0)))
  
# END ------
  
  insert_tail()