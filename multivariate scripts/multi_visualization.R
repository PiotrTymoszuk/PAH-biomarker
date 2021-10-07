# This script generates plots with the results of multivariate modeling 

  insert_head()
  
# data container -----
  
  multi_plots <- list()
  
# analysis and plotting tables -----
  
  insert_msg('Analysis tables')
  
  ## all stats
  
  multi_plots$model_stats <- list(data = list(multi_modeling$cv_pass$stats, 
                                              multi_modeling$cv_pass$cv_results, 
                                              multi_modeling$validation$stats$LZ_0), 
                                  suffix = c('IBK_0', 
                                             'cv', 
                                             'LZ_0')) %>% 
    pmap(function(data, suffix) mutate(data, 
                                       !!sym(paste0('c_', suffix)) := c_index, 
                                       !!sym(paste0('lower_', suffix)) := lower_ci, 
                                       !!sym(paste0('upper_', suffix)) := upper_ci, 
                                       !!sym(paste0('rsq_', suffix)) := rsq, 
                                       !!sym(paste0('mse_', suffix)) := mse) %>% 
           select(model_id, 
                  ends_with(suffix))) %>% 
    reduce(left_join, 
           by = 'model_id') %>% 
    mutate(model_lab = stri_replace(model_id, fixed = 'sign_', replacement = 'Sign '))
  
  ## C indexes, appending with the comparator C index values

  multi_plots$model_c <- multi_plots$model_stats %>% 
    select(model_id, 
           starts_with('c_'), 
           starts_with('lower'), 
           starts_with('upper'))
  
  multi_plots$model_c <- list(prefix = c('c_', 'lower_', 'upper_'), 
                              name = c('c_index', 'lower_ci', 'upper_ci')) %>% 
    pmap(function(prefix, name) multi_plots$model_c %>% 
           select(model_id, 
                  starts_with(prefix)) %>% 
           gather(key = 'dataset', 
                  value = 'value', 
                  starts_with(prefix)) %>% 
           mutate(dataset = stri_extract(dataset, regex = 'IBK_0|cv|LZ_0')) %>% 
           set_names(c('model_id', 'dataset', name))) %>% 
    reduce(left_join, 
           by = c('model_id', 'dataset')) %>% 
    mutate(sign_type = 'new')
  
  multi_plots$model_c <- multi_modeling$comparators$stats %>% 
    map2(., c('IBK_0', 'LZ_0'), 
         ~mutate(.x, dataset = .y, sign_type = 'comparator')) %>% 
    map_dfr(select, 
            model_id, 
            dataset, 
            c_index, 
            lower_ci, 
            upper_ci, 
            sign_type) %>% 
    rbind(multi_plots$model_c, .) %>% 
    mutate(model_lab = ifelse(sign_type == 'comparator', 
                              globals$comp_labs[model_id], 
                              stri_replace(model_id, fixed = 'sign_', replacement = 'Sign ')), 
           dataset = factor(dataset, c('IBK_0', 'LZ_0', 'cv')), 
           sign_type = factor(sign_type, c('new', 'comparator')))
  
  ## 5-year mortality ROC: cutpoint stats
  
  multi_plots$model_cutstats <- list(data = five_surv$roc_modeling$stats, 
                                        suffix = c('IBK_0', 'LZ_0')) %>% 
    pmap(function(data, suffix) data %>% 
           mutate(!!sym(paste0('Se_', suffix)) := Se, 
                  !!sym(paste0('Sp_', suffix)) := Sp, 
                  !!sym(paste0('J_', suffix)) := Optimal.criterion) %>%
           select(model_id, 
                  ends_with(suffix))) %>% 
    reduce(left_join, by = 'model_id') %>% 
    mutate(model_lab = ifelse(model_id %in% pah_study$comparators$variable, 
                              globals$comp_labs[model_id], 
                              stri_replace(model_id, fixed = 'sign_', replacement = 'Sign ')), 
           sign_type = ifelse(model_id %in% pah_study$comparators$variable, 
                              'comparator', 'new'))
    
  ## 5-year mortality ROC: AUC values
  
  multi_plots$model_auc <- five_surv$roc_modeling$stats %>%
    map2_dfr(., names(.), ~mutate(.x, dataset = .y)) %>% 
    mutate(lower_ci = lowerCI, 
           upper_ci = upperCI, 
           sign_type = ifelse(model_id %in% pah_study$comparators$variable, 
                              'comparator', 'new'), 
           model_lab = ifelse(sign_type == 'comparator', 
                              globals$comp_labs[model_id], 
                              stri_replace(model_id, fixed = 'sign_', replacement = 'Sign ')), 
           dataset = factor(dataset, c('IBK_0', 'LZ_0', 'cv')), 
           sign_type = factor(sign_type, c('new', 'comparator'))) %>% 
    select(model_id, 
           dataset, 
           model_lab, 
           sign_type, 
           AUC, 
           lower_ci, 
           upper_ci)
  
# Plotting the Rsq and MSE values -----
  
  insert_msg('Plots of the model Rsq and MSE')
  
  multi_plots$r_error_plots <- list(x_var = c('rsq_IBK_0', 'mse_IBK_0', 'rsq_IBK_0', 'mse_IBK_0'), 
                                    y_var = c('rsq_cv', 'mse_cv', 'rsq_LZ_0', 'mse_LZ_0'), 
                                    plot_title = c(rep('Risk signature development: OS', 2), 
                                                   rep('Risk signature testing: OS', 2)), 
                                    plot_subtitle = c(rep('Training IBK cohort and cross-validation', 2), 
                                                      rep('Training IBK and test LZ/W cohort', 2)), 
                                    x_lab = list(expression('R'^2*', IBK'), 
                                                 'MSE IBK', 
                                                 expression('R'^2*', IBK'), 
                                                 'MSE IBK'), 
                                    y_lab = list(expression('R'^2*', CV'), 
                                                 'MSE CV', 
                                                 expression('R'^2*', LZ/W'), 
                                                 'MSE LZ/W')) %>% 
    pmap(plot_model_stats, 
         inp_tbl = multi_plots$model_stats, 
         highlight = multi_plots$model_stats$model_id, 
         label_var = 'model_lab', 
         default_fill = globals$signature_color['new']) %>% 
    set_names(c('development_rsq', 
                'development_mse', 
                'testing_rsq', 
                'testing_mse'))
  
# Plotting the C indexes with confidence intervals -----
  
  insert_msg('C index plots')
  
  ## base plots
  
  multi_plots$c_plots <- list(inp_tbl = list(multi_plots$model_c %>% 
                                               filter(dataset %in% c('IBK_0', 'cv'), 
                                                      sign_type == 'new'), 
                                             multi_plots$model_c %>% 
                                               filter(dataset %in% c('IBK_0', 'LZ_0'))), 
                              plot_title = c('Risk signature development: OS', 
                                             'Risk signature testing: OS'), 
                              plot_subtitle = c('Training IBK cohort and cross-validation', 
                                                'Training IBK and test LZ/W cohort')) %>% 
    pmap(plot_c) %>% 
    set_names(c('development_c', 
                'testing_c'))
  
  ## splitting by the cohort/comparator
  
  multi_plots$c_plots$development_c <- multi_plots$c_plots$development_c + 
    facet_grid(. ~ dataset, 
               labeller = as_labeller(globals$center_labs)) + 
    guides(color = F)
  
  multi_plots$c_plots$testing_c <- multi_plots$c_plots$testing_c + 
    facet_grid(sign_type ~ dataset, 
               labeller = labeller(.rows = globals$signature_labels, 
                                   .cols = globals$center_labs), 
               space = 'free', 
               scales = 'free_y') + 
    guides(color = F)
  
# Plotting the Forest plots for particular multivariable risk signature Cox models in the IBK cohort ----
  # Extracting the score formula texts
  
  insert_msg('Forest plots with the signature Cox model estimates, score formulas')
  
  multi_plots$signature_forests <- list(cox_model = multi_modeling$cv_pass$models, 
                                        plot_title = names(multi_modeling$cv_pass$models) %>% 
                                          stri_replace(fixed = 'sign_', replacement = 'Signature ') %>% 
                                          paste(., ': IBK')) %>% 
    pmap(plot_score_forest)
  
  multi_plots$signature_score_formulas <- multi_modeling$cv_pass$models %>% 
    map(get_score_formula)
    
  
# 5-year survival: sensitivity, specificity, J ------
  
  insert_msg('Plotting five-year survival AUC values')
  
  multi_plots$five_surv_cutpoint <- list(x_var = c('Se_IBK_0', 'Sp_IBK_0', 'J_IBK_0'), 
                                         y_var = c('Se_LZ_0', 'Sp_LZ_0', 'J_LZ_0'), 
                                         x_lab = c('Sensitivity, IBK', 
                                                   'Specificity, IBK', 
                                                   'J, IBK'), 
                                        y_lab = c('Sensitivity, LZ/W', 
                                                   'Specificity, LZ/W', 
                                                   'J, LZ/W')) %>% 
    pmap(plot_model_stats, 
         inp_tbl = multi_plots$model_cutstats, 
         highlight = five_surv$score_vars, 
         label_var = 'model_lab', 
         fill_var = 'sign_type', 
         plot_title = 'Risk signature testing: 5-year mortality', 
         plot_subtitle = 'Training IBK and test LZ/W cohort') %>% 
    map(function(x) x + 
          scale_fill_manual(values = globals$signature_color,
                                          labels = globals$signature_labels, 
                                          name = '') + 
          scale_color_manual(values = globals$signature_color, 
                             labels = globals$signature_labels, 
                             name = '')) %>% 
    set_names(c('roc_sensitivity', 
                'roc_specificity', 
                'roc_j'))
  
# Plotting the 5-year mortality AUC values for the signatures and comparators -----
  
  insert_msg('Plotting AUC for 5 year motality')
  
  multi_plots$five_surv_auc <- plot_c(inp_tbl = multi_plots$model_auc, 
                                      c_var = 'AUC', 
                                      plot_title = 'Risk signature testing: 5-year mortality', 
                                      plot_subtitle = 'Training IBK and test LZ/W cohort', 
                                      x_lab = 'AUC') + 
    facet_grid(sign_type ~ dataset, 
               labeller = labeller(.rows = globals$signature_labels, 
                                   .cols = globals$center_labs), 
               space = 'free', 
               scales = 'free_y') + 
    guides(color = F)
  
# Plotting the results of 5-year mortality modeling with logistic regression ----
  
  insert_msg('Plotting logis estimates for the 5 year mortality')
  
  multi_plots$five_surv_logis <- five_surv$logis_modeling$summaries %>% 
    map(filter, 
        is.na(level)) %>% 
    map(mutate, 
        level = '', 
        mod_lab = ifelse(variable %in% pah_study$comparators$variable, 
                         globals$comp_labs[variable], 
                         stri_replace(variable, fixed = 'sign_', replacement = 'Sign ')), 
        var_type = ifelse(variable %in% pah_study$comparators$variable, 
                          'comparator', 
                          'new') %>% 
          factor(c('new', 'comparator'))) %>% 
    map2_dfr(., names(.), ~mutate(.x, cohort = .y)) %>% 
    plot_summ_forest(variable = 'mod_lab', 
                     level = 'level', 
                     estimate = 'estimate', 
                     cohort_var = 'cohort', 
                     p_value = 'p_adj', 
                     facet = F, 
                     plot_title = 'Risk signature testing: 5-year mortality', 
                     plot_subtitle = 'Training IBK and test LZ/W cohort', 
                     x_lab = 'OR, normalized variable', 
                     x_trans = 'log2') + 
    theme(strip.text.y = element_text(size = 8, angle = -90), 
          strip.background.y = element_rect(fill = 'gray95', color = 'gray80')) + 
    facet_grid(var_type ~ cohort, 
               scales = 'free_y', 
               space = 'free', 
               labeller = labeller(.rows = globals$signature_labels, 
                                   .cols = globals$center_labs))

  
# END ----
  
  insert_tail()