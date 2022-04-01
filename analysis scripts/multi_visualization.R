# This script generates plots with the results of multivariate modeling 

  insert_head()
  
# data container -----
  
  multi_plots <- list()
  
# globals: lasso coefficients, fits, stats ------
  
  insert_msg('Globals setup')
  
  ## lasso coefficients
  
  multi_plots$lasso_coefs <- coef(multi_cox$model_IBK_0) %>% 
    as.matrix %>% 
    as.data.frame %>% 
    filter(s0 != 0) %>% 
    rownames_to_column('parameter') %>% 
    mutate(variable = stri_extract(parameter, regex = paste(multi_cox$variables, collapse = '|')), 
           level = stri_replace(parameter, regex = paste(multi_cox$variables, collapse = '|'), replacement = ''), 
           order = ifelse(level == '_sec', '\u00B2', 
                          ifelse(level == '', '', NA)), 
           level = ifelse(level == '_sec' | level == '', NA, level), 
           regulation = ifelse(s0 < 0, 'negative', 'positive'), 
           hr = exp(s0), 
           y_txt = ifelse(!is.na(level), 
                          paste0(translate_vars(variable), ': ', level), 
                          paste0(translate_vars(variable), order)), 
           y_txt = stri_replace_all(y_txt, regex = 'yes|no|\\s{1}$', replacement = ''), 
           y_txt = stri_replace(y_txt, fixed = 'cm2', replacement = 'cm\u00B2'), 
           y_txt = stri_replace(y_txt, fixed = ' ,  ', replacement = ''))

  ## survfits
  
  multi_plots[c('train_fits', 'test_fits')] <- map2(multi_cox[c('score_IBK_0', 'score_LZ_0')], 
                                                    multi_cox[c('coxph_IBK_0', 'coxph_LZ_0')], 
                                                    ~list(outcome = survfit(Surv(surv_months, death_study) ~ 1, data = .x), 
                                                          predicted = survfit(as_coxph(.y$model), data = .x)))

  ## survfits and difference in survival for the score tertiles
  
  multi_plots[c('quant_train_fits', 'quant_test_fits')] <- multi_cox[c('score_IBK_0', 'score_LZ_0')] %>% 
    map(~survfit(Surv(surv_months, death_study) ~ score_quart, data = .x))
  
  multi_plots[c('quant_train_diffs', 'quant_test_diffs')] <- multi_cox[c('score_IBK_0', 'score_LZ_0')] %>% 
    map(~survdiff(Surv(surv_months, death_study) ~ score_quart, data = .x, rho = 0))
  
  multi_plots$diff_summary <- multi_plots[c('quant_train_diffs', 'quant_test_diffs')] %>% 
    map2_dfr(., c('IBK_0', 'LZ_0'), ~tibble(cohort = .y, 
                                            chisq = .x[['chisq']], 
                                            df = length(.x[['obs']]) - 1)) %>% 
    mutate(p_value = 1 - pchisq(chisq, df), 
           p_adjusted = p.adjust(p_value, 'BH'))
  
  ## KM plot captions with the R-squared and C index values
  
  multi_plots$plot_lab <- multi_cox$coxph_stats %>% 
    mutate(plot_cap = paste0('R\u00B2 = ', signif(rsq_mev, 2), 
                             ', C = ', signif(c_index, 2), ' [', 
                             signif(c_lower_ci, 2), ' - ', 
                             signif(c_upper_ci, 2), ']')) %>% 
    .$plot_cap
  
  ## n numbers
  
  multi_plots$n_numbers <- multi_cox[c('score_IBK_0', 'score_LZ_0')] %>% 
    map(count, death_study)
  
  multi_plots$n_tags <- multi_plots$n_numbers %>% 
    map(~paste('\ntotal: n = ', sum(.x$n), 
               'events: n = ', .x$n[2]))

# Forest plot with the coefficients of the LASSO model ----
  
  insert_msg('Coefficients of the LASSO model')
  
  multi_plots$lasso_hr_plot <- multi_plots$lasso_coefs %>% 
    ggplot(aes(x = hr, 
               y = reorder(y_txt, hr), 
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
    facet_grid(variable ~ ., scales = 'free', space = 'free') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    labs(title = 'LASSO model estimates', 
         plot_tag = paste0('Total: n = ', multi_cox$n_numbers$IBK_0['total'], 
                           ', Events: n = ', multi_cox$n_numbers$LZ_0['events']), 
         x = expression('HR'[LASSO]))
  
# Kaplan-Maier plots with the outcome and fitted values -----
  
  insert_msg('KM plots with the model fits and the actual survival')
  
  multi_plots$km_plots <- list(fit = multi_plots[c('train_fits', 'test_fits')], 
                               data = multi_cox$analysis_tbl, 
                               title = c('Training: IBK', 'Test: LZ/W')) %>% 
    pmap(ggsurvplot_combine, 
         conf.int = TRUE, 
         conf.int.alpha = 0.15, 
         palette = c('steelblue4', 'coral4'), 
         legend.labs = c('actual', 'predicted'), 
         legend.title = '') %>% 
    map(~.x$plot) %>% 
    map(~.x + globals$common_theme) %>% 
    map2(., multi_plots$plot_lab, 
         ~.x + labs(subtitle = .y, 
                    x = 'Overall survival, months')) %>% 
    map2(., multi_plots$n_tags, 
         ~.x + labs(tag = .y))
  
# Kaplan-Meier plots with the score quantiles in the IBK and LZ/W cohorts -------
  
  insert_msg('KM plots for the score quantiles')

  multi_plots$km_quart_plots <- list(fit = multi_plots[c('quant_train_fits', 'quant_test_fits')], 
                                     data = multi_cox[c('score_IBK_0', 'score_LZ_0')], 
                                     title = c('Training: IBK', 'Test: LZ/W'), 
                                     pval = signif(multi_plots$diff_summary$p_adjusted, 2)) %>% 
    pmap(ggsurvplot, 
         palette = c('darkolivegreen', 'steelblue', 'coral3'), 
         legend.labs = c('T1', 'T2', 'T3'), 
         legend.title = '', 
         pval.size = 2.75) %>% 
    map(~.x$plot) %>% 
    map(~.x + globals$common_theme) %>% 
    map(~.x + labs(subtitle = 'Linear predictor score quartiles', 
                   x = 'Overall survival, months')) %>% 
    map2(., multi_plots$n_tags, 
         ~.x + labs(tag = .y))
  
# END -----
  
  insert_tail()
  