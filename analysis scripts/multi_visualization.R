# This script generates plots with the results of multivariate modeling 

  insert_head()
  
# data container -----
  
  multi_plots <- list()
  
# globals: lasso coefficients, fits, stats ------
  
  insert_msg('Globals setup')
  
  ## variable extraction regex
  
  multi_plots$var_regex <- multi_cox$variables %>% 
    sort(decreasing = TRUE) %>% 
    paste(collapse = '|')
  
  ## lasso coefficients
  
  multi_plots$coefs <- multi_cox$glmnet_IBK_0 %>% 
    coef %>% 
    as.matrix %>% 
    as.data.frame %>% 
    filter(s0 != 0) %>% 
    rownames_to_column('parameter') %>% 
    mutate(variable = stri_extract(parameter, 
                                   regex =  multi_plots$var_regex), 
           level = stri_replace(parameter, 
                                regex =  multi_plots$var_regex, 
                                replacement = ''), 
           order = ifelse(level == '_sec', '\u00B2', 
                          ifelse(level == '', '', NA)), 
           level = ifelse(level == '_sec' | level == '', 
                          NA, level), 
           regulation = ifelse(s0 < 0, 
                               'negative', 'positive'), 
           hr = exp(s0), 
           y_txt = ifelse(!is.na(level), 
                          paste(exchange(variable, 
                                          dict = globals$var_labs), 
                                 level, sep = ': '), 
                          paste0(exchange(variable,
                                          dict = globals$var_labs), 
                                 order)), 
           y_txt = stri_replace_all(y_txt, 
                                    regex = 'yes|no|\\s{1}$', 
                                    replacement = ''), 
           y_txt = stri_replace(y_txt, 
                                fixed = 'cm2', 
                                replacement = 'cm\u00B2'), 
           y_txt = stri_replace(y_txt, 
                                fixed = ' ,  ', 
                                replacement = ''))

  ## survfits and difference in survival for the score tertiles

  multi_plots$quant_diffs <- multi_cox$dn_calibration %>% 
    map(~.x$surv_fit) %>% 
    map(surv_pvalue, method = 'survdiff') %>% 
    set_names(c('IBK_0', 'LZ_0')) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(p_adjusted = p.adjust(pval, 'BH'))

  ## n numbers
  
  multi_plots$n_numbers <- multi_cox$lp_scores %>% 
    map(count, death_study)
  
  multi_plots$n_tags <- multi_plots$n_numbers %>% 
    map(~paste('\ntotal: n = ', sum(.x$n), 
               ', events: n = ', .x$n[2]))

# Forest plot with the coefficients of the LASSO model ----
  
  insert_msg('Coefficients of the LASSO model')
  
  multi_plots$coef_plot <- multi_plots$coefs %>% 
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
    facet_grid(variable ~ ., scales = 'free', space = 'free') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank()) + 
    labs(title = 'ElasticNet model estimates', 
         plot_tag = paste0('Total: n = ', 
                           multi_cox$n_numbers$IBK_0['total'], 
                           ', Events: n = ', 
                           multi_cox$n_numbers$LZ_0['events']), 
         x = expression('HR'[ElasticNet]))
  
# Kaplan-Maier plots with the outcome and fitted values -----
  
  insert_msg('KM plots with the model fits and the actual survival')
  
  multi_plots$km_plots <- 
    list(x = multi_cox$coxph_models %>% 
           map(~.x$model), 
         title = c('Training: IBK', 'Test: LZ/W')) %>% 
    pmap(plot, 
         conf.int = TRUE, 
         conf.int.alpha = 0.15, 
         palette = c('steelblue4', 'coral4')) %>% 
    map(~.x$plot) %>% 
    map2(., multi_plots$n_tags, 
         ~.x + 
           scale_color_manual(values = c('steelblue4', 'coral4'), 
                              labels = c('actual', 'predicted')) + 
           scale_fill_manual(values = c('steelblue4', 'coral4'), 
                             labels = c('actual', 'predicted')) + 
           globals$common_theme + 
           labs(tag = .y))

# Kaplan-Meier plots with the score quantiles in the IBK and LZ/W cohorts -------
  
  insert_msg('KM plots for the score quantiles')
  
  ## log-rank test 
  
  multi_plots$km_quart_test <- multi_cox$dn_calibration %>% 
    map(~.x$surv_fit) %>% 
    map(surv_pvalue, method = 'survdiff') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(p_adjusted = p.adjust(pval, 'BH'), 
           significance = ifelse(p_adjusted >= 0.05, 
                                 paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                                 ifelse(p_adjusted < 0.01, 
                                        'p < 0.001', 
                                        paste('p =', signif(p_adjusted, 2)))))
  
  ## Kaplan-Meier plots

  multi_plots$km_quart_plots <- 
    list(x = multi_cox$dn_calibration, 
         title = c('Training: IBK', 'Test: LZ/W')) %>% 
    pmap(plot, 
         palette = c('darkolivegreen', 'steelblue', 'coral3'), 
         xlab = 'Overall survival, months', 
         show_cox = FALSE) %>% 
    map2(., multi_plots$km_quart_test$significance, 
         ~.x + 
           annotate('text', 
                    label = .y, 
                    size = 2.75, 
                    x = 10, 
                    y = 0.1, 
                    hjust = 0, 
                    vjust = 0) + 
           globals$common_theme)

# END -----
  
  insert_tail()
  