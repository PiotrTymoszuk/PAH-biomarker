# Univariate cox modeling. Numeric variables are splinned with rms::rcs

  insert_head()
  
# container list -----
  
  uni_cox <- list()
  
# globals: median-centered analysis table -----
  
  insert_msg('Analysis table')
  
  uni_cox$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, - surv_months, - death_study) %>% 
    map(~map_dfc(.x, function(var) if(is.numeric(var)) scale(var, center = median(var))[, 1] else var)) %>% 
    map2(., 
         map(pah_study[c('IBK_0', 'LZ_0')], ~.x[c('ID', 'surv_months', 'death_study')]), 
         left_join, by = 'ID')
  
  ## n numbers
  
  uni_cox$n_numbers <- uni_cox$analysis_tbl %>% 
    map(count, death_study) %>% 
    map(~paste0('total: n = ', sum(.x$n), 
                ', events: n = ', sum(.x$n[2])))  %>% 
    map2_chr(c('IBK', 'LZ/W'), ., paste, sep = ': ') %>% 
    paste(collapse = '\n') %>% 
    paste0('\n', .)
  
# serial modeling ------
  
  insert_msg('Serial univariable Cox modeling')
  
  uni_cox$test_results <- uni_cox$analysis_tbl %>% 
    map(function(cohort) pah_study$mod_variables$variable %>% 
          map(model_spline_cox, 
              data = cohort, 
              event_variable = 'death_study', 
              time_variable = 'surv_months') %>% 
          set_names(pah_study$mod_variables$variable))

# model summaries, identification of the significantly regulated factors in each cohort -----
  
  insert_msg('Model summaries')
  
  uni_cox$summary <- uni_cox$test_results %>% 
    map(~map_dfr(.x, ~.x$summary)) %>% 
    map(mutate, 
        p_adjusted = p.adjust(p_value, 'BH')) %>% 
    map(group_by, variable) %>% 
    map(mutate, significant = ifelse(any(p_adjusted < 0.05), 'significant', 'ns')) %>% 
    map(ungroup) %>% 
    map2(., names(.), ~mutate(.x, cohort = .y))
    
  uni_cox$signif_fct <- uni_cox$summary %>% 
    map(filter, significant == 'significant') %>% 
    map(~.x$variable) %>% 
    reduce(union)

# forest plot with the factors significant in at least one of the cohorts -----
  
  insert_msg('Forest plot')
  
  uni_cox$forest_plot <- plot_summ_forest(inp_tbl = uni_cox$summary %>% 
                                            reduce(rbind) %>% 
                                            filter(variable %in% uni_cox$signif_fct), 
                                          plot_title = 'Univariable modeling', 
                                          plot_subtitle = 'Cox proportional hazard model', 
                                          x_lab = 'HR, median-centered', 
                                          x_trans = 'log2', 
                                          plot_tag = uni_cox$n_numbers)
  
# END ----
  
  insert_msg()