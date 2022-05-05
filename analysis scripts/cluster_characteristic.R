# Characteristic of the participant clusters. The features are compared with tow-tailed T test

  insert_head()
  
# container list ----
  
  cl_chara <- list()
  
# globals: analysis tables with the cluster assigment schemes ------
  
  insert_msg('Analysis tables and variables')
  
  ## analysis tables
  
  cl_chara$analysis_tbl <- clust[c('clust_obj_train', 
                                   'clust_obj_test')] %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('ID', 'clust_id')) %>% 
    map2(., pah_study[c('IBK_0', 'LZ_0')], left_join, by = 'ID') %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
  cl_chara$comparator_tbl <- pah_study$data_master %>% 
    filter(timepoint %in% c('IBK_0', 'LZ_0')) %>% 
    select(ID, all_of(pah_study$comparators$variable))
  
  cl_chara$analysis_tbl <- cl_chara$analysis_tbl %>% 
    map(left_join, cl_chara$comparator_tbl, by = 'ID')

  ## analysis responses
  
  ## a table with variable names and comparison types
  
  cl_chara$var_tbl <- tibble(variable = c(pah_study$mod_variables$variable, 
                                         'event3', 'event5', 
                                         filter(pah_study$comparators, 
                                                !stri_detect(variable, fixed = 'Reveal'))$variable))
  
  ## identifying the numeric features
  
  cl_chara$numeric_variables <- cl_chara$var_tbl$variable
  
  cl_chara$numeric_variables <- cl_chara$analysis_tbl$IBK_0[cl_chara$numeric_variables] %>% 
    map_lgl(is.numeric) %>% 
    cl_chara$var_tbl$variable[.]
  
  ## a table with variable names and comparison types
  
  cl_chara$var_tbl <- cl_chara$var_tbl %>% 
    mutate(var_type = ifelse(variable %in% cl_chara$numeric_variables, 
                             'numeric', 'factor'), 
           eff_size_type = ifelse(variable %in% cl_chara$numeric_variables, 
                                  'wilcoxon_r', 'cramer_v'), 
           plot_lab = paste(translate_vars(variable), 
                            translate_vars(variable, 'unit'), 
                            sep = ', '), 
           plot_lab = stri_replace(plot_lab, regex = '\\,\\s{1}$', replacement = ''))
  
  ## n numbers
  
  cl_chara$n_numbers <- clust[c('clust_obj_train', 
                                'clust_obj_test')] %>% 
    map(ngroups)
  
  cl_chara$n_tags <- cl_chara$n_numbers  %>% 
    map(~map2_chr(.x$clust_id, .x$n, ~paste0(.x, ': n = ', .y))) %>% 
    map(paste, collapse = ', ') %>% 
    map(~paste0('\n', .x))
  
# Exploration: normality and EOV -----
  
  insert_msg('Normality and EOV')
  
  cl_chara$normality <- cl_chara$analysis_tbl %>% 
    map(~explore(.x, 
                 split_factor = 'clust_id', 
                 variables = cl_chara$numeric_variables, 
                 what = 'normality', 
                 pub_styled = TRUE) %>% 
          map2_dfr(., names(.), ~mutate(.x, clust_id = .y)))
  
  cl_chara$eov <- cl_chara$analysis_tbl %>% 
    map(~compare_variables(.x, 
                           split_factor = 'clust_id', 
                           variables = cl_chara$numeric_variables, 
                           what = 'variance', 
                           pub_styled = TRUE))
  
# Descriptive stats -----
  
  insert_msg('Descriptive stats')
  
  cl_chara$desc_stats <- cl_chara$analysis_tbl %>% 
    map(~explore(.x, 
                 split_factor = 'clust_id', 
                 variables = cl_chara$var_tbl$variable, 
                 what = 'table', 
                 pub_styled = TRUE) %>% 
          reduce(left_join, by = 'variable') %>% 
          set_names(c('variable', 'clust_#1', 'clust_#2')))
  
  cl_chara$desc_stats <- map2(cl_chara$desc_stats, 
                              cl_chara$n_numbers, 
                              ~rbind(tibble(variable = 'n_number', 
                                            `clust_#1` = .y$n[1], 
                                            `clust_#2` = .y$n[2]), 
                                     .x))
    
# Testing ------
  
  insert_msg('Testing')
  
  cl_chara$test_results <- cl_chara$analysis_tbl %>% 
    map(~compare_variables(.x, 
                           split_factor = 'clust_id', 
                           variables = cl_chara$var_tbl$variable, 
                           what = 'eff_size', 
                           types = cl_chara$var_tbl$eff_size_type, 
                           pub_styled = TRUE, 
                           ci = FALSE, 
                           adj_method = 'BH'))
  
# Common result table -----
  
  insert_msg('Common result table')
  
  cl_chara$result_tbl <- map2(cl_chara$desc_stats, 
                              map(cl_chara$test_results, ~.x[c('variable', 'significance', 'eff_size')]), 
                              left_join, by = 'variable')
  
  cl_chara$result_tbl <- cl_chara$result_tbl %>% 
    map(~map_dfc(.x, stri_replace, regex = '\\nComplete:\\s{1}.*', replacement = '') %>% 
          map_dfc(stri_replace, regex = '^no.*\\nyes:\\s{1}', replacement = '') %>% 
          map_dfc(stri_replace, regex = 'Mean.*\\n', replacement = '') %>% 
          map_dfc(stri_replace_all, fixed = '% (', replacement = '% (n = ') %>% 
          map_dfc(stri_replace, fixed = 'Median =', replacement = 'median:') %>% 
          map_dfc(stri_replace, fixed = 'Range', replacement = 'range'))
  
# Single plots of the numeric variables ------
  
  insert_msg('Single violin plots')
  
  cl_chara$plots <- cl_chara$analysis_tbl %>% 
    map(function(cohort) list(variable = cl_chara$numeric_variables, 
                              plot_title = translate_vars(cl_chara$numeric_variables), 
                              y_lab = translate_vars(cl_chara$numeric_variables, 
                                                     value = 'plot_lab', 
                                                     lexicon = cl_chara$var_tbl)) %>% 
          pmap(plot_variable, 
               cohort, 
               split_factor = 'clust_id', 
               type = 'violin', 
               x_lab = 'Cluster', 
               cust_theme = globals$common_theme) %>%
          set_names(cl_chara$numeric_variables)) %>% 
    unlist(recursive = FALSE)
  
  ## adding the colors and p values in the sub-captions
  
  cl_chara$plot_caps <- cl_chara$test_results %>% 
    map(filter, variable %in% cl_chara$numeric_variables) %>% 
    map(~.x$significance) %>% 
    reduce(c)
  
  cl_chara$plots <- map2(cl_chara$plots, 
                         cl_chara$plot_caps, 
                         ~.x + labs(subtitle = .y)) %>% 
    map(~.x + scale_fill_manual(values = globals$cluster_colors, 
                                name = 'Cluster'))
  
# Summary scatter plot -----  
  
  insert_msg('Summary scatter plots')
  
  ## plotting table
  
  cl_chara$test_plot <- cl_chara$analysis_tbl %>% 
    map(~compare_variables(.x, 
                           split_factor = 'clust_id', 
                           variables = cl_chara$var_tbl$variable, 
                           what = 'eff_size', 
                           types = cl_chara$var_tbl$eff_size_type, 
                           pub_styled = FALSE, 
                           ci = FALSE, 
                           adj_method = 'BH') %>% 
          mutate(significant = ifelse(p_adjusted < 0.05, 'significant', 'ns'), 
                 var_lab = translate_vars(variable)))
  
  ## common significant factors
  
  cl_chara$common_signif <- cl_chara$test_plot %>% 
    map(filter, significant == 'significant') %>% 
    map(~.x$variable) %>% 
    reduce(intersect)
  
  cl_chara$test_plot <- cl_chara$test_plot %>% 
    map(mutate, plot_lab = ifelse(variable %in% cl_chara$common_signif, var_lab, NA))
  
  ## plotting
  
  cl_chara$summary_plots <- list(data = cl_chara$test_plot, 
                                 title = globals$center_labs[c('IBK_0', 'LZ_0')], 
                                 tag = cl_chara$n_tags) %>% 
    pmap(function(data, title, tag) data %>% 
           ggplot(aes(x = estimate, 
                      y = -log10(p_adjusted), 
                      fill = significant)) + 
           geom_hline(yintercept = -log10(0.05), 
                      linetype = 'dashed') +  
           geom_point(shape = 21, 
                      size = 2, 
                      alpha = 0.8, 
                      position = position_jitter(width = 0.002, 
                                                 height = 0.01, 
                                                 seed = 1234)) + 
           geom_text_repel(aes(label = plot_lab), 
                           size = 2.75) + 
           scale_fill_manual(values = c(ns = 'gray70', 
                                        significant = 'coral3'), 
                             name = '') + 
           globals$common_theme + 
           labs(title = title, 
                subtitle = 'Cluster comparison, Mann-Whitney or \u03C7\u00B2 test', 
                tag = tag, 
                x = 'Effect size, Wilcoxon r or Cramer V', 
                y = expression('-log'[10]*' pFDR')))

# END -----
  
  insert_tail()