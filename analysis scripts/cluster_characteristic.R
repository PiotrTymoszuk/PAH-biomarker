# Characteristic of the participant clusters. The features are compared with tow-tailed T test

  insert_head()
  
# container list ----
  
  cl_chara <- list()
  
# analysis tables ------
  
  insert_msg('Analysis tables')
  
  ## analysis tables
  
  cl_chara$analysis_tbl <- clust$clust_obj %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('ID', 'clust_id')) %>% 
    map2(., 
         pah_study[c('IBK_0', 'LZ_0')], 
         left_join, by = 'ID')

  cl_chara$analysis_tbl <- 
    map2(cl_chara$analysis_tbl, 
         pah_study$data_master %>% 
           blast(timepoint) %>% 
           map(select, ID, all_of(pah_study$comparators$variable)), 
         left_join, by = 'ID')
  
  cl_chara$analysis_tbl$LZ_0 <- cl_chara$analysis_tbl$LZ_0 %>% 
    select(- starts_with('Reveal'))

# analysis variables --------
  
  insert_msg('Analysis variables')

  ## variable lexicons with variable names, labels, test and plot types
  
  cl_chara$var_tbl <- cl_chara$analysis_tbl %>% 
    map(select, -ID, -clust_id, -death_study, -surv_months) %>% 
    map(~tibble(variable = names(.x), 
                factor = map_lgl(.x, is.factor))) %>% 
    map(mutate, 
        test_type = ifelse(factor, 'cramer_v', 'wilcoxon_r'), 
        plot_type = ifelse(factor, 'stack', 'violin'))
  
  ## variable labels
  
  cl_chara$var_tbl <- cl_chara$var_tbl %>% 
    map2(., c('IBK', 'LZ/W'), 
         ~mutate(.x, 
                 label = exchange(variable, 
                                  dict = pah_study$legend), 
                 plot_title = paste(label, .y, sep = ', '), 
                 y_lab = ifelse(factor, 
                                '% of cluster', 
                                paste(label, 
                                      exchange(variable,
                                               dict = pah_study$legend, 
                                               value = 'unit'), 
                                      sep = ', ')), 
                 tab_lab =  ifelse(factor, 
                                   label, 
                                   paste(label, 
                                         exchange(variable,
                                                  dict = pah_study$legend, 
                                                  value = 'unit'), 
                                         sep = ', '))))
  
  ## numeric variables: identical for both cohorts
  
  cl_chara$numeric_variables <- cl_chara$var_tbl$IBK_0 %>% 
    filter(!factor) %>% 
    .$variable

# N numbers -------
  
  insert_msg('Cluster N numbers')
  
  cl_chara$n_numbers <- clust$clust_obj %>% 
    map(ngroups)
  
  cl_chara$n_tags <- cl_chara$n_numbers  %>% 
    map(~map2_chr(.x$clust_id, .x$n, ~paste0(.x, ': n = ', .y))) %>% 
    map(paste, collapse = ', ') %>% 
    map(~paste0('\n', .x))
  
# Exploration: normality and EOV -----
  
  insert_msg('Normality and EOV')
  
  ## normality is violated seriously for multiple variables
  ## statistical testing with Wilcoxon test!
  
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
  
  cl_chara$desc_stats <- 
    list(data = cl_chara$analysis_tbl, 
         variables = map(cl_chara$var_tbl, ~.x$variable)) %>% 
    pmap(explore, 
         split_factor = 'clust_id', 
         what = 'table', 
         pub_styled = TRUE) %>% 
    map(reduce, left_join, by = 'variable') %>% 
    map(set_names, c('variable', 'clust_#1', 'clust_#2'))

# Testing ------
  
  insert_msg('Testing')
  
  cl_chara$test <- 
    list(x = cl_chara$analysis_tbl, 
         y = map(cl_chara$var_tbl, ~.$variable), 
         z = map(cl_chara$var_tbl, ~.x$test_type)) %>% 
    pmap(function(x, y, z) x %>% 
           compare_variables(split_factor = 'clust_id', 
                             variables = y, 
                             what = 'eff_size', 
                             types = z, 
                             pub_styled = TRUE, 
                             ci = FALSE, 
                             exact = FALSE, 
                             adj_method = 'BH'))
  
# Common significant differences between the clusters ------
  
  insert_msg('Common significant differences between the clusters')

  cl_chara$cmm_significant <- cl_chara$test %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$variable) %>% 
    reduce(intersect)
    
# Common result table -----
  
  insert_msg('Common result table')
  
  cl_chara$result_tbl <- 
    map2(cl_chara$desc_stats, 
         map(cl_chara$test, 
             ~.x[c('variable', 'significance', 'eff_size')]), 
         left_join, by = 'variable') %>% 
    map2(., cl_chara$var_tbl, 
        ~mutate(.x, variable = exchange(variable, 
                                        dict = .y, 
                                        value = 'tab_lab')))
  
  cl_chara$result_tbl <- 
    map2(cl_chara$result_tbl, 
         cl_chara$n_numbers, 
         ~full_rbind(tibble(variable = 'Participants, n', 
                           `clust_#1` = .y$n[1], 
                           `clust_#2` = .y$n[2]), 
                     .x))

# Summary scatter plot -----  
  
  insert_msg('Summary scatter plots')
  
  ## plotting object
  
  cl_chara$test_obj <- 
    list(x = cl_chara$analysis_tbl, 
         y = map(cl_chara$var_tbl, ~.$variable), 
         z = map(cl_chara$var_tbl, ~.x$test_type)) %>% 
    pmap(function(x, y, z) x %>% 
           compare_variables(split_factor = 'clust_id', 
                             variables = y, 
                             what = 'eff_size', 
                             types = z, 
                             pub_styled = FALSE, 
                             ci = FALSE, 
                             exact = FALSE, 
                             adj_method = 'BH')) %>% 
    map2(., cl_chara$var_tbl, 
         ~mutate(.x, variable = exchange(variable, dict = .y)))

  ## plotting
  
  cl_chara$summary_plots <- 
    list(x = cl_chara$test_obj, 
         plot_title = globals$center_labs[c('IBK_0', 'LZ_0')]) %>% 
    pmap(plot, 
         cust_theme = globals$common_theme, 
         show_labels = 'none', 
         point_alpha = 1) %>% 
    map2(., cl_chara$n_tags, 
         ~.x + 
           geom_hline(yintercept = -log10(0.05), 
                      linetype = 'dashed') + 
           labs(tag = .y, 
                subtitle = 'Cluster comparison, Mann-Whitney or \u03C7\u00B2 test', 
                x = 'Effect size, Wilcoxon r or Cramer V', 
                y = expression('-log'[10]*' pFDR')))
  
  ## adding labels for common significant factors
  
  for(i in names(cl_chara$summary_plots)) {
    
    cl_chara$summary_plots[[i]]$data <- cl_chara$summary_plots[[i]]$data %>% 
      mutate(plot_label = ifelse(variable %in% exchange(cl_chara$cmm_significant, 
                                                        dict = cl_chara$var_tbl[[1]]), 
                                 plot_label, NA))
    
  }
  
  cl_chara$summary_plots <- cl_chara$summary_plots %>% 
    map(~.x + 
          geom_text_repel(aes(label = plot_label), 
                          size = 2.75))

# Plots for single variables -------
  
  insert_msg('Plots for single variables')
  
  cl_chara$plots <- c(IBK_0 = 'IBK_0', 
                      LZ_0 = 'LZ_0') %>% 
    map(function(cohort) list(variable = cl_chara$var_tbl[[cohort]]$variable, 
                              type = cl_chara$var_tbl[[cohort]]$plot_type, 
                              plot_title = cl_chara$var_tbl[[cohort]]$plot_title, 
                              y_lab = cl_chara$var_tbl[[cohort]]$y_lab, 
                              plot_subtitle = cl_chara$test[[cohort]]$significance) %>% 
          pmap(plot_variable, 
               cl_chara$analysis_tbl[[cohort]], 
               split_factor = 'clust_id', 
               scale = 'percent', 
               x_lab = 'PAH cluster', 
               txt_size = 2.5, 
               x_n_labs = TRUE, 
               cust_theme = globals$common_theme) %>% 
          set_names(cl_chara$var_tbl[[cohort]]$variable))
  
# END -----
  
  rm(i)

  insert_tail()