# Explorative data analysis and comparison of the study features 
# between the cohorts

# tools -------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(soucer)

  library(exda)
  library(survival)
  library(survminer)

# container list ----
  
  data_ex <- list()
  
# globals -----
  
  insert_msg('Globals setup')
  
  ## a table with variable names and comparison types
  
  data_ex$var_tbl <- 
    tibble(variable = c(pah_study$mod_variables$variable, 
                        'event3', 'event5', 'death_study_fct', 
                        'surv_months', 
                        pah_study$comparators$variable)) %>% 
    filter(!stri_detect(variable, fixed = 'Reveal'))
  
  ## the analysis tables
  
  data_ex$comparator_tbl <- pah_study$data_master %>% 
    filter(timepoint %in% c('IBK_0', 'LZ_0')) %>% 
    select(ID, all_of(pah_study$comparators$variable))
  
  data_ex$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(~left_join(.x, data_ex$comparator_tbl, by = 'ID'))
  
  ## identifying the numeric features

  data_ex$numeric_variables <- data_ex$var_tbl$variable
  
  data_ex$numeric_variables <- 
    data_ex$analysis_tbl$IBK_0[data_ex$numeric_variables] %>% 
    map_lgl(is.numeric) %>% 
    data_ex$var_tbl$variable[.]
  
  ## a table with variable names and comparison types
  
  data_ex$var_tbl <- data_ex$var_tbl %>% 
    mutate(var_type = ifelse(variable %in% data_ex$numeric_variables, 
                                  'numeric', 'factor'), 
           eff_size_type = ifelse(variable %in% data_ex$numeric_variables, 
                                  'wilcoxon_r', 'cramer_v'), 
           plot_lab = paste(exchange(variable, 
                                     dict = pah_study$legend), 
                            exchange(variable, 
                                     dict = pah_study$legend, 
                                     value = 'unit'), 
                            sep = ', '), 
           plot_lab = stri_replace(plot_lab, regex = '\\,\\s{1}NA$', 
                                   replacement = ''))
  
# Normality and EOV (cohort comparison) of the numeric variables ------
  
  insert_msg('Normality and EOV')
  
  data_ex$normality <- data_ex$analysis_tbl %>% 
    map(~explore(.x, 
                 variables = data_ex$numeric_variables, 
                 what = 'normality', 
                 pub_styled = TRUE))
  
  data_ex$eov <- compare_variables(data_ex$analysis_tbl$IBK_0, 
                                   data_ex$analysis_tbl$LZ_0, 
                                   variables = data_ex$numeric_variables, 
                                   what = 'variance', 
                                   pub_styled = TRUE)
  
# descriptive statistics -----
  
  insert_msg('Descriptive stats')
  
  data_ex$desc_stats <- data_ex$analysis_tbl %>% 
    map(~explore(.x, 
                 variables = data_ex$var_tbl$variable, 
                 what = 'table', 
                 pub_styled = TRUE)) %>% 
    reduce(left_join, by = 'variable') %>% 
    set_names(c('variable', 'IBK_0', 'LZ_0'))
  
  data_ex$desc_stats <- rbind(tibble(variable = 'n_number', 
                                     IBK_0 = nrow(data_ex$analysis_tbl$IBK_0), 
                                     LZ_0 = nrow(data_ex$analysis_tbl$LZ_0)), 
                              data_ex$desc_stats)
  
# testing for the differences between the cohorts: non-parametric ----
  
  insert_msg('Testing for the differences between the cohorts')
  
  data_ex$test_results <- 
    compare_variables(data_ex$analysis_tbl$IBK_0, 
                      data_ex$analysis_tbl$LZ_0, 
                      variables = data_ex$var_tbl$variable, 
                      what = 'eff_size', 
                      types = data_ex$var_tbl$eff_size_type, 
                      ci = FALSE, 
                      pub_styled = TRUE, 
                      adj_method = 'BH')

# Violin plots of the numeric variables -----
  
  insert_msg('Violin plots with the numeric features')
  
  data_ex$plots <- 
    list(variable = data_ex$numeric_variables, 
         plot_title = exchange(data_ex$numeric_variables, 
                               dict = pah_study$legend), 
         y_lab = exchange(data_ex$numeric_variables, 
                          dict = data_ex$var_tbl, 
                          value = 'plot_lab'), 
         plot_subtitle = filter(data_ex$test_results, 
                                variable %in% data_ex$numeric_variables) %>% 
           .$significance) %>% 
    pmap(plot_variable, 
         data_ex$analysis_tbl$IBK_0, 
         data_ex$analysis_tbl$LZ_0,
         data_names = globals$center_labs[c('IBK_0', 'LZ_0')], 
         x_lab = 'Cohort',
         cust_theme = globals$common_theme) %>% 
    set_names(data_ex$numeric_variables) %>% 
    map(~.x + 
          scale_fill_manual(values = unname(globals$center_colors[c('IBK_0', 'LZ_0')])))
  
# Comparison of the survival between the cohorts -----
  
  insert_msg('Survival differences between the cohorts')
  
  ## survival stats
  
  data_ex$surv_fits <- 
    surv_fit(Surv(surv_months, death_study) ~ cohort, 
             data = data_ex$analysis_tbl %>% 
               compress(names_to = 'cohort'))
  
  data_ex$surv_summary <- data_ex$surv_fits %>% 
    surv_pvalue(method = 'survdiff') %>% 
    mutate(significance = ifelse(pval < 0.05, 
                                 paste('p =', signif(pval, 2)), 
                                 paste0('ns (p = ', signif(pval, 2), ')')))

  ## n numbers
  
  data_ex$n_numbers <- data_ex$analysis_tbl %>% 
    map(count, death_study)
  
  data_ex$n_tag <- 
    paste0('\nIBK: total: n = ', sum(data_ex$n_numbers$IBK_0$n), 
           ', events: n = ', data_ex$n_numbers$IBK_0$n[2], 
           '\nLZ/W: total: n = ', sum(data_ex$n_numbers$LZ_0$n), 
           ', events: n = ', data_ex$n_numbers$LZ_0$n[2])
  
  ## plotting
  
  data_ex$surv_plot <- 
    ggsurvplot(fit = data_ex$surv_fits, 
               palette = unname(globals$center_colors[c('IBK_0', 'LZ_0')]), 
               title = 'IBK and LZ/W survival differences', 
               xlab = 'Overall survival, months', 
               legend.title = '', 
               legend.labs = unname(globals$center_labs[c('IBK_0', 'LZ_0')]), 
               conf.int = TRUE, 
               conf.int.alpha = 0.15, 
               pval = data_ex$surv_summary$significance, 
               pval.size = 2.75)$plot + 
    globals$common_theme + 
    labs(tag = data_ex$n_tag)
  
# Common result table -----
  
  insert_msg('Common result table')
  
  data_ex$result_table <- 
    left_join(data_ex$desc_stats, 
              data_ex$test_results[c('variable', 'significance', 'eff_size')], 
              by = 'variable')
  
  ## appending with the survival testign results (Mentel-Henszel test)

  data_ex$result_table[data_ex$result_table$variable == 'surv_months', 'significance'] <- 
    data_ex$surv_summary[['significance']][1]
  
  data_ex$result_table[data_ex$result_table$variable == 'surv_months', 'eff_size'] <- NA
  
# END ----
  
  insert_tail()