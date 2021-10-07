# This script performs univariate survival modeling with the selected modeling variables present in both cohorts as independent 
# variables and survival objects (time, event) as responses with Cox regression. Timepoints: 0 and 1 for IBK, 0 for Graz
# Variables: the continuous versions of the the stratified variables are skipped, as agreed with Thomas

# data and toolbox ----

  c('./tools/sys_tools.R', 
    './tools/km_toolbox.R', 
    './tools/modeling_tools.R', 
    './tools/project_tools.R') %>% 
    walk(source)
  
  insert_head()
  
# data container -----
  
  uni_modeling <- list()

# serial modeling and extracting the summaries -----
  
  insert_msg('Serial univariate Cox modeling with the variables present in both cohorts')
  
  ## models
  
  uni_modeling$models <- list(inp_tbl = pah_study[c('IBK_0', 
                                                    'LZ_0')], 
                              surv_object = pah_study$surv_obj[c('IBK_0', 
                                                                 'LZ_0')]) %>% 
    pmap(model_coxph_lst, 
         var_vector = pah_study$mod_variables$variable) %>% 
    map(function(cohort) map(cohort, 
                             ~.$cox_object))
  
  ## summaries, adding the level information
  
  uni_modeling$summary <- uni_modeling$models %>% 
    map(get_cox_results_lst, 
        correct = 'BH') %>% 
    map(mutate, 
        level = pah_study$levels[pah_study$mod_variables$variable] %>% 
          unlist) %>% 
    map2_dfr(., names(.), ~mutate(.x, cohort = .y)) %>% 
    mutate(correlation = ifelse(significant == 'no', 
                                'ns', 
                                ifelse(estimate > 1, 
                                       'positive', 
                                       'negative')), 
           var_lab = translate_vars(variable))
  
  ## model stats
  
  uni_modeling$stats <- uni_modeling$models %>% 
    map(extract_cox_info) %>% 
    map(~.x$stats)

# identifying the variables found to significantly correlate with survival in any of the IBK_0 and LZ_0 data sets ----
  
  insert_msg('Identifying variables significantly correlating with survival in any of the IBK_0 and LZ_0 cohorts')
  
  uni_modeling$top_variables <- uni_modeling$summary %>% 
    dlply(.(cohort)) %>% 
    map(select, 
        variable, 
        significant) %>% 
    map(filter, 
        significant == 'yes') %>% 
    map(~.x$variable) %>% 
    reduce(union)
  
# Generating a Forest plot with the significant survival correlations ----
  
  insert_msg('Displaying the results in a Forest plot')
  
  uni_modeling$forest_plot <- plot_summ_forest(inp_tbl = uni_modeling$summary %>% 
                                                 filter(variable %in% uni_modeling$top_variables), 
                                               plot_title = 'Univariable modeling', 
                                               plot_subtitle = 'Cox proportional hazard model', 
                                               plot_tag = paste('n =', 
                                                                paste(range(uni_modeling$summary$n_complete), 
                                                                      collapse = ' - ')), 
                                               x_lab = 'OR, normalized variables')
# END -----
  
  insert_tail()