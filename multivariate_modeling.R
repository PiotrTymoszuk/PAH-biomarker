# this script performs the tasks of multivariate modeling:
# 1) A family of Cox models each with 2 - 4 variables is generated with the IBK_0 cohort
# 2) With all significant estimates, significant LRT and Wald test in the test cohort and significant in 20-fold CV are selected
# 3) The significance of the models is verified in the LZ_0 cohort
# 4) Their performance in predicting overall survival is compared with the established risk scales

# data and toolbox ----

  c('./tools/sys_tools.R', 
    './tools/km_toolbox.R', 
    './tools/modeling_tools.R', 
    './tools/model_space_tools.R',
    './tools/roc_toolbox.R', 
    './tools/lm_qc_tools.R') %>% 
    walk(source)
  
  insert_head()
  
# data containers ----
  
  multi_modeling <- list() ## for the results

# globals, options, list, modeling variables and variable combinations (2 to 4) ----
  
  insert_msg('Setting the globals and making the 2 - 4 variable combination and the comparator list')

  options(future.globals.maxSize= 25000000000) ## to accommodate modeling tasks

# generation of the model space for the training IBK cohort -----
  
  insert_msg('Generation of the trainig data set model space')
  
  multi_modeling$training_models <- gen_model_space(mod_fun = coxph, 
                                                    variable_vect = pah_study$mod_variables$variable, 
                                                    response = 'pah_study$surv_obj$IBK_0', 
                                                    inp_data = pah_study$IBK_0, 
                                                    parallel = T, 
                                                    min_var = 2, 
                                                    max_var = 4, 
                                                    globals = c('pah_study', 
                                                                'make_model', 
                                                                'make_formula', 
                                                                'coxph'))  %>% 
    map(function(x) set_names(x, paste0('sign_', 1:length(x))))
  
# Extracting training model summaries, identification of the models with all-significant betas and LRT  -----
  
  insert_msg('Model summaries and re-distribution stats')
  
  ## model summary, re-distribution stats
  
  multi_modeling$training_summary <- multi_modeling$training_models$models %>% 
    extract_cox_info(.parallel = F)
  
  ## significant models 
  
  multi_modeling$training_signif$stats <- multi_modeling$training_summary$stats %>% 
    filter(signif_estimates == 'yes', 
           p_lrt_adj < 0.05, 
           p_wald_adj < 0.05)
  
  multi_modeling$training_signif[c('vars', 
                                   'models')] <- multi_modeling$training_models %>% 
    map(~.x[multi_modeling$training_signif$stats$model_id])

# Cross-validation of the significant training models, selection of the signifincat models ----
  
  insert_msg('Cross-validation errors')

  multi_modeling$training_signif$cv_results <- multi_modeling$training_signif$models %>% 
    cv_cox(n_folds = 20, 
           seed = 123, 
           new_data = pah_study$IBK_0)
  
  ## conversion to a handy wide table
  
  multi_modeling$training_signif$cv_results <- multi_modeling$training_signif$cv_results %>% 
    dlply(.(stat)) %>% 
    map(select, 
        model_id, 
        mean, 
        lower_ci, 
        upper_ci) %>% 
    map2(., names(.), ~set_names(.x, c('model_id', 
                                       paste(.y, c('mean', 
                                                   'lower', 
                                                   'upper'), 
                                             sep = '_')))) %>% 
    reduce(left_join, 
           by = 'model_id') %>% 
    select(model_id, 
           starts_with('c_index'), 
           mae_mean, 
           mse_mean, 
           rsq_mean) %>% 
    set_names(c('model_id', 
                'c_index', 
                'lower_ci', 
                'upper_ci', 
                'mae', 
                'mse', 
                'rsq')) %>% 
    as_tibble
  
# Filtering out the models failing to pass the cross-validation ----
  
  insert_msg('Filtering out the models failing in the CV')
  
  ## CV passing models
  
  multi_modeling$cv_pass$cv_results <- multi_modeling$training_signif$cv_results %>% 
    filter(lower_ci > 0.5)
  
  ## ... and their stats, summaries and Cox objects
  
  multi_modeling$cv_pass$stats <- multi_modeling$training_signif$stats %>% 
    filter(model_id %in% multi_modeling$cv_pass$cv_results$model_id)
  
  multi_modeling$cv_pass[c('vars', 
                           'models')] <- multi_modeling$training_models %>% 
    map(~.x[multi_modeling$cv_pass$cv_results$model_id])

# External validation of the significant models -----
  
  insert_msg('External validation')
  
  ## score tables
  
  multi_modeling$score_tbl[c('IBK_0', 
                             'LZ_0')] <- pah_study[c('IBK_0', 
                                                     'LZ_0')] %>% 
    map(function(x)  multi_modeling$cv_pass$models %>% 
          map_dfc(predict, 
                  newdata = x, 
                  type = 'lp'))
  
  ## serial modeling
  
  multi_modeling$validation$models <-  list(inp_tbl = multi_modeling$score_tbl, 
                                            surv_object = pah_study$surv_obj) %>% 
    pmap(model_coxph_lst, 
         var_vector = multi_modeling$cv_pass$stats$model_id) %>% 
    map(function(cohort) map(cohort, 
                             ~.$cox_object))
  
  ## summaries
  
  multi_modeling$validation$summary <- multi_modeling$validation$models %>% 
    map(get_cox_results_lst, 
        correct = 'BH')
   
  ## stats
  
  multi_modeling$validation$stats <- multi_modeling$validation$models %>% 
    map(extract_cox_info) %>% 
    map(~.x$stats)

# comparator modeling ----
  
  insert_msg('Cox modeling')
  
  ## serial modeling
  
  multi_modeling$comparators$models <-  list(inp_tbl = pah_study[c('IBK_0', 
                                                                   'LZ_0')], 
                                             surv_object = pah_study$surv_obj) %>% 
    pmap(model_coxph_lst, ## no Reveal scoring for the LZ cohort available 
         var_vector = pah_study$comparators$variable, 
         safely = T) %>% 
    map(function(cohort) map(cohort, 
                             ~.$result)) %>% 
    map(function(cohort) map(cohort, 
                             ~.$cox_object)) %>% 
    map(compact)

  ## summaries
  
  multi_modeling$comparators$summary <- multi_modeling$comparators$models %>% 
    map(get_cox_results_lst, 
        correct = 'BH')
  
  ## stats
  
  multi_modeling$comparators$stats <- multi_modeling$comparators$models %>% 
    map(extract_cox_info) %>% 
    map(~.x$stats)
  
# execution of the downstream analysis and visualization scripts ------
  
  insert_msg('Downstream analysis and visualization')
  
  c('./multivariate scripts/km_analysis_scores.R', 
    './multivariate scripts/five_year_survival.R', 
    './multivariate scripts/multi_visualization.R') %>% 
    walk(source)
  
# END ----
  
  insert_tail()