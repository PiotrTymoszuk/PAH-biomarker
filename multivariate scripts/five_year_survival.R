# This script checks the performance of the newly established risk signatures to predict 5-year survival
# and compares it with the established PAH risk scales

  insert_head()
  
# data container ----
  
  five_surv <- list()
  
# globals: analysis table ----
  
  insert_msg('Globals setup')
  
  ## variables
  
  five_surv$score_vars <- c(multi_modeling$cv_pass$stats$model_id, 
                            pah_study$comparators$variable)
  
  ## analysis table
  
  five_surv$analysis_tbl <- map2(multi_modeling$score_tbl, 
                                 pah_study[c('IBK_0', 'LZ_0')] %>% 
                                   map(select, 
                                       all_of(pah_study$comparators$variable), 
                                       death_acute), 
                                 cbind) %>% 
    map(as_tibble)
  
  ## n numbers of observations and cases for each of the models
  
  five_surv$n_numbers <- five_surv$analysis_tbl %>% 
    map(function(cohort) five_surv$score_vars %>% 
          map_dfr(function(sign_id) tibble(model_id = sign_id, 
                                           total = sum(complete.cases(cohort[c(model_id, 'death_acute')])), 
                                           cases = cohort[c(model_id, 'death_acute')] %>% 
                                             filter(complete.cases(.)) %>% 
                                             count(death_acute) %>% 
                                             .$n %>% 
                                             .[2])))
  
  five_surv$n_tags <- five_surv$n_numbers %>% 
    map(function(cohort) map2(cohort$total, 
                              cohort$cases, 
                              ~paste0('\nTotal: n = ', .x, ', cases: n = ', .y))) %>% 
    map2(., five_surv$n_numbers, 
         ~set_names(.x, .y[['model_id']]))
  
# ROC modeling ----
  
  insert_msg('ROC modeling of the 5-year survival')
  
  ## optimal cutpoints 
  
  five_surv$roc_modeling$cutoffs <- five_surv$analysis_tbl %>% 
    map(function(cohort) five_surv$score_vars %>% 
          map(safely(find_optimal_cutoff), 
              inp_table = cohort, 
              status_variable = 'death_acute') %>% 
          map(~.x$result) %>% 
          set_names(five_surv$score_vars) %>% 
          compact)
  
  ## stats
  
  five_surv$roc_modeling$stats <- five_surv$roc_modeling$cutoffs %>% 
    map(function(cohort) cohort %>% 
          map(get_optimal_cutpoint_coefs) %>% 
          map(~.x[[1]]) %>% 
          map(t) %>% 
          map2_dfr(., names(.), ~mutate(data.frame(.x), model_id = .y)) %>% 
          as_tibble)
  
  ## AUC
  
  five_surv$roc_modeling$stats <- five_surv$roc_modeling$cutoffs %>% 
    map(function(cohort) cohort %>% 
          map(get_auc) %>% 
          reduce(rbind) %>% 
          as_tibble) %>% 
    map2(., five_surv$roc_modeling$cutoffs, 
         ~mutate(.x, model_id = names(.y))) %>% 
    map2(five_surv$roc_modeling$stats, ., 
         left_join, 
         by = 'model_id')
  
# displaying the scores as ROC plots -----
  
  insert_msg('ROC curves')
  
  five_surv$roc_plots <- list(x = five_surv$analysis_tbl, 
                              y = five_surv$roc_modeling$stats, 
                              z = five_surv$n_tags, 
                              cohort = globals$center_labs[c('IBK_0', 'LZ_0')]) %>% 
    pmap(function(x, y, z, cohort) list(signature = multi_modeling$cv_pass$stats$model_id, 
                                        plot_tag = z[multi_modeling$cv_pass$stats$model_id], 
                                        plot_title = stri_replace(multi_modeling$cv_pass$stats$model_id, 
                                                                  fixed = 'sign_', 
                                                                  replacement = 'Signature ') %>% 
                                          paste(., cohort, sep = ': ')) %>% 
           pmap(plot_com_roc, 
                score_tbl = x, 
                stat_tbl = y, 
                show_auc = F) %>% 
           set_names(multi_modeling$cv_pass$stats$model_id))

# logistic modeling, normalized scores -----
  
  insert_msg('Serial logistic modeling of the 5-year mortality')

  five_surv$logis_modeling$models <- five_surv$analysis_tbl %>% 
    map(select, 
        - death_acute) %>% 
    map(function(cohort) cohort %>% map_dfc(~scale(.x)[, 1])) %>% 
    map2(map(five_surv$analysis_tbl, ~.x['death_acute']), ., cbind) %>% 
    map(make_lm_model, 
        response = 'death_acute', 
        indep_variable = five_surv$score_vars, 
        mod_fun = glm, 
        family = 'binomial', 
        est_transf = exp)
  
  five_surv$logis_modeling$summaries <- five_surv$logis_modeling$models %>% 
    map(get_model_summary)

# END -----
  
  insert_tail()