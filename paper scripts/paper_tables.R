# This script generates paper tables

  insert_head()
  
# data container ----
  
  suppl_tables <- list()
  
# Supplementary Table S1: study variables -----
  
  insert_msg('Table S1: study variables')
  
  suppl_tables$study_vars$legend <- pah_study$legend
  
  suppl_tables$study_vars$levels <- pah_study$mod_variables$variable %>% 
    map(function(x) levels(pah_study$data_master[[x]])) %>% 
    set_names(pah_study$mod_variables$variable) %>% 
    map_chr(paste, collapse = ', ') %>% 
    tibble(variable = names(.), 
           levels = .)
  
  suppl_tables$study_vars <- suppl_tables$study_vars %>% 
    reduce(left_join, by = 'variable') %>% 
    mutate(variable_type = car::recode(variable_type, "'ignore' = 'numeric'")) %>% 
    set_names(c('Variable', 
                'Description', 
                'Variable label', 
                'Unit', 
                'Variable type', 
                'Stratification'))
  
# Supplementary Table S2: results of univariable Cox modeling -----
  
  insert_msg('Table S1: univariable Cox results')
  
  suppl_tables$univariable_cox <- uni_modeling$summary %>% 
    mutate(estimate = paste0(signif(estimate, 3), 
                             '[', 
                             signif(lower_ci, 3), 
                             ' - ', 
                             signif(upper_ci, 3), 
                             ']'), 
           significance = ifelse(p_adjusted < 0.05, 
                                 paste('p =', signif(p_adjusted, 2)), 
                                 'ns'), 
           variable = paste(var_lab, level, sep = ': '), 
           cohort = globals$center_labs[cohort]) %>% 
    select(cohort, 
           variable, 
           estimate, 
           n_complete, 
           significance) %>% 
    set_names(c('Cohort', 
                'Variable', 
                'HR', 
                'N', 
                'pFDR'))
  
# Supplementary Table S3: complete results of multivariate cox modeling for the candidate signatures in the training cohort -----
  
  insert_msg('Table S2: multivariate Cox results')
  
  ## variables
  
  suppl_tables$multivariable_cox$vars <-  
    tibble(model_id = names(multi_modeling$training_models$vars), 
           variables = map(multi_modeling$training_models$vars, 
                           ~translate_vars(.x)) %>% 
             map_chr(paste, 
                     collapse = ', '))
  
  ## model stats
  
  suppl_tables$multivariable_cox$stats <- multi_modeling$training_summary$stats %>% 
    mutate(p_lrt_adj = signif(p_lrt_adj, 2), 
           p_wald_adj = signif(p_wald_adj, 2),
           c_index = paste0(signif(c_index, 2), 
                            ' [', 
                            signif(lower_ci, 2), 
                            ' - ', 
                            signif(upper_ci, 2), 
                            ']'), 
           cohort = 'IBK') %>% 
    select(model_id, 
           signif_estimates,
           p_lrt_adj, 
           p_wald_adj, 
           c_index, 
           n, 
           cohort)
  
  ## complete table
  
  suppl_tables$multivariable_cox <- reduce(suppl_tables$multivariable_cox, 
                                           left_join, 
                                           by = 'model_id') %>% 
    mutate(model_id = stri_replace(model_id, fixed = 'sign_', replacement = 'Signature ')) %>% 
    set_names(c('Signature', 
                'Variables', 
                'Significant HRs',
                'pLRT FDR', 
                'pWald FDR', 
                'C', 
                'N', 
                'Cohort'))
  
# Supplementary Table S4: Formulas of the best signature scores passing the CV validation -----
  
  insert_msg('Table S4: formulas of the best signature scores')
  
  suppl_tables$score_formulas <- tibble(Signature = names(multi_plots$signature_score_formulas), 
                                        `Score formula` = unlist(multi_plots$signature_score_formulas)) %>% 
    mutate(Signature = stri_replace(Signature, fixed = 'sign_', replacement = 'Signature '))
  
# Supplementary Table S5: AUC values of the signature scores and comparators at predicting the 5-year survival -----
  
  insert_msg('Table S5: ROC AUC values for the 5-year survival')
  
  suppl_tables$roc_auc <- five_surv$roc_modeling$stats %>% 
    map2_dfr(., c('IBK', 'LZ/W'), ~mutate(.x, cohort = .y)) %>% 
    mutate(AUC = paste0(signif(AUC, 3), 
                        ' [', 
                        signif(lowerCI, 3), 
                        ' - ', 
                        signif(upperCI, 3), 
                        ']'),
           model_id = ifelse(model_id %in% pah_study$comparators$variable, 
                             globals$comp_labs[model_id], 
                             stri_replace(model_id, fixed = 'sign_', replacement = 'Signature ')), 
           cutoff = signif(cutoff, 2), 
           Se = signif(Se, 2), 
           Sp = signif(Sp, 2), 
           Optimal.criterion = signif(Optimal.criterion, 2)) %>% 
    select(cohort, 
           model_id, 
           cutoff, 
           Se, 
           Sp, 
           Optimal.criterion, 
           AUC) %>% 
    set_names(c('Cohort', 
                'Signature', 
                'Cutoff', 
                'Sensitivity', 
                'Specificity', 
                'J', 
                'AUC'))
  
# Saving on the disc -----
  
  insert_msg('Saving the tables on the disc')
  
  suppl_tables %>% 
    set_names(paste0('Table S', 
                     1:length(suppl_tables))) %>% 
    write_xlsx(path = './paper/supplementary_tables.xlsx')
  
# END -----
  
  insert_tail()
  


    
  
  
  
  