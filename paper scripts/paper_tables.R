# This script generates paper tables

  insert_head()
  
# data container ----
  
  paper_tables <- list()
  suppl_tables <- list()
  
# Table 1: cohort characteristic ------
  
  insert_msg('Table 1: cohort characteristic')
  
  paper_tables$cohort_features <- data_ex$result_table %>% 
    mutate(variable = translate_vars(variable, value = 'plot_lab', lexicon = data_ex$var_tbl), 
           variable = ifelse(variable == '', 'N participants', variable)) %>% 
    set_names(c('Variable', 'IBK', 'LZ/W', 'Significance', 'Effect size'))
  
# Supplementary Table S1: study variables -----
  
  insert_msg('Table S1: study variables')
  
  suppl_tables$study_vars <- pah_study$legend %>% 
    filter(variable_type != 'ignore')
  
  suppl_tables$study_vars  <- pah_study$level_dict %>% 
    group_by(variable) %>% 
    summarise(level = paste(level, collapse = '; ')) %>% 
    left_join(suppl_tables$study_vars , ., by = 'variable') %>% 
    select(variable, description, label, unit, level, variable_type) %>% 
    mutate(variable_type = ifelse(variable_type == 'independent', 'yes', 'no')) %>% 
    set_names(c('Variable', 'Description', 'Label', 'Unit', 'Stratification', 'Used in risk modeling'))

# Supplementary Table S2: results of univariable Cox modeling -----
  
  insert_msg('Table S1: univariable Cox results')
  
  suppl_tables$univariable_cox <- uni_cox$summary %>% 
    map_dfr(mutate, 
            estimate = paste0(signif(estimate, 2), ' [', 
                              signif(lower_ci, 2), ' - ', 
                              signif(upper_ci, 2), ']'), 
            significance = ifelse(p_adjusted >= 0.05, 
                                  paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                                  ifelse(p_adjusted < 0.001, 'p < 0.001', 
                                         paste('p =', signif(p_adjusted, 2)))), 
            c_index = paste0(signif(c_index, 2), ' [', signif(c_lower_ci, 2), ' - ', 
                             signif(c_upper_ci, 2), ']'), 
            variable = translate_vars(variable), 
            rsq_mev = as.character(signif(rsq_mev, 2)), 
            cohort = globals$center_labs[cohort]) %>% 
    select(cohort, variable, level, order, estimate, significance, c_index, rsq_mev) %>% 
    set_names(c('Cohort', 'Variable', 'Level', 'Model order', 'HR', 'Significance', 'C index', 'R\u00B2'))

# Table S2 and S3: characteristic of the participant clusters -----
  
  insert_msg('Tables S2 and S3: characteristic of the participant clusters')
  
  suppl_tables[c('clust_IBK', 'clust_LZ')] <- cl_chara$result_tbl %>% 
    map(mutate, 
        variable = translate_vars(variable, value = 'plot_lab', lexicon = data_ex$var_tbl), 
           variable = ifelse(variable == '', 'N participants', variable), 
        variable = stri_replace(variable, fixed = 'cm2', replacement = 'cm\u00B2')) %>% 
    map(set_names, 
        c('Variable', 'Cluster #1', 'Cluster #2', 'Significance', 'Effect size'))
  
# Saving on the disc -----
  
  insert_msg('Saving the tables on the disc')
  
  suppl_tables %>% 
    set_names(paste0('Table S', 
                     1:length(suppl_tables))) %>% 
    write_xlsx(path = './paper/supplementary_tables.xlsx')
  
# END -----
  
  insert_tail()