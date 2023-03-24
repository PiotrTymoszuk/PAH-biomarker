# This script generates paper tables

  insert_head()
  
# data container ----
  
  paper_tables <- list()
  suppl_tables <- list()
  
# Table 1: cohort characteristic ------
  
  insert_msg('Table 1: cohort characteristic')
  
  paper_tables$cohort_features <- data_ex$result_table %>% 
    filter(variable %in% c('n_number', 
                           'age_fc', 
                           'sex', 
                           'anemia', 
                           'renal_ins', 
                           'percardial_effusion', 
                           'WHOFc_class', 
                           'SMWD', 
                           'mPAP', 
                           'PVR', 
                           'event5', 
                           'surv_months')) %>% 
    mutate(variable = factor(variable, 
                             c('n_number', 
                               'age_fc', 
                               'sex', 
                               'anemia', 
                               'renal_ins', 
                               'percardial_effusion', 
                               'WHOFc_class', 
                               'SMWD', 
                               'mPAP', 
                               'PVR', 
                               'event5', 
                               'surv_months'))) %>% 
    arrange(variable) %>% 
    format_table
  
  paper_tables$cohort_features <- paper_tables$cohort_features %>% 
    mutate(variable = exchange(as.character(variable), 
                               value = 'plot_lab', 
                               dict = data_ex$var_tbl), 
           variable = ifelse(is.na(variable), 
                             'Participants, n', 
                             variable)) %>% 
    set_names(c('Variable', 
                'IBK', 'LZ/W', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'cohort_features', 
            ref_name = 'cohort_features', 
            caption = paste('Characteristic of the study cohorts.', 
                            'Numeric variables are presented as medians with', 
                            'interquartile ranges (IQR) and ranges.', 
                            'Categorical variables are presented as percents', 
                            'and counts withing the complete observation set.'))
  
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
    set_names(c('Name in R', 'Description', 
                'Name in the report', 'Unit', 
                'Stratification', 'Used in risk modeling'))
  
  suppl_tables$study_vars <- suppl_tables$study_vars %>% 
    filter(`Name in R` %in% unique(c(data_ex$var_tbl$variable, 
                                     pah_study$comparators$variable, 
                                     multi_cox$variables, 
                                     clust$variables, 
                                     cl_chara$var_tbl$IBK_0$variable, 
                                     cl_chara$var_tbl$LZ_0$variable, 
                                     'death_study'))) %>% 
    mdtable(label = 'study_vars', 
            ref_name = 'study_vars', 
            caption = 'Study variables.')

# Supplementary Table S2: supplementary cohort characteristic -----
  
  insert_msg('Table S2: supplementary')

  suppl_tables$cohort_features <- data_ex$result_table %>% 
    filter(variable %in% c('n_number', 'mRAP', 'SO2_RL_class', 
                           'NTproBNP_log', 'cardiac_index', 'RA_area', 
                           'MCV', 'RDW_log', 'FT_log', 'TSAT_log', 
                           'mRASP', 'Compera', 'SPAHR', 
                           'FRENCH3p', 'FRENCH4p', 
                           'event3', 'death_study_fct')) %>% 
    mutate(variable = factor(variable, 
                             c('n_number', 'mRAP', 'SO2_RL_class', 
                               'NTproBNP_log', 'cardiac_index', 'RA_area', 
                               'MCV', 'RDW_log', 'FT_log', 'TSAT_log', 
                               'mRASP', 'Compera', 'SPAHR', 
                               'FRENCH3p', 'FRENCH4p', 
                               'event3', 'death_study_fct'))) %>% 
    arrange(variable) %>% 
    mutate(variable = exchange(as.character(variable), 
                               value = 'plot_lab', 
                               dict = data_ex$var_tbl), 
           variable = ifelse(is.na(variable), 
                             'Participants, n', variable), 
           variable = ifelse(variable %in% c('mRASP', 'COMPERA', 
                                             'SPAHR', 'Reveal Lite', 
                                             'Reveal 2.0'), 
                             paste0(variable, ', risk strata'), 
                             ifelse(stri_detect(variable, fixed = 'FPHR'), 
                                    paste0(variable, ', number of risk factors'), 
                                    variable)), 
           variable = stri_replace(variable, 
                                   fixed = 'cm2', 
                                   replacement = 'cm\u00B2')) %>% 
    format_table %>% 
    set_names(c('Variable', 'IBK', 'LZ/W', 'Significance', 'Effect size'))
  
  suppl_tables$cohort_features <- 
    mdtable(suppl_tables$cohort_features, 
            label = 'cohort_features', 
            ref_name = 'cohort_features', 
            caption = paste('Supplementary characteristic of the study cohorts.', 
                            'Numeric variables are presented as medians with', 
                            'interquartile ranges (IQR) and ranges.', 
                            'Categorical variables are presented as percents', 
                            'and counts withing the complete observation set.'))
  
# Supplementary Table S3: results of univariable Cox modeling -----
  
  insert_msg('Table S3: univariable Cox results')
  
  suppl_tables$univariable_cox <- uni_cox$summary %>% 
    map_dfr(mutate, 
            estimate = paste0(signif(estimate, 2), ' [', 
                              signif(lower_ci, 2), ' - ', 
                              signif(upper_ci, 2), ']'), 
            significance = ifelse(p_adjusted >= 0.05, 
                                  paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                                  ifelse(p_adjusted < 0.001, 'p < 0.001', 
                                         paste('p =', signif(p_adjusted, 2)))), 
            c_index = paste0(signif(c_index, 2), ' [', 
                             signif(c_lower_ci, 2), ' - ', 
                             signif(c_upper_ci, 2), ']'), 
            variable = exchange(variable, 
                                dict = globals$var_labs), 
            rsq_mev = as.character(signif(rsq_mev, 2)), 
            ibs_model = as.character(signif(ibs_model, 2)), 
            cohort = globals$center_labs[cohort]) %>% 
    select(cohort, variable, level, order, 
           estimate, significance, c_index, 
           rsq_mev, ibs_model) %>% 
    set_names(c('Cohort', 'Variable', 'Level', 'Order', 
                'HR, 95% CI', 'Significance', 
                'C index, 95% CI', 'R\u00B2', 
                'IBS'))
  
  suppl_tables$univariable_cox <- 
    mdtable(suppl_tables$univariable_cox, 
            label = 'univariable_cox', 
            ref_name = 'univariable_cox', 
            caption = 'Results of univariable Cox modeling.')

# Supplementary Table S4: performance of the Elastic Net score and comparators -------  
  
  insert_msg('Table S4: elastic net score and established risk tools')
  
  suppl_tables$risk_tools <- 
    surv_tools$fit_stats %>% 
    map_dfr(filter, 
            !stri_detect(variable, fixed = 'RF'), 
            !stri_detect(variable, fixed = 'reveal')) %>% 
    transmute(cohort = globals$center_labs[cohort], 
           variable = exchange(variable, 
                               dict = surv_tools$var_lexicon), 
           c_index = paste0(signif(c_index, 2), ' [', 
                            signif(lower_ci, 2), ' - ', 
                            signif(upper_ci, 2), ']'), 
           rsq_mev = as.character(signif(raw_rsq, 2)), 
           ibs_model = as.character(signif(ibs_model, 2))) %>% 
    set_names(c('Cohort', 'Variable', 
                'C index, 95% CI', 'R\u00B2', 
                'IBS')) %>% 
    mdtable(label = 'risk_tools', 
            ref_name = 'risk_tools', 
            caption = paste('Performance of the Elastic Net', 
                            'score and established risk assessment', 
                            'tools at predicting overall PAH survival.'))
  
# Supplementary Table S4 and S5: characteristic of the participant clusters -----
  
  insert_msg('Tables S4 and S5: characteristic of the participant clusters')
  
  suppl_tables[c('clust_IBK', 'clust_LZ')] <- cl_chara$result_tbl %>% 
    map(filter, 
        variable != 'PCWP') %>% 
    map(mutate, 
        variable = stri_replace(variable, 
                                fixed = 'cm2', 
                                replacement = 'cm\u00B2'), 
        variable = ifelse(variable %in% c('mRASP', 'COMPERA', 
                                          'SPAHR', 'Reveal Lite', 
                                          'Reveal 2.0'), 
                          paste0(variable, ', risk strata'), 
                          ifelse(stri_detect(variable, fixed = 'FPHR'), 
                                 paste0(variable, ', number of risk factors'), 
                                 variable))) %>% 
    map(format_table) %>% 
    map(set_names, 
        c('Variable', 'Cluster #1', 'Cluster #2', 
          'Significance', 'Effect size'))
  
  suppl_tables[c('clust_IBK', 'clust_LZ')] <- 
    list(x = suppl_tables[c('clust_IBK', 'clust_LZ')], 
         label = c('clust_IBK', 'clust_LZ'), 
         ref_name = c('clust_IBK', 'clust_LZ'), 
         caption = c(paste('Characteristic of the participant clusters', 
                           'in the Innsbruck cohort.', 
                           'Numeric variables are presented as medians with', 
                           'interquartile ranges (IQR) and ranges.', 
                           'Categorical variables are presented as percents', 
                           'and counts withing the complete observation set.'), 
                     paste('Characteristic of the participant clusters', 
                           'in the Linz/Vienna cohort.', 
                           'Numeric variables are presented as medians with', 
                           'interquartile ranges (IQR) and ranges.', 
                           'Categorical variables are presented as percents', 
                           'and counts withing the complete observation set.'))) %>% 
    pmap(mdtable)
  
# Saving on the disc -----
  
  insert_msg('Saving the tables on the disc')
  
  ## cover sheet with the table number and caption
  
  suppl_tables$cover <- suppl_tables %>% 
    map_chr(attr, 'caption') %>% 
    tibble(Table = paste0('Table S', 1:length(suppl_tables)), 
           Caption = .)
  
  suppl_tables <- 
    suppl_tables[c('cover', names(suppl_tables)[names(suppl_tables) != 'cover'])]
  
  suppl_tables %>% 
    set_names(c('Cover', 
                paste0('Table S', 
                     1:(length(suppl_tables) - 1)))) %>% 
    write_xlsx(path = './paper/supplementary_tables.xlsx')
  
# END -----
  
  insert_tail()