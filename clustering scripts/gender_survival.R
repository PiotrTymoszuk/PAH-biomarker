# This script checks for gender-specific differences in survival (overall and 5-year)
# motility (SMWD), daily functioning (WHO FC) and risk scoring in the age classes

  insert_head()
  
# data container -----
  
  gender <- list()
  
# globals: variables, analysis tables and survival objects -----
  
  ## variables, their modeling functions and families and transformation functions
  
  gender$variables <- c(pah_study$comparators$variable, 
                        c('sign_309', 'sign_2525'), 
                        'WHOFc', 
                        'SMWD', 
                        'mPAP', 
                        'cardiac_index', 
                        'NTproBNP')
  
  gender$var_labs <- c(pah_study$comparators$label, 
                       c('Signature 309', 'Signature 2252'), 
                       'WHO FC', 
                       'SMWD', 
                       'mPAP', 
                       'Cardiac index', 
                       'NT-pro-BNP')
  
  gender$axis_labs <- c(rep('Risk score', nrow(pah_study$comparators) + 2), 
                        'WHO FC', 
                        'SMWD, m', 
                        'mPAP, mmHg', 
                        'Cardiac index', 
                        'NT-pro-BNP, pg/ml')

  ## analysis tables, split by the age
  
  gender$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select,
        timepoint, 
        ID, 
        sex, 
        age_class, 
        WHOFc, 
        SMWD, 
        mPAP, 
        cardiac_index, 
        death_study, 
        Survival_time_from_FD_months, 
        any_of(gender$variables)) %>% 
    map2(., 
         map(multi_modeling$score_tbl, 
             select, 
             any_of(gender$variables)), 
         cbind) %>% 
    map(mutate, 
        study_group = interaction(sex, age_class)) %>% 
    map(~dlply(.x, .(age_class), as_tibble)) %>% 
    unlist(recursive = F) %>% 
    set_names(c('IBK_young', 
                'IBK_elderly', 
                'LZ_young', 
                'LZ_elderly'))
  
  ## initial n numbers and plot tags
  
  gender$n_numbers <- gender$analysis_tbl %>% 
    map(count, 
        sex)
  
  gender$plot_tags <- gender$n_numbers %>% 
    map(function(n) paste0('\nMale: n = ', n$n[1], ', Female: n = ', n$n[2]))
  
  ## age-specific survival objects
  
  gender$surv_obj <- gender$analysis_tbl %>% 
    map(create_surv,
        time_variable = 'Survival_time_from_FD_months', 
        event_variable = 'death_study')
  
  ## colors and labels
  
  gender$gender_labs <- c('male' = 'Male', 
                          'female' = 'Female')
  
  gender$gender_colors <- c('male' = 'steelblue', 
                            'female' = 'indianred')
  
  gender$time_labs <- c('IBK, \u2264 65 years', 
                        'IBK, > 65 years', 
                        'LZ/W, \u2264 65 years', 
                        'LZ/W, > 65 years') 
  
  gender$group_labs <- c('M \u226465', 
                         'F \u226465', 
                         'M >65', 
                         'F >65') %>% 
    set_names(levels(gender$analysis_tbl$IBK_young$study_group))
  
  gender$group_colors <- c('steelblue1', 
                           'indianred1', 
                           'steelblue4', 
                           'indianred4') %>% 
    set_names(levels(gender$analysis_tbl$IBK_young$study_group))

# survival differences investigated by KM -----
  
  insert_msg('Overall survival')
  
  ## analysis
  
  gender$os$analysis_obj <- list(inp_table = gender$analysis_tbl, 
                             surv_obj = gender$surv_obj) %>% 
    pmap(model_km_cutpoint, 
         indep_variable = 'sex', 
         strata_labs = gender$gender_labs)
  
  ## summaries
  
  gender$os$summary <- gender$os$analysis_obj %>% 
    map_dfr(km_summary) %>% 
    mutate(subset = names(gender$os$analysis_obj), 
           p_adj = p.adjust(P_value, 'BH'))
  
  ## plots
  
  gender$os$plots <- list(km_model_object = gender$os$analysis_obj, 
                          p_value = signif(gender$os$summary$p_adj, 2), 
                          plot_title = gender$time_labs) %>% 
    pmap(plot_km, 
         palette = unname(gender$gender_colors), 
         legend.labs = unname(gender$gender_labs), 
         legend.title = '') %>% 
    map2(., gender$plot_tags, ~.x + labs(tag = .y))

  
# 5-year mortality, WHO FC and SMWD and comparator risk scoring, plotting -----
  
  insert_msg('Investigating and plotting the differences between the age and gender groups')
  
  ## analyses
  
  gender$funct$analyses <- gender$analysis_tbl %>% 
    reduce(rbind) %>% 
    dlply(.(timepoint)) %>% 
    map(function(x) gender$variables %>% 
          map(safely(analyze_feature), 
              inp_tbl = x, 
              split_var = 'study_group') %>% 
          set_names(gender$variables) %>% 
          map(~.x$result))
  
  ## summaries
  
  gender$funct$summaries <- gender$funct$analyses %>% 
    map(function(x) compact(x) %>% 
          map_dfr(extract_test_summary) %>% 
          filter(test == 'kruskal') %>% 
          mutate(p_adj = p.adjust(p_value, 'BH')))
  
  gender$funct$adj_p_values <- map(gender$funct$summaries, ~.x$p_adj) %>% 
    map(~paste('p =', signif(.x, 2)))
  
  ## plots
  
  gender$funct$plots <- c('IBK_0', 'LZ_0') %>% 
    map(function(x) list(analysis_object = gender$funct$analyses[[x]], 
                         y_lab = gender$axis_labs, 
                         label = paste(gender$var_labs, globals$center_labs[x], sep = ': ')) %>% 
          pmap(safely(plot_analysis), 
               violin = T, 
               labeller = gender$group_labs, 
               fill_colors = unname(gender$group_colors), 
               cust_theme = globals$common_theme, 
               point_alpha = 0.8) %>% 
          map(~.$result) %>% 
          compact) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
  gender$funct$plots <- map2(gender$funct$plots, 
                             gender$funct$adj_p_values, 
                             function(x, y) list(plot = x, 
                                                 sub = y) %>% 
                               pmap(function(plot, sub) plot + labs(subtitle = sub)))
  
  gender$funct$plots$IBK_0$NTproBNP <- gender$funct$plots$IBK_0$NTproBNP + 
    scale_y_continuous(trans = 'log10')
  
  gender$funct$plots$LZ_0$NTproBNP <- gender$funct$plots$LZ_0$NTproBNP + 
    scale_y_continuous(trans = 'log10')
  
# END -----
  
  insert_msg()