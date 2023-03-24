# Random Forest ensemble of the established risk assessment tools
# and of the newly developed Elastic Net score

  insert_head()
  
# container ------
  
  rf_tools <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals ------
  
  insert_msg('Analysis globals')

  ## analysis variables
  
  rf_tools$var_lexicon <- 
    rbind(tibble(variable = 'score', 
                 label = 'ElasticNet'), 
          pah_study$comparators[c('variable', 'label')]) %>% 
    filter(!stri_detect(variable, fixed = 'Reveal'))
  
  ## analysis tables: the established risk scales are treated as factors

  rf_tools$analysis_tbl <- 
    map2(multi_cox$lp_scores, 
         pah_study$data_master %>% 
           blast(timepoint) %>% 
           map(select, 
               ID, any_of(rf_tools$var_lexicon$variable)), 
         left_join, by = 'ID') %>% 
    map(column_to_rownames, 'ID')
  
  for(i in rf_tools$var_lexicon$variable) {
    
    rf_tools$analysis_tbl <- rf_tools$analysis_tbl %>% 
      map(mutate,
          !!i := factor(.data[[i]]))
    
  }
  
  ## model formula
  
  rf_tools$formula <- 
    rf_tools$var_lexicon$variable %>% 
    paste(collapse = '+') %>% 
    paste('Surv(surv_months, death_study) ~', .) %>% 
    as.formula

# tuning of the classifier -------
  
  insert_msg('Tuning')
  
  ## tuning grid
  
  set.seed(1234)

  rf_tools$tuning$grid <- 
    expand.grid(seed = sample(1:1000, 5, replace = FALSE), 
                splitrule = c('logrank', 'logrankscore', 'bs.gradient'), 
                mtry = 2:5, 
                nodesize = 1:20) %>% 
    mutate(splitrule = as.character(splitrule))
  
  ## tuned models

  rf_tools$tuning$models <- 
    rf_tools$tuning$grid[c('splitrule', 'mtry', 'nodesize')] %>% 
    as.list %>% 
    future_pmap(rfsrc, 
                formula = rf_tools$formula, 
                data = rf_tools$analysis_tbl$IBK_0, 
                .options = furrr_options(seed = TRUE))
  
  ## OOB prediction errors (1 - C)
  
  rf_tools$tuning$results <- rf_tools$tuning$models %>% 
    map_dbl(~.x$err.rate[!is.na(.x$err.rate)])
  
  rf_tools$tuning$results <- rf_tools$tuning$grid %>% 
    mutate(c_index = 1 - rf_tools$tuning$results, 
           optimum = ifelse(c_index == max(c_index),
                            'yes', 'no')) %>% 
    arrange(- c_index) %>% 
    as_tibble
  
  rf_tools$tuning$models <- NULL
  
  rf_tools <- compact(rf_tools)
  
  plan('sequential')
  
  gc()
  
# Plotting the tuning results -------
  
  insert_msg('Plotting the tuning results')
  
  ## split rule and mtry
  
  rf_tools$tuning$plots$splitrule <- rf_tools$tuning$results %>% 
    ggplot(aes(x = splitrule, 
               y = c_index)) + 
    facet_grid(. ~ mtry, 
               labeller = as_labeller(set_names(paste('mtry =', 2:5), 
                                                2:5))) + 
    geom_violin(aes(fill = factor(mtry)), 
                alpha = 0.25) + 
    geom_point(aes(color = factor(mtry)), 
               shape = 16, 
               alpha = 0.8, 
               position = position_jitter(width = 0.1, 
                                          height = 0)) + 
    scale_fill_manual(values = c('steelblue4', 
                                 'steelblue2', 
                                 'coral2', 
                                 'coral4')) + 
    scale_color_manual(values = c('steelblue4', 
                                  'steelblue2', 
                                  'coral2', 
                                  'coral4')) + 
    globals$common_theme + 
    labs(title = 'Splitting rule and mtry', 
         subtitle = 'IBK cohort', 
         x = 'Splitting rule', 
         y = 'C-index, out-of-bag')
  
  ## minimal node size and mtry
  
  rf_tools$tuning$plots$nodesize <- rf_tools$tuning$results %>% 
    ggplot(aes(x = nodesize, 
               y = c_index)) + 
    facet_grid(. ~ mtry, 
               labeller = as_labeller(set_names(paste('mtry =', 2:5), 
                                                2:5))) + 
    geom_point(aes(color = factor(mtry)), 
               shape = 16, 
               alpha = 0.8, 
               position = position_jitter(width = 0.1, 
                                          height = 0)) + 
    geom_smooth() +
    scale_fill_manual(values = c('steelblue4', 
                                 'steelblue2', 
                                 'coral2', 
                                 'coral4')) + 
    scale_color_manual(values = c('steelblue4', 
                                  'steelblue2', 
                                  'coral2', 
                                  'coral4')) + 
    globals$common_theme + 
    labs(title = 'Miminal node size and mtry', 
         subtitle = 'IBK cohort', 
         x = 'Minimal node size', 
         y = 'C-index, out-of-bag')
  
# Training the model ------
  
  insert_msg('Training the model')
  
  set.seed(1234)

  rf_tools$rf_models$IBK_0 <- rfsrc(formula = rf_tools$formula, 
                                     data = rf_tools$analysis_tbl$IBK_0, 
                                     mtry = 5, 
                                     nodesize = 2, 
                                     splitrule = 'logrankscore', 
                                     ntree = 1000)
  
# Predictions --------
  
  insert_msg('Predictions')
  
  rf_tools$rf_predictions <- rf_tools$analysis_tbl %>% 
    map(~predict(rf_tools$rf_models$IBK_0, 
                 newdata = .x))
  
# Fit stats -------
  
  insert_msg('Fit stats')
  
  ## C-indexes
  
  rf_tools$fit_stats$c_index <- 
    rf_tools$rf_predictions %>% 
    map_dbl(~get.cindex(time = .x$yvar[, 1], 
                    censoring = .x$yvar[, 2], 
                    predicted = .x$predicted)) %>% 
    compress(names_to = 'cohort', 
             values_to = 'c_index') %>% 
    mutate(c_index = 1 - c_index)
  
  ## brier scores
  
  rf_tools$fit_stats$ibs <- 
    rf_tools$rf_predictions %>% 
    map(get.brier.survival) %>% 
    map_dbl(~.x$crps.std) %>% 
    compress(names_to = 'cohort', 
             values = 'ibs_model')
  
  ## a common table
  
  rf_tools$fit_stats <- rf_tools$fit_stats %>% 
    reduce(left_join, by = 'cohort')

# Variable importance -------
  
  insert_msg('Variable importance')
  
  ## re-fitting the model with the permutation importance metric
  
  rf_tools$importance$model <- 
    rfsrc(formula = rf_tools$formula, 
          data = rf_tools$analysis_tbl$IBK_0, 
          mtry = rf_tools$rf_models$IBK_0$mtry, 
          nodesize = rf_tools$rf_models$IBK_0$nodesize, 
          splitrule = rf_tools$rf_models$IBK_0$splitrule, 
          importance = 'permute')
  
  ## importance plots
  
  rf_tools$importance$plot <- rf_tools$importance$model$importance %>% 
    compress(names_to = 'variable', 
            values_to = 'error') %>% 
    ggplot(aes(x = error, 
               y = reorder(variable, error))) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             fill = 'steelblue') + 
    scale_y_discrete(labels = exchange(rf_tools$var_lexicon$variable, 
                                       dict = rf_tools$var_lexicon)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Risk tool ensemble, variable importance', 
         x = 'Variable importance')

# END -------
  
  rm(i)

  insert_tail()