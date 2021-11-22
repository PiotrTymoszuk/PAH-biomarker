# A medley of functional project tools

  require(plyr)
  require(tidyverse)
  require(rlang)
  require(cowplot)
  require(survival)
  require(stringi)
  require(caret)
  
# import helpers ----
  
  det_mortality <- function(inp_tbl, update = T) {
    
    # identifies if the patient survived in the study period: 0 - survivor, 1 - deceased
    
    mort_tbl <- inp_tbl %>% 
      ddply(.(ID), summarize, death_study = if(sum(dead_alive, na.rm = T) > 0) 1 else 0) %>% 
      as_tibble
    
    if(!update) {
      
      return(mort_tbl)
      
    } else {
      
      return(left_join(inp_tbl, mort_tbl, by = 'ID'))
      
    }
    
  }
  
  recode_01 <- function(vector) {
    
    car::recode(vector, "0 = 'no'; 1 = 'yes'") %>% 
      factor(c('no', 'yes'))
    
    
  }
  
  binarize <- function(vector, cutoff, right = T, ...) {
    
    if(right) {
      
      vec_labels <- paste0(c('\u2264', '>'), signif(cutoff, 3))
      
    } else {
      
      vec_labels <- paste0(c('<', '\u2265'), signif(cutoff, 3))
      
    }
    
    cutpoints <- c(-Inf, cutoff, Inf)
    
    cut(vector, 
        breaks = cutpoints, 
        labels = vec_labels, ...)
    
  }
  
  cut_median <- function(vector, right = T, ...) {
    
    binarize(vector = vector, 
             cutoff = median(vector, na.rm = T), 
             right = right, ...)
    
  }
  
  cut_quartile <- function(vector, ...) {
    
    quart <- list(prime = c(0.25, 0.5, 0.75), 
                  fall1 = c(0.25, 0.5),
                  fall2 = c(0.5, 0.75), 
                  fall3 = 0.5) %>% 
      map(quantile, 
          x = vector, 
          na.rm = T) %>% 
      map(function(x) c(-Inf, x, Inf))
    
    labs <- list(prime = c('Q1', 'Q2', 'Q3', 'Q4'), 
                 fall1 = c('Q1', 'Q2', 'Q3/Q4'), 
                 fall2 = c('Q1/Q2', 'Q3', 'Q4'), 
                 fall3 = c('Q1/Q2', 'Q3/Q4'))
    
    for(i in names(quart)) {
      
      cut_vector <- try(cut(vector, 
                            quart[[i]], 
                            labs[[i]], ...), silent = T)
      
      if(class(cut_vector) != 'try-error') {
        
        return(cut_vector)
        
      }

    }

    return(i)
    
  }
  
  min_max <- function(vector) {
    
    stopifnot(is.numeric(vector))
    
    rec_vector <- (vector - min(vector, na.rm = T))/(max(vector, na.rm = T) - min(vector, na.rm = T))
    
    return(rec_vector)
    
  }
  
  set_rownames <- function(inp_tbl, vector) {
    
    inp_tbl <- data.frame(inp_tbl)
    
    rownames(inp_tbl) <- vector

    return(inp_tbl)    
    
  }
  
# variable translation -----
  
  translate_vars <- function(variable_vector, show_units = F) {
    
    ## translates the variable vector into a vector of corresponding labels
    
    if(class(variable_vector) == 'list') {
      
      return(variable_vector %>% 
               map(globals$translate_vars, 
                   show_units = show_units))
      
    }
    
    lab_vec <- globals$var_labs[variable_vector] %>% 
      unname
    
    if(show_units) {
      
      return(lab_vec)
      
    } else {
      
      return(lab_vec %>% 
               stri_replace(regex = ',.*', replacement = ''))
      
    }
    
  }
  
# wrappers for Cox modeling and model validation -----
  
  model_coxph_lst <- function(inp_tbl, surv_object, var_vector, safely = F) {
    
    ## a wrapper for the model_coxph() function from modeling_tool.R to handle variable vectors
    
    if(!safely) {
      
      cox_res <- var_vector %>% 
        map(model_coxph, 
            surv_object = surv_object, 
            inp_table = inp_tbl) %>% 
        set_names(var_vector)
      
    } else {
      
      cox_res <- var_vector %>% 
        map(safely(model_coxph), 
            surv_object = surv_object, 
            inp_table = inp_tbl) %>% 
        set_names(var_vector)
      
    }
    
    return(cox_res)
    
  }
  
  get_cox_results_lst <- function(cox_model_list, 
                                  exponentiate = T, 
                                  r_stats = T, 
                                  tibble_output = T, 
                                  correct = 'BH') {
    
    ## a wrapper for the get_cox_results() function handling lists
    
    cox_summary <- cox_model_list %>% 
      map(get_cox_results, 
          exponentiate = exponentiate, 
          r_stats = r_stats)
    
    if(!tibble_output) {
      
      return(cox_summary)
      
    } else {
      
      cox_summary <- suppressWarnings(map2_dfr(cox_summary, 
                                               names(cox_summary),  
                                               function(x, y) mutate(x, variable = y)) %>% 
                                        as_tibble %>% 
                                        mutate(level = stri_replace(parameter, fixed = variable, replacement = '')))
      
      
      
      
      if(correct == 'none') {
        
        return(cox_summary)
        
      } else {
        
        cox_summary <- cox_summary %>% 
          mutate(p_adjusted = p.adjust(p_value, method = correct), 
                 significant = ifelse(p_adjusted < 0.05, 'yes', 'no'), 
                 plot_label = paste(signif(estimate, 2), 
                                    ' (', 
                                    signif(lower_ci, 2), 
                                    ', ', 
                                    signif(upper_ci, 2), 
                                    '), p = ', 
                                    signif(p_adjusted, 2), 
                                    sep = ''))
        
        return(cox_summary)
        
      }
      
    }
    
  }
  
  get_cox_stats <- function(cox_model) {
    
    mod_summary <- summary(cox_model)
    
    mod_stats <- tibble(n = cox_model$n, 
                        n_event = cox_model$nevent, 
                        p_lrt = mod_summary$logtest[['pvalue']], 
                        p_wald = mod_summary$waldtest[['pvalue']], 
                        c_index = cox_model$concordance[['concordance']], 
                        lower_ci = cox_model$concordance[['concordance']] + cox_model$concordance[['std']] * qnorm(0.025), 
                        upper_ci = cox_model$concordance[['concordance']] + cox_model$concordance[['std']] * qnorm(0.975), 
                        rsq = mod_summary$rsq[[1]], 
                        aic = AIC(cox_model), 
                        mse = mean(cox_model$residuals^2), 
                        mae = mean(abs(cox_model$residuals)))
    
    return(mod_stats)
    
  }
  
  extract_cox_info <- function(cox_model, .parallel = F) {
    
    ## extracts estimates and measures of model performance in the training data set
    
    if(class(cox_model) == 'list') {
      
      start_time <- Sys.time()
      message(paste('Extracting stats for', length(cox_model), 'models'))
      on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
      
      if(.parallel) {
        
        plan('multisession')
        
        stat_results <- cox_model %>% 
          future_map(extract_cox_info) %>% 
          set_names(names(cox_model)) %>% 
          transpose
        
        plan('sequential')
        
      } else {
        
        stat_results <- cox_model %>% 
          map(extract_cox_info) %>% 
          set_names(names(cox_model)) %>% 
          transpose
        
      }
      
      stat_results$estimates <- stat_results$estimates %>% 
        map2(., names(.), ~mutate(.x, model_id = .y))
      
      stat_results$stats <- stat_results$stats %>% 
        map2(., names(.), ~mutate(.x, model_id = .y)) %>% 
        do.call('rbind', .) %>% 
        mutate(p_lrt_adj = p.adjust(p_lrt, 'BH'), 
               p_wald_adj = p.adjust(p_wald, 'BH'))
      
      return(stat_results)
      
    }
    
    mod_est <- get_cox_results(cox_model = cox_model, 
                               exponentiate = T)
    
    mod_stats <- get_cox_stats(cox_model = cox_model) %>% 
      mutate(signif_estimates = ifelse(all(mod_est$p_value < 0.05), 'yes', 'no'))
     
    return(list(estimates = mod_est, 
                stats = mod_stats))
    
  }
  
  get_c_fold <- function(cox_model, new_data = NULL, fold_obs = NULL) {
    
    ## gets prediction stats for the data fold
    
    if(is.null(fold_obs)) {
      
      return(get_cox_stats(cox_model))
      
    }
    
    ## modifying the survival objects and data sets
    
    train_surv <- cox_model$y[-fold_obs, ]
    test_surv <- cox_model$y[fold_obs, ]
    
    if(is.null(new_data)) {
      
      new_data <- model.frame(cox_model)
      
    }
    
    train_data <- new_data[-fold_obs, ]
    test_data <- new_data[fold_obs, ]
    
    ## re-fitting the input model to the training data
    
    train_formula <- paste('train_surv ~', paste(names(cox_model$xlevels), collapse = '+')) %>% 
      as.formula
    
    train_model <- coxph(formula = train_formula, 
                         data = train_data)
    
    ## calculating the lp score for the test data set, fitting a cox model
    
    test_formula <- 'test_surv ~ lp_score' %>% 
      as.formula
    
    test_data <- test_data %>% 
      mutate(lp_score = predict(train_model, 
                                newdata = test_data, 
                                type = 'lp'))
    
    test_model <- coxph(formula = test_formula, 
                        data = test_data)
    
    return(get_cox_stats(test_model))
    
  }
  
  cv_cox <- function(cox_model, 
                     new_data = NULL, 
                     n_folds = 20, 
                     seed = 1234, 
                     detailed = T, 
                     hide_errors = T) {
    
    if(class(cox_model) == 'list') {
      
      res_tbl <- cox_model %>% 
        map(cv_cox, 
            n_folds = n_folds, 
            new_data = new_data, 
            seed = seed, 
            detailed = F, 
            hide_errors = T) %>% 
       map2_dfr(., names(.), ~mutate(.x, model_id = .y))
      
      return(res_tbl)
      
    }
    
    ## performs n-fold CV of a Cox model
    
    set.seed(seed = seed)
    
    folds <- createFolds(cox_model$y, k = n_folds) %>% 
      map(as.numeric)
    
    fold_stats <- folds %>% 
      map(safely(get_c_fold), 
          cox_model = cox_model, 
          new_data = new_data) %>% 
      transpose
    
    fold_stats$result <- fold_stats$result %>% 
      compact %>% 
      map2_dfr(., names(.), ~mutate(.x, fold_id = .y))
    
    fold_stats$error <- compact(fold_stats$error)
    
    if(length(fold_stats$error) > 0) {
      
      warning(paste(length(fold_stats$error), 'modeling faults occured'))
      
    }
    
    if(detailed) {
      
      if(hide_errors) {
        
        return(fold_stats$result)
        
      } else {
        
        return(fold_stats)
        
      }
     
      
    }
    
    cv_stats <- c('c_index', 
                  'rsq', 
                  'mse', 
                  'mae') %>% 
      map_dfr(~c(mean(fold_stats$result[[.x]], na.rm = T), 
                 quantile(fold_stats$result[[.x]], c(0.025, 0.975), na.rm = T)) %>% 
                set_names(c('mean', 'lower_ci', 'upper_ci'))) %>% 
      mutate(stat =  c('c_index', 
                       'rsq', 
                       'mse', 
                       'mae'))
    
    return(cv_stats)
    
  }
  
# displaying the modeling results -----
  
  plot_summ_forest <- function(inp_tbl, 
                               variable = 'var_lab', 
                               level = 'level', 
                               estimate = 'estimate', 
                               lower_ci = 'lower_ci', 
                               upper_ci = 'upper_ci', 
                               cohort_var = 'cohort', 
                               p_value = 'p_adjusted', 
                               x_trans = 'log2', 
                               show_label = T, 
                               plot_title = NULL, 
                               plot_subtitle = NULL, 
                               plot_tag = NULL, 
                               x_lab = 'HR', 
                               color_lab = 'Mortality\ncorrelation', 
                               cutpoint = 1, 
                               signif_digits = 2, 
                               facet = T) {
    
    ## plots a two-cohort Forest plot
    
    ## meta
    
    grid_formula <- expr(!!ensym(variable) ~ !!ensym(cohort_var)) %>% 
      as.formula
    
    ## plot
    
    inp_tbl <- inp_tbl %>% 
      mutate(axis_lab = ifelse(.data[[level]] %in% c('no', 'yes', ''), 
                               .data[[variable]], 
                               paste(.data[[variable]], .data[[level]], sep = ': ')), 
             est_lab = paste0(signif(.data[[estimate]], signif_digits),
                              ' [', 
                              signif(.data[[lower_ci]], signif_digits), 
                              ' - ', 
                              signif(.data[[upper_ci]], signif_digits), 
                              ']'), 
             significant = ifelse(.data[[p_value]] < 0.05, 'yes', 'no'), 
             corr_sign = ifelse(significant == 'no', 
                                'ns', 
                                ifelse(.data[[estimate]] > cutpoint, 
                                       'positive', 'negative')))
    
    forest_plot <- inp_tbl %>% 
      ggplot(aes(x = .data[[estimate]], 
                 y = reorder(axis_lab, .data[[estimate]]), 
                 color = corr_sign, 
                 shape = .data[[cohort_var]])) + 
      geom_vline(xintercept = cutpoint, 
                 linetype = 'dashed') + 
      geom_errorbarh(aes(xmin = .data[[lower_ci]], 
                         xmax = .data[[upper_ci]]), 
                     height = 0) + 
      geom_point(size = 2) + 
      scale_shape_manual(values = c(IBK_0 = 15, 
                                    LZ_0 = 16), 
                         labels = globals$center_labs, 
                         name = 'Cohort') + 
      scale_x_continuous(trans = x_trans) + 
      scale_color_manual(values = globals$pos_neg_scale, 
                         name = color_lab) + 
      guides(shape = F) +
      globals$common_theme + 
      theme(axis.title.y = element_blank(), 
            panel.grid.major = element_line(color = 'gray90'), 
            strip.background.y = element_blank(), 
            strip.text.y = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = x_lab)
    
    if(show_label) {
      
      forest_plot <- forest_plot + 
        geom_text(aes(label = est_lab), 
                  size = 2.5, 
                  hjust = 0.5, 
                  vjust = -1)
      
    }
    
    if(facet) {
      
      forest_plot <- forest_plot + 
        facet_grid(grid_formula, 
                   labeller = labeller(.cols = globals$center_labs), 
                   space = 'free', 
                   scales = 'free_y')
      
    }
    
    return(forest_plot)
    
  }
  
  plot_model_stats <- function(inp_tbl, 
                               x_var, 
                               y_var, 
                               highlight = NULL, 
                               label_var = NULL, 
                               lower_x = NULL, 
                               upper_x = NULL, 
                               lower_y = NULL, 
                               upper_y = NULL, 
                               fill_var = NULL, 
                               default_fill = 'steelblue', 
                               plot_title = NULL, 
                               plot_subtitle = NULL, 
                               plot_tag = NULL, 
                               x_lab = x_var, 
                               y_lab = y_var) {
    
    ## makes a custom point plot with the chosen model statistics
    
    if(!is.null(highlight) & !is.null(label_var)) {
      
      inp_tbl <- inp_tbl %>% 
        mutate(plot_lab = ifelse(model_id %in% highlight, 
                                 .data[[label_var]], 
                                 NA))
      
    }
    
    ## plot
    
    if(is.null(fill_var)) {
      
      point_plot <- inp_tbl %>% 
        ggplot(aes(x = .data[[x_var]], 
                   y = .data[[y_var]]))
      
    } else {
      
      point_plot <- inp_tbl %>% 
        ggplot(aes(x = .data[[x_var]], 
                   y = .data[[y_var]], 
                   fill = .data[[fill_var]], 
                   color = .data[[fill_var]]))
      
    }
    
    if(!is.null(lower_x) & !is.null(lower_y) & !is.null(upper_x) & !is.null(upper_y)) {
      
      point_plot <- point_plot + 
        geom_rect(aes(xmin = .data[[lower_x]], 
                      xmax = .data[[upper_x]], 
                      ymin = .data[[lower_y]], 
                      ymax = .data[[upper_y]]), 
                  alpha = 0.15)
      
    } 
    
    if(is.null(fill_var)) {
      
      point_plot <- point_plot + 
        geom_point(size = 2, 
                   shape = 21, 
                   fill = default_fill, 
                   color = 'black')
      
    } else {
      
      point_plot <- point_plot + 
        geom_point(size = 2, 
                   shape = 21, 
                   color = 'black')
      
    }
    
    if(!is.null(highlight) & !is.null(label_var)) {
      
      point_plot <- point_plot + 
        geom_text_repel(aes(label = .data[[label_var]]), 
                        size = 2.5, 
                        box.padding = 0.1)
      
    }
    
    point_plot <- point_plot + 
      globals$common_theme + 
      theme(panel.grid.major = element_line(color = 'gray90')) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = x_lab, 
           y = y_lab)
    
    return(point_plot)
    
  }
  
  plot_c <- function(inp_tbl, 
                     model_var = 'model_lab', 
                     c_var = 'c_index', 
                     lower_ci = 'lower_ci', 
                     upper_ci = 'upper_ci', 
                     cohort_var = 'dataset', 
                     show_label = T, 
                     plot_title = NULL, 
                     plot_subtitle = NULL, 
                     plot_tag = NULL, 
                     x_lab = 'C-index', 
                     color_lab = 'Cohort', 
                     signif_digits = 2) {
    
    ## presents C index values in a forest plot
    
    inp_tbl <- inp_tbl %>% 
      mutate(c_lab = paste0(signif(.data[[c_var]], signif_digits), 
                            ' [', 
                            signif(.data[[lower_ci]], signif_digits), 
                            ' - ', 
                            signif(.data[[upper_ci]], signif_digits), 
                            ']'))
    
    
    forest_plot <- inp_tbl %>% 
      ggplot(aes(x = .data[[c_var]], 
                 y = reorder(.data[[model_var]], .data[[c_var]]), 
                 color = .data[[cohort_var]], 
                 shape = .data[[cohort_var]])) + 
      geom_vline(xintercept = 0.5, 
                 linetype = 'dashed') + 
      geom_errorbarh(aes(xmin = lower_ci, 
                         xmax = upper_ci), 
                     height = 0) + 
      geom_point(size = 2) + 
      scale_color_manual(values = globals$center_colors, 
                         labels = globals$center_labs, 
                         name = color_lab) + 
      scale_shape_manual(values = c(IBK_0 = 16, 
                                    cv = 17, 
                                    LZ_0 = 15)) + 
      guides(shape = F) + 
      globals$common_theme + 
      theme(axis.title.y = element_blank(), 
            panel.grid.major = element_line(color = 'gray90')) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle,
           tag = plot_tag, 
           x = x_lab)
    
    if(show_label) {
      
      forest_plot <- forest_plot + 
        geom_text(aes(label = c_lab), 
                  size = 2.5, 
                  hjust = 0.5, 
                  vjust = -0.8)
      
    }
    
    return(forest_plot)
    
  }
  
  
  get_score_formula <- function(cox_model, 
                                signif_digits = 2, 
                                level_dictionary = pah_study$level_dict) {
    
    ## extracts the score formula from a cox object
    
    dic_vector <- set_names(level_dictionary$param_label, 
                            level_dictionary$parameter)
    
    cox_model_terms <- dic_vector[names(cox_model$coefficients)]
    cox_model_estimates <- signif(cox_model$coefficients, signif_digits)
    
    score_formula <- map2(cox_model_estimates, 
                          cox_model_terms, 
                          ~paste(.x, .y, sep = ' \u00D7 ')) %>%  
      paste(collapse = ' + ')
    
    return(score_formula)
    
  }
  
  plot_score_forest <- function(cox_model, 
                                signif_digits = 2, 
                                level_dictionary = pah_study$level_dict, 
                                x_trans = 'log2', 
                                plot_title = NULL) {
    
    ## plots cox model estimates in a forest plot
    
    summ_tbl <- get_cox_results(cox_model = cox_model, 
                                exponentiate = T) %>% 
      mutate(est_lab = paste0(signif(estimate, signif_digits), 
                              ' [', 
                              signif(lower_ci, signif_digits), 
                              ' - ', 
                              signif(upper_ci, signif_digits), 
                              ']'), 
             significant = ifelse(p_value < 0.05, 'yes', 'no'),
             corr_sign = ifelse(significant == 'yes', 
                                ifelse(estimate > 1, 'positive', 'negative'), 
                                'ns'))
    
    dic_vector <- set_names(level_dictionary$param_label, 
                            level_dictionary$parameter)
    
    y_scale <- dic_vector[summ_tbl[['parameter']]]
    
    plot_tag <- paste0('\nTotal: n = ', 
                       cox_model$n, 
                       ', cases: n = ', 
                       cox_model$nevent)
    
    model_stats <- get_cox_stats(cox_model)
    
    plot_subtitle <- paste0('MSE = ', 
                           signif(model_stats$mse[1], 2), 
                           ', C = ', 
                           signif(model_stats$c_index[1], 2), 
                           ' [', 
                           signif(model_stats$lower_ci, 2), 
                           ' - ', 
                           signif(model_stats$upper_ci, 2), 
                           '], pLRT = ', 
                           signif(model_stats$p_lrt[1], 2))
    
    ## plotting
    
    forest_plot <- summ_tbl %>% 
      ggplot(aes(x = estimate, 
                 y = reorder(parameter, estimate), 
                 color = corr_sign)) + 
      geom_vline(xintercept = 1, 
                 linetype = 'dashed') + 
      geom_errorbarh(aes(xmin = lower_ci, 
                         xmax = upper_ci), 
                     height = 0) + 
      geom_point(shape = 16, 
                 size = 2) + 
      geom_text(aes(label = est_lab), 
                size = 2.5, 
                hjust = 0.5, 
                vjust = -0.8) + 
      scale_x_continuous(trans = x_trans) + 
      scale_y_discrete(labels = y_scale) + 
      scale_color_manual(values = globals$pos_neg_scale) + 
      guides(color = F) + 
      globals$common_theme + 
      theme(axis.title.y = element_blank(), 
            panel.grid.major = element_line(color = 'gray90')) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = 'HR')
    
    return(forest_plot)
    
  }

# ROC plotting and analysis wrappers -----
  
  plot_com_roc <- function(score_tbl, 
                           stat_tbl, 
                           signature, 
                           status_var = 'death_acute', 
                           plot_title = NULL, 
                           plot_subtitle = NULL, 
                           plot_tag = NULL, 
                           show_auc = T) {
    
    ## plots the ROC for the given signature
    
    ## meta
    
    score_tbl <- score_tbl %>% 
      select(.data[[status_var]], 
             all_of(c(signature, 
                      pah_study$comparators$variable))) %>% 
      gather(key = 'model_id', 
             value = 'score', 
             all_of(c(signature, 
                      pah_study$comparators$variable))) %>% 
      mutate(model_id = factor(model_id, 
                               c(signature, 
                                 pah_study$comparators$variable)))
    
    stat_tbl <- stat_tbl %>% 
      filter(model_id %in% c(signature, 
                             pah_study$comparators$variable)) %>% 
      mutate(model_lab = ifelse(model_id %in% pah_study$comparators$variable, 
                                globals$comp_labs[model_id], 
                                stri_replace(model_id, fixed = 'sign_', replacement = 'Sign ')), 
             AUC_lab = paste0(signif(AUC, 2), 
                              ' [', 
                              signif(lowerCI, 2), 
                              ' - ', 
                              signif(upperCI, 2), 
                              ']'))
    
    if(show_auc) {
      
      stat_tbl <- stat_tbl %>% 
        mutate(model_lab = paste(model_lab, 
                                 AUC_lab, sep = ', '))
      
    }
    
    model_labs <- set_names(stat_tbl$model_lab, 
                            stat_tbl$model_id)
    
    model_colors <- c(globals$comp_colors, 
                      set_names('black', signature))
    
    ## plot
    
    roc <- plot_roc(inp_table = score_tbl, 
                    m_variable = 'score', 
                    d_variable = status_var, 
                    marker_variable = 'model_id', 
                    labels = F, 
                    pointsize = NA) + 
      geom_abline(slope = 1, 
                  intercept = 0, 
                  linetype = 'dashed', 
                  color = 'gray50') + 
      scale_color_manual(values = model_colors, 
                         labels = model_labs, 
                         name = '') + 
      theme(plot.title = element_text(size = 8, face = 'bold'), 
            plot.subtitle = globals$common_text, 
            plot.tag = globals$common_text) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag)
    
    return(roc)
    
  }
  
  
  
# Clustering accessories -----
  
  plot_radial <- function(clust_stats, 
                          plot_title = NULL, 
                          plot_subtitle = NULL, 
                          plot_tag = NULL) {
    
    ## makes a base radial plot for the selected clustering features
    
    radial_plot <- clust_stats %>% 
      ggplot(aes(x = variable,
                 y = mean * 100, 
                 fill = variable)) +
      geom_bar(stat = 'identity', 
               color = 'black') + 
      coord_polar(theta = 'x') + 
      globals$common_theme + 
      theme(panel.grid.major = element_line(color = 'gray80'), 
            legend.position = 'bottom', 
            axis.line = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.y = element_text(size = 6)) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           y = '% maximum')
    
    return(radial_plot)
    
  }
  
# END ----