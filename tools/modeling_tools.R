# This script provides functional tools for modeling tasks in christas Sepsis project

# libraries -----

  library(plyr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lme4)
  library(lmerTest)
  library(stringi)

# globals ----

  default_signif = 2

# simple model generation functions -----

  make_lmer <- function(inp_table, response = NULL, fix_variables = NULL, random_variables = NULL, interaction_term = '', ...) {
    
    # generaes a simple no-interaction lm (if random variables not given) or lmer
    # Interaction can be included by explicitly stating the interaction term
    
    inp_formula <- paste(response, '~', paste(fix_variables, collapse = '+'), interaction_term, sep = '+')
    
    if(is.null(random_variables)) {
      
      out_model <- lm(formula = as.formula(inp_formula), data = inp_table, ...)
      
    } else {
      
      for(var in random_variables) {
        
        inp_formula <- paste(inp_formula, '+ (1|', var, ')')
        
      }
      
      out_model <- lmer(formula = as.formula(inp_formula), data = inp_table, ...)
      
    }
    
    return(out_model)
    
  }
  
  make_glm <- function(inp_table, response = NULL, fix_variables = NULL, ...) {
    
    # generaes a simple no-interaction glm
    
    inp_formula <- paste(response, '~', paste(fix_variables, collapse = '+'))
    
    return(glm(formula = as.formula(inp_formula), data = inp_table, ...))
    
  }
  
# model summary functions ----

  get_lme_results <- function(model_mix, L_matrix = NULL, anova = F, ...){
    
    require(lmerTest)
    
    ## a function returning confidence intervals and significances for particular coefficents of the mixed model with one fixed
    ## effect. The coefficients to display are defined by the L_matrix. Uses the contest function from the lmerTest package
    ## If the L_matrix is not given, stats for all coefficients are calculated
    
    output <- model_mix %>% 
      as_lmerModLmerTest
    
    if(is.null(L_matrix)) {
      
      L_matrix <- identity_matrix(length(model_mix@beta))
      
    }
    
    if(anova) {
      
      output <- contest(output, L_matrix, ...)
      
    } else {
      
      output <- contest(output, L_matrix, joint = F, ...) %>% 
        mutate(P_BH = p.adjust(`Pr(>|t|)`, method = 'BH')) %>% 
        mutate(CI_error_bar = (upper - lower)/2)
      as_tibble ## multiple testing adjustment by Benjamini-Hochberg
      ## returning explicit error CI error bars for plot annotation
      
      output <- output %>% 
        mutate(plot_label = paste(signif(Estimate, 2), '\u00B1', signif(CI_error_bar, 2), 
                                  ', p = ', signif(P_BH, 2))) ## ready to use significance labels
      
      coef_names <- rownames(summary(model_mix)$coefficients) %>% 
        as_tibble
      
      output <- cbind(coef_names, output)
      
    }
    
    return(output)
  }
  
  save_lme_anova_coef_tests <- function(model_mix, file_prefix = NA, ...) {
    
    ## a wrapper for the get_lme_results and anova function for a given lme object and contrast matrix
    ## returns an list with the output of the anova (as described in lmerTest) and the get_lme_results
    ## is file_prefix is given, results are saved on the disc
    
    anova_results <- anova(model_mix, ...)
    
    if (class(model_mix) == 'lm') {
      
      lme_results <- ext_summary_fixed(model_mix)
      
    } else {
      
      lme_results <- get_lme_results(model_mix, ...)
      
    }
    
    if(!is.na(file_prefix)) {
      
      write_tsv(anova_results, paste(file_prefix, '_anova.txt', sep = ''))
      
      write_tsv(lme_results, paste(file_prefix, '_coefs.txt', sep = ''))
      
    }
    
    return(list(anova = anova_results, coefs = lme_results))
    
  }
  
  ext_summary_fixed <- function(lm_model, correction = 'BH', ...) {
    
    ## the output is an extended summary for lm coefficients.
    ## in addition to the summary output, 95% CI is calculated
    ## and the p values are corrected with the p.adjust function
    ## (the correction method is specified by the correction argument)
    ## ... stands for additional parameters passed to the confint function
    
    output <- summary(lm_model)$coefficients %>% 
      as_tibble()
    
    output$p_adjusted <- p.adjust(output[[ncol(output)]], method = correction)
    
    ci <- confint(lm_model, ...)
    
    output <- cbind(output, ci)
    
    output$CI <- (output[['97.5 %']] - output[['2.5 %']])/2
    
    output$plot_label = paste(signif(output$Estimate, default_signif), 
                              ' \u00B1 ', 
                              signif(output$CI, default_signif), 
                              ', p = ', 
                              signif(output$p_adjusted, default_signif), sep = '')
    output$parameter <- rownames(output)
    
    return(output %>% set_names(c('estimate', 'se', 't', 'p_value', 'p_adjusted', 'lower_ci', 'upper_ci', 'ci', 
                                  'plot_label', 'paramater')))
    
  }
  
  get_glm_results <- function(glm_model, exponentiate = F) {
    
    # retrieves regression estimate stats from a glm
    
    ci <- glm_model %>% 
      confint %>% 
      as_tibble %>% 
      set_names(c('lower_ci', 'upper_ci'))
    
    coefs <- glm_model %>% 
      summary %>% 
      coefficients %>% 
      data.frame()
    
    estimates <- coefs[, c('Estimate', 'Std..Error')] %>% 
      as_tibble %>% 
      set_names(c('estimate', 'se'))
    
    parameters <- data.frame(parameter = rownames(coefs))

    stats <- coefs[, c('z.value', 'Pr...z..')] %>% 
      as_tibble %>% 
      set_names(c('z', 'p_value'))
    
    if(any(class(glm_model) == 'glmerMod')) {
      
      ci <- ci[-1, ] ## dropping out the CI estimates for sigma
      n_number <- nrow(glm_model@frame)
      
    } else {
      
      n_number <- nrow(glm_model$model)
      
    }
    
    if(!exponentiate) {
      
      summary_table <- cbind(parameters, estimates, ci, stats)
      
    } else {
      
      summary_table <- cbind(parameters, exp(estimates), exp(ci), stats)
      
    }
    
    summary_table <- summary_table %>% 
      mutate(n_number = n_number, 
             plot_label = paste(signif(estimate, default_signif), 
                                ' [', 
                                signif(lower_ci, default_signif), 
                                ' - ', 
                                signif(upper_ci, default_signif), 
                                '], p = ', 
                                signif(p_value, default_signif), sep = ''))
    
    return(summary_table)
    
    
  }
  
  save_glm_tests <- function(model_mix, file_prefix = NA, ...) {
    
    write_tsv(get_glm_results(glm_model = model_mix, ...), paste(file_prefix, '_glm.txt', sep = ''))
    
  }
  
  make_anova_sum_table <- function(model_list) {
    
    ## for a given model_fit, a summary table with Effect names, F and p values is returned
    
    sum_tables <- model_list %>% 
      map(anova) %>% 
      map(function(x) mutate(x, Effect = rownames(x))) %>% 
      map(as_tibble)
    
    sum_tables <- sum_tables %>% 
      map(function(x) mutate(x, ANOVA = paste('F = ', 
                                              signif(`F value`, default_signif), 
                                              ' (', 
                                              signif(NumDF, default_signif), 
                                              ', ', 
                                              signif(DenDF, default_signif), 
                                              ')\n', 
                                              'p = ', 
                                              signif(`Pr(>F)`, default_signif), 
                                              sep = ''))) %>% 
      map(function(x) select(x, Effect, ANOVA))
    
    return(sum_tables)
    
  }  
  
  get_cox_results <- function(cox_model, exponentiate = F, r_stats = F){
    
    ci <- cox_model %>% 
      confint %>% 
      as_tibble %>% 
      set_names(c('lower_ci', 'upper_ci'))
    
    coefs <- cox_model %>% 
      summary %>% 
      coefficients %>% 
      data.frame()
    
    estimates <- coefs[, c('coef', 'se.coef.')] %>% 
      as_tibble %>% 
      set_names(c('estimate', 'se'))
    
    parameters <- data.frame(parameter = rownames(coefs))
    
    stats <- coefs[, c('z', 'Pr...z..')] %>% 
      as_tibble %>% 
      set_names(c('z', 'p_value'))
    
    if(!exponentiate) {
      
      summary_table <- cbind(parameters, estimates, ci, stats)
      
    } else {
      
      summary_table <- cbind(parameters, exp(estimates), exp(ci), stats)
      
    }
    
    summary_table <- summary_table %>% 
      mutate(plot_label = paste(signif(estimate, default_signif), 
                                ' (', 
                                signif(lower_ci, default_signif), 
                                ', ', 
                                signif(upper_ci, default_signif), 
                                '), p = ', 
                                signif(p_value, default_signif), sep = ''), 
             n_complete = cox_model$n)
    
    if(r_stats) {
      
      ## adds some of the R square stats returned by survMisc::rsq
      
      require(survMisc)
      
      r_stat_tbl <- cox_model %>% 
        rsq
      
      r_stat_tbl <- r_stat_tbl %>% 
        reduce(cbind) %>% 
        as.data.frame %>% 
        set_names(paste('rsq', names(r_stat_tbl), sep = '_')) %>% 
        as_tibble
      
      summary_table <- cbind(summary_table, r_stat_tbl)
      
    }
    
    
    return(summary_table)
    
  }
  
# multiple testing ----
  
  t_tester <- function(inp_table, response, test_variable, 
                       level1, level2, 
                       test_fun = function(x, y) t.test(x, y)$p.value) {
    
    ## a function testing differences in response value in data subsets defined by
    ## level1 and level2 of the test variable
    
    vector1 <- inp_table %>% 
      filter(.[[test_variable]] == level1) %>% 
      select(response) %>% 
      unlist
    
    vector2 <- inp_table %>% 
      filter(.[[test_variable]] == level2) %>% 
      select(response) %>% 
      unlist
    
    return(test_fun(vector1, vector2))
    
  }
  
  test_multiple <- function(inp_table, response, test_variable, 
                            test_fun = function(x, y) t.test(x, y)$p.value, 
                            adjust.method = NULL, round_signif = default_signif) {
    
    ## performs multiple testing for the given response stratified by levels of the test variable
    
    t_test_levels <- levels(factor(inp_table[[test_variable]]))
    
    ### result container
    
    tester_matrix <- matrix(NA, nrow = length(t_test_levels), ncol = length(t_test_levels)) %>% 
      set_rownames(t_test_levels) %>% 
      set_colnames(t_test_levels)
    
    ### testing loop
    
    for(row_index in rownames(tester_matrix)) {
      
      for(col_index in colnames(tester_matrix)) {
        
        tester_matrix[row_index, col_index] <- t_tester(inp_table = inp_table, 
                                                        response = response, 
                                                        test_variable = test_variable, 
                                                        test_fun = test_fun, 
                                                        level1 = row_index, 
                                                        level2 = col_index)
        
      }
    }
    
    ### adjusting for multiple testing
    
    if(!is.null(adjust.method)) {
      
      adj_vector <- c()
      
      for(row_index in rownames(tester_matrix)) {
        
        adj_vector <- c(adj_vector, tester_matrix[row_index, ])
        
      }
      
      adj_vector <- p.adjust(adj_vector, method = adjust.method)
      
      if(!is.null(round_signif)) {
        
        adj_vector <- signif(adj_vector, round_signif)
        
      }
      
      adj_matrix <- matrix(data = adj_vector, 
                           nrow = length(t_test_levels), 
                           ncol = length(t_test_levels), 
                           byrow = T) %>% 
        set_rownames(rownames(tester_matrix)) %>% 
        set_colnames(colnames(tester_matrix))
      
      return(adj_matrix)
      
    }
    
    return(tester_matrix)
    
  }
  
# n number determination ----
  
  get_group_n_numbers <- function(inp_table, grouping_variable, response) {
    
    n_number_info <- inp_table %>% 
      filter(!is.na(.[[grouping_variable]]), 
             !is.na(.[[response]])) %>% 
      dlply(grouping_variable) %>% 
      map_dfc(function(x) nrow(x)) %>% 
      set_names(levels(factor(inp_table[[grouping_variable]]))) %>% 
      mutate(parameter = 'n')
    
    return(n_number_info)
    
  }
  
  get_n_numbers <- function(inp_table, variables) {
    
    ## gets n numbers (non-NA filtered variables) from a table
    
    out_table <- inp_table
    
    for(var in variables) {
      
      out_table <- filter(out_table, !is.na(out_table[[var]]))
      
    }
    
    return(nrow(out_table))
    
  }
  
# prediction tools ----
  
  predict_glm <- function(glm_model, update_table = NULL, type = 'response', ...) {
    
    ## predicts the outcome of the glm model in the response scale.
    ## If an update tabel is provided, the predicted are merged with it
    
    if(is.null(update_table)) {
      
      return(predict(glm_model, type = type, ...))
      
    }
    
    pred_vals <- predict(glm_model, type = type, newdata = update_table, ...)
    
    return(cbind(update_table, pred_vals) %>% as_tibble)
    
  }
  
  calculate_exp_score <- function(){}
  
# significance testing -----
  
  identify_significant <- function(inp_table, param_variable, p_val_variable, est_variable, 
                                   alpha = 0.05, est_limit = 1, update_tbl = T) {
    
    ## identifies significantly regulated paramaters in the given table.
    ## param_variable specifies a column with paramater names, 
    ## p_val_variable with p_values, est_variable with estimate variables
    ## Alpha specifies significance limit, est_limit the middle point for regulation sign
    ## Update: updates the given inp_table
    
    out_table <- inp_table[, c(param_variable, p_val_variable, est_variable)]
    
    out_table <- out_table %>% 
      mutate(regulation = ifelse(.[[est_variable]] < est_limit, 'down', 'up')) %>% 
      mutate(regulation = ifelse(.[[p_val_variable]] < alpha, regulation, 'ns'))
    
    if(update_tbl) {
      
      return(left_join(inp_table, 
                       out_table[, c(param_variable, 'regulation')], 
                       by = param_variable))
      
    } else {
      
      return(out_table)
      
    }
    
  }
  
# varia ----
  
  identity_matrix <- function(dimension) {
    
    ## creates an identity matrix with a given dimension
    
    out_matrix <- matrix(0, dimension, dimension)
    
    for (i in 1:dimension) {
      
      out_matrix[i, i] <- 1
      
    }
    
    return(out_matrix)
    
  }
  
# END ----