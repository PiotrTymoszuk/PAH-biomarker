# This script provides tools for summaries and QC of linear (lm) and generalized linear models

# libraries ------

  require(tidyverse)
  require(ggrepel)
  require(broom)
  require(furrr)
  require(stringi)

# auxiliary functions ----

  point_plot_ <- function(data, x_var, y_var, x_lab = x_var, y_lab = y_var, 
                          plot_title = NULL, smooth = T, silent = T, ...) {
    
    ## draws a simple point plot for diagnostic purposes, takes the output of get_qc_tbl() as data argument
    ## color-codes model missfits
    
    ## table for plotting 
    
    data <- data %>% 
      mutate(misslab = ifelse(.candidate_missfit == 'yes',
                              .rownames, 
                              NA))
    
    ## fill colors
    
    fill_colors <- c(no = 'cornflowerblue', 
                     yes = 'firebrick4')
    
    ## point plot
    
    point_plot <- data %>% 
      ggplot(aes(x = .data[[x_var]], 
                 y = .data[[y_var]], 
                 fill = .candidate_missfit)) + 
      geom_point(size = 2, 
                 shape = 21) + 
      geom_text_repel(aes(label = misslab), 
                      show.legend = F) + 
      scale_fill_manual(values = fill_colors, 
                        name = 'Candidate missfit') + 
      labs(x = x_lab, 
           y = y_lab, 
           title = plot_title)
    
    if(smooth) {
      
      if(silent) {
        
        suppressWarnings(point_plot <- point_plot + 
                           geom_smooth(show.legend = F, 
                                       color = 'black', 
                                       fill = 'dodgerblue2', ...))
        
      } else {
        
        point_plot <- point_plot + 
          geom_smooth(show.legend = F, 
                      color = 'black', 
                      fill = 'dodgerblue2', ...)
        
      }
      
    }
      
    return(point_plot)
    
  }
  
  calc_expected_ <- function(inp_data, observed) {
    
    ## calculates expected normal distribution of a variable observed
    ## credits to: https://stackoverflow.com/questions/43217104/coloring-points-in-a-geom-qq-plot
    
    
    inp_data <- inp_data[order(inp_data[[observed]]), ] %>% 
      mutate(.expect.norm = qnorm(ppoints(nrow(.))))
  
    return(inp_data)
    
  }

# model summary -----

  get_estimates <- function(linear_model, transf_fun = NULL, fallback = F, silent_messages = T, ...) {
    
    ## extract coefficients from a linear or generalized linear model (betas)
    ## together with confidence intervals. The argument trans_fun allows for 
    ## transformation of the coefficients and CI (e.g. to convert them to OR in logistic regression)
    ## ... specifies other arguments to the CI-calculating function confint()
    ## the fallback option: calculates the CI based on the normality assumption: i. e. SEM * critical norm distr
    
    ## transforming function
    
    if(is.null(transf_fun)) {
      
      transf_fun <- function(x) x
      
    }
    
    ## model summary: to get se and p values
    
    mod_summary <- summary(linear_model)
    
    ## model estimates, error and CI, transforming
    
    model_coefs <- coefficients(linear_model)
    
    model_se <- mod_summary$coefficients[, 2]
    
    if(fallback) {
      
      model_ci <- tibble(lower_ci = model_coefs + qnorm(0.025)*model_se, 
                         upper_ci = model_coefs + qnorm(0.975)*model_se)
      
      
    } else {
      
      if(silent_messages) {
        
        model_ci <- suppressMessages(confint(linear_model, ...))

        
      } else {
        
        model_ci <- confint(linear_model, ...)
        
      }
      
      model_ci <- model_ci %>% 
        as_tibble %>% 
        set_names(c('lower_ci', 
                    'upper_ci'))
      
    }
    
    #return(model_ci)
    
    est_tibble <- model_ci %>% 
      mutate(estimate = unname(model_coefs)) %>% 
      map_dfc(transf_fun) %>% 
      mutate(parameter = names(model_coefs), 
             se = unname(model_se)) %>% 
      select(parameter, 
             estimate, 
             se, 
             lower_ci, 
             upper_ci)
    
    ## p values extracted from model summary
    ## n number of complete observations extracted from the model frame
    
    model_p <- mod_summary$coefficients[, 4]
    
    est_tibble <- est_tibble %>% 
      mutate(p_value = unname(model_p), 
             n_complete = nrow(model.frame(linear_model)))

    return(est_tibble)
    
  }
  
# model QC -----
  
  get_qc_tbl <- function(linear_model, ...) {
    
    ## generates a table with residuals of the model by calling 
    ## augment() from package broom. Finds potential model missfits
    ## by calculating normal distr. p for .std.resid (standardized residuals)
    ## 95% confidence region for the .fitted is estimated with normal distribution
    ## ... are further arguments to augment() function provided by broom
    ## handles both the ready-to-use analysis objects and linear model objects
    
    if(any(class(linear_model) == 'lm_analysis')) {
      
      linear_model <- linear_model$model
      
    }
    
    qc_tbl <- augment(linear_model, 
                      se_fit = T, ...)
    
    qc_tbl <- qc_tbl %>% 
      mutate(.sq.std.resid = .std.resid^2, 
             .lower_ci.fit = .fitted + .se.fit * qnorm(0.025), 
             .upper_ci.fit = .fitted + .se.fit * qnorm(0.975), 
             .candidate_missfit = ifelse(abs(.std.resid) > qnorm(0.975), 'yes', 'no'))
    
    ## adding a variable holding expected normal distribution for the standardized residuals
    
    qc_tbl <- calc_expected_(qc_tbl, '.std.resid')
    
    return(qc_tbl)
    
  }
  
  get_qc_plots <- function(linear_model, ...) {
    
    ## draws standard model qc plots with missfits labeled by obs. numbers
    ## ... are additional arguments passed to get_qc_tbl()
    ## handles both the ready-to-use analysis objects and linear model objects
    
    if(any(class(linear_model) == 'lm_analysis')) {
      
      linear_model <- linear_model$model
      
    }
    
    ## QC table
    
    qc_tbl <- get_qc_tbl(linear_model, ...)

    ## QC plots
    
    qc_plotting_lst <- list(x_var = c('.fitted', '.fitted', '.fitted', '.expect.norm', '.sigma'), 
                            y_var = c('.resid', '.std.resid', '.sq.std.resid', '.std.resid', '.cooksd'), 
                            plot_title = c('Residuals vs. fitted', 
                                           'Standardized residuals vs. fitted', 
                                           'Sqared residuals vs. fitted', 
                                           'QQ standardized residuals vs expected normal', 
                                           'Cook distance vs dropout sigma'),
                            method = c('loess', 'loess', 'loess', 'lm', 'loess'), 
                            smooth = c(T, T, T, T, F))
    
    qc_plots <- qc_plotting_lst %>% 
      pmap(point_plot_, 
           data = qc_tbl) %>% 
      set_names(c('resid_fitted', 
                  'std.resid_fitted', 
                  'sq.resid_fitted', 
                  'qq.std.resid', 
                  'cook_sigma'))
    

    
    return(qc_plots)
    
  }
  
# a combi linear model function - e.g. serial linear model generation ------
  
  make_lm_model <- function(data, response, indep_variable, weight_variable = NULL, mod_fun = glm, 
                            family = 'quasibinomial', est_transf = exp, confounder = NULL, 
                            .parallel = F, error_resistant = T, verbose = F, silent_messages = T, ...) {
    
    ## makes a lm or glm model with a summary and fit goodness measures
    ## for serial modeling: handles vectors of independent variables by recursion
    
    if(length(indep_variable) > 1) {
      
      start_time <- Sys.time()
      
      message(paste('Fitting', length(indep_variable), 'models'))
      
      if(.parallel) {
        
        plan('multisession')
        
        analysis_lst <- indep_variable %>% 
          future_map(make_lm_model, 
                     data = data, 
                     response = response, 
                     weight_variable = weight_variable, 
                     mod_fun = mod_fun, 
                     family = family, 
                     est_transf = est_transf, 
                     confounder = confounder, 
                     verbose = verbose, 
                     silent_messages = silent_messages, 
                     .options = furrr_options(packages = 'broom'), ...) %>% 
          set_names(indep_variable)
        
        plan(sequential)
        
      } else {
        
        analysis_lst <- indep_variable %>% 
          map(make_lm_model, 
              data = data, 
              response = response, 
              weight_variable = weight_variable, 
              mod_fun = mod_fun, 
              family = family, 
              est_transf = est_transf, 
              confounder = confounder, 
              verbose = verbose, 
              silent_messages = silent_messages, ...) %>% 
          set_names(indep_variable)
        
      }
      
      message(paste('Elapsed time:', Sys.time() - start_time))
      
      return(analysis_lst %>% 
               compact)
      
    }
    
    if(verbose) {
      
      message(paste('Modeling response:', 
                    response, 
                    ', indep. variable:', 
                    indep_variable))
      
    }
    
    ## model formula
    
    mod_formula <- paste(response, 
                         '~', 
                         paste(c(indep_variable, 
                                 confounder), 
                               collapse = '+')) %>% 
      as.formula
    
    ## weighting vector
    
    if(!is.null(weight_variable)) {
      
      w_vector <- data[[weight_variable]]
      
    } else {
      
      w_vector <- NULL
      
    }
    
    ## lm model
    
    if(error_resistant) {
      
      model <- try(mod_fun(formula = mod_formula, 
                           family = family, 
                           data = data, 
                           weights = w_vector),
                   silent = T)
      
      if(any(class(model) == 'try-error')) {
        
        return(NULL)
        
      }
      
    } else {
      
      model <- mod_fun(formula = mod_formula, 
                       family = family, 
                       data = data, 
                       weights = w_vector)
      
      
    }
    
    ## model summary
    
    if(error_resistant) {
      
      mod_summary <- try(model %>% 
                           get_estimates(transf_fun = est_transf), 
                         silent = T)
      
      if(any(class(mod_summary) == 'try-error')) {
        
        mod_summary <- try(model %>% 
                             get_estimates(transf_fun = est_transf, 
                                           fallback = T, 
                                           silent_messages = silent_messages), 
                           silent = T)
        
      }
      
    } else {
      
      mod_summary <- model %>% 
        get_estimates(transf_fun = est_transf, 
                      silent_messages = silent_messages)
      
    }
    
    mod_summary <- mod_summary %>% 
      mutate(response = response, 
             variable = ifelse(parameter == '(Intercept)', 
                               indep_variable, 
                               stri_extract(parameter, 
                                            fixed = indep_variable)), 
             level = ifelse(!is.na(variable), 
                            stri_replace(parameter, 
                                         fixed = indep_variable, 
                                         replacement  = ''), 
                            stri_replace(parameter, 
                                         fixed = confounder, 
                                         replacement  = '')),
             variable = ifelse(is.na(variable), 
                               confounder, 
                               indep_variable), 
             level = ifelse(level == '', 
                            NA, 
                            ifelse(level == '(Intercept)', 
                                   'baseline', 
                                   level)))
    
    mod_summary <- mod_summary %>% 
      select(all_of(c('response', 
                      'variable', 
                      'level', 
                      'parameter', 
                      'n_complete', 
                      'estimate', 
                      'se', 
                      'lower_ci', 
                      'upper_ci', 
                      'p_value')))
    
    ## model residual table
    
    mod_resid_tbl <- get_qc_tbl(model)
    
    ## model goodness measures: pseudo Rsq, AIC and MSE
    
    model_measures <- tibble(mse = mean(mod_resid_tbl$.resid^2), 
                             aic = AIC(model))
    
    if(class(model)[1] == 'lm') {
      
      model_measures$rsq <- summary(model)$r.squared
      
    } else {
      
      model_measures$rsq <- 1 - model$deviance/model$null.deviance
      
    }

    
    ## output
    
    out_list <- list(response = response, 
                     indep_variable = indep_variable, 
                     confounder = confounder, 
                     mod_fun = class(model)[1], 
                     family = family, 
                     model = model,
                     summary = mod_summary, 
                     measures = model_measures)
    
    attr(out_list, 'class') <- 'lm_analysis'
    
    return(out_list)
    
  }
  
# extractor functions -----
  
  get_model_summary <- function(lm_analysis, adj_method = 'BH') {
    
    ## extracts a model summary, handles lists of analyses by recursion
    ## for multiple objects: calculation of adjusted p
    
    if(all(class(lm_analysis) == 'list')) {
      
      summ_tbl <- lm_analysis %>% 
        map_dfr(get_model_summary) %>% 
        mutate(p_adj = p.adjust(p_value, method = adj_method))

      return(summ_tbl)
      
      
    }
    
    return(lm_analysis$summary)
    
  }
  
# END ----