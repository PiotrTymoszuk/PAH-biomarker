# This script contains a set of tools used for Kaplan-Meier analysis of survival

# libraries ----

  library(plyr)
  library(tidyverse)
  library(purrr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggrepel)
  library(furrr)

# graphic theme -----

  common_theme <- theme_classic() + 
    theme(legend.text = element_text(size = 8, face = 'plain', color = 'black'), 
          legend.title = element_text(size = 8, face = 'plain', color = 'black'), 
          axis.text = element_text(size = 8, face = 'plain', color = 'black'), 
          axis.title = element_text(size = 8, face = 'plain', color = 'black'), 
          strip.text = element_text(size = 8, face = 'plain', color = 'black'), 
          strip.background = element_rect(color = 'black', fill = 'gray95'), 
          plot.tag = element_text(size = 8, face = 'plain', color = 'black', 
                                  margin = margin(8, 5, 5, 5), hjust = 0), 
          plot.tag.position = 'bottom', 
          plot.title = element_text(size = 8, face = 'plain', color = 'black'), 
          plot.subtitle = element_text(size = 8, face = 'plain', color = 'black'), 
          plot.margin = margin(t = 3, l = 3, r = 3, b = 3, unit = 'mm'))
    
# survival toolbox ----
    
    create_surv <- function(inp_table, time_variable, event_variable) {
      
      ## creates a survival object
      
      return(Surv(inp_table[[time_variable]], inp_table[[event_variable]]))
      
    }
    
    model_coxph <- function(inp_table, surv_object, indep_variable) {
      
      ## models time to the event for as a function of a dependent variable
      ## with a one parameter Cox model

      mod_form <- paste('surv_object', 
                        paste(indep_variable, collapse = '+'), 
                        sep = '~') %>% 
        as.formula
      
      cox_object <- coxph(formula = mod_form, 
                          data = inp_table)
      
      return(list(surv_object = surv_object, 
                  cox_object = cox_object, 
                  summary = summary(cox_object)))
      
    }
    
    model_km_cutpoint <- function(inp_table, surv_object, indep_variable, 
                                  cutpoint = NULL, strata_labs = NULL, ...) {
      
      ## performs KM modeling for the given table, survival object
      ## and an independent variable stratified by the given cutpoint and optionally custom-labeled.
      ## Alternatively a binary variable or factor may be supplied.
      ## vectors of independent variables/cutoffs are handled by recursion.
      ## ... specifies further arguments for the survdiff() function.
      ## Returns one or a list of 'km_model_object' instances with the following components:
      ## test_table: modeling table with the stratified independent variable of interest
      ## surv_object: survival object used in modeling
      ## km_fit: survfit object which may be used for plotting
      ## test: results of significance testing
      ## stats: statistics returned by survdiff()
      ## p_value: p_value of significance in survival differnces between strata of the independent value
      
      if(length(indep_variable) > 1 & length(cutpoint) > 1) {
        
        ## handling of vactors
        
        res_lst <- list(indep_variable = indep_variable, 
                        cutpoint = cutpoint) %>% 
          pmap(model_km_cutpoint, 
               inp_table = inp_table, 
               surv_object = surv_object, 
               strata_labs = strata_labs, 
               ...) %>% 
          set_names(indep_variable)
        
        return(res_lst)
        
      }
      
      modeling_table <- inp_table %>% 
        mutate(indep_variable = .[[indep_variable]])
      
      if(!is.null(cutpoint)) {
        
        min_var <- min(modeling_table$indep_variable, na.rm = T) - 1
        
        max_var <- max(modeling_table$indep_variable, na.rm = T)
        
        modeling_table <- modeling_table %>% 
          mutate(var_strata = cut(indep_variable, 
                                  c(min_var, cutpoint, max_var), 
                                  labels = strata_labs))
        
      } else {
        
        modeling_table <- modeling_table %>% 
          mutate(var_strata = indep_variable)
        
      }
      
      km_fit <- survfit(surv_object ~ var_strata, data = modeling_table)
      
      test <- survdiff(surv_object ~ var_strata, data = modeling_table, ...)
      
      stats <- list(chi = test$chisq, df = length(test$n) - 1, n = test$n)
      
      p_value <- 1 - pchisq(test$chisq, length(test$n) - 1)
      
      km_model_obj <- list(test_table = modeling_table, 
                           surv_object = surv_object, 
                           km_fit = km_fit, 
                           test = test, 
                           cutpoint = cutpoint, 
                           stats = stats, 
                           p_value = p_value)
      
      attr(km_model_obj, 'class') <- 'km_model_object'
      
      return(km_model_obj)
      
    }
    
    optimize_km <- function(inp_table, surv_object, indep_variable, step = 0.01, min_n = 0, 
                            return_first = T, detailed = F, .parallel = F, ...) {
      
      ## finds am optimal cutpoint for the given independent value with the largest significance
      ## for the survival difference between the strata
      ## step specifies the step of expression range screening, 
      ## min_n argument defines the minimal group size for a indep_variable strata
      ## return_first: return the first maximum (lowest cutpoint value) or all maxima
      ## detailed: should the p_values for all tested cutpoints be returned of just the maximum
      ## handles multiple independent variables by recursion, serial or parallel
      ## ... specifies additional arguments passed to the model_km_cutpoint() function
      
      if(length(indep_variable) > 1) {
        
        if(.parallel) {
          
          plan('multisession')
          
          res_lst <- indep_variable %>% 
            future_map(optimize_km, 
                       inp_table = inp_table, 
                       surv_object = surv_object, 
                       step = step, 
                       min_n = min_n, 
                       return_first = return_first, 
                       detailed = detailed, 
                       ...) %>% 
            set_names(indep_variable)
          
          plan('sequential')
      
        } else {
          
          res_lst <- indep_variable %>% 
            map(optimize_km, 
                inp_table = inp_table, 
                surv_object = surv_object, 
                step = step, 
                min_n = min_n, 
                return_first = return_first, 
                detailed = detailed, 
                ...) %>% 
            set_names(indep_variable)
          
        }
        
        return(res_lst)
        
      }
      
      end_cutoff <- max(inp_table[[indep_variable]], 
                        na.rm = T)
      
      min_cutoff <- min(inp_table[[indep_variable]], 
                            na.rm = T)
      
      cutoff_seq <- seq(min_cutoff, end_cutoff, by = step)
      
      km_models <- cutoff_seq %>% 
        map(function(x) try(model_km_cutpoint(inp_table = inp_table, 
                                              surv_object = surv_object, 
                                              indep_variable = indep_variable, 
                                              cutpoint = x), 
                            silent = T)) ## needs to be error resistant, cause some cutpoints give strata with n = 0
      
      km_models <- km_models %>% ## removal of the errors from the resul list
        map(function(x) if(class(x) == 'km_model_object') x else NULL) %>% 
        compact
        
      result_frame <- km_models %>% 
        map_dfr(function(x) tibble(cutoff = x$cutpoint, 
                               p_value = x$p_value, 
                               n1 = x$stats$n[1], 
                               n2 = x$stats$n[2]))
      
      optimal_cutoff <- result_frame %>% 
        filter(n1 > min_n, 
               n2 > min_n) %>% 
        filter(p_value == min(p_value)) %>% 
        .$cutoff

      if(return_first) {
        
        optimal_cutoff <- optimal_cutoff[1]
        
      }
      
      if(detailed) {
        
        return(list(result_table = result_frame, 
                    optimal_cutoff = optimal_cutoff))
        
      } else {
        
        return(optimal_cutoff)
        
      }
      
    }
    
    km_optimization_plot <- function(optimize_km_results) {
      
      ## takes the detailed results of optimize_km() and returns a basic diagnostic plot
      ## to assess the quality of KM optimization: group sizes (n) and p value as a function
      ## of the cutpoint
      
      opt_plot <- optimize_km_results$result_table %>% 
        mutate(neg_log10_p_value = -log10(p_value)) %>% 
        gather(key = 'parameter', 
               value = 'par_value', 
               p_value, 
               n1, 
               n2, 
               neg_log10_p_value) %>% 
        ggplot(aes(x = cutoff, 
                   y = par_value)) + 
        geom_line() + 
        facet_grid(parameter ~ ., 
                   scales = 'free', 
                   labeller = labeller(.cols = as_labeller(c(n1 = 'n1',
                                                             n2 = 'n2', 
                                                             p_value = 'p value', 
                                                             neg_log10_p_value = '-log10 p value')))) + 
        theme(axis.title.y = element_blank())
      
      return(opt_plot)
      
    }
    
    km_summary <- function(km_model_object) {
      
      ## returns the summary of the km_model_object: stats and p value
      ## of the difference in survival between the independent variable strata
      ## handles lists of km_model_objects by recursion
      
      if(class(km_model_object) == 'list') {
        
        res_lst <- km_model_object %>% 
          map_dfr(km_summary) %>% 
          mutate(feature = names(km_model_object))
        
        return(res_lst)
        
      }
      
      output <- km_model_object$stats[c('chi', 'df')] %>% 
        reduce(cbind)
      
      output <- km_model_object$stats$n %>% 
        reduce(cbind) %>% 
        cbind(output, .)
      
      output <- cbind(output, km_model_object$p_value)
      
      if(!is.null(km_model_object$cutpoint)) {
        
        output <- cbind(output, km_model_object$cutpoint) %>% 
          as.data.frame
        
        names(output) <- c('ChiSq', 
                           'Df', 
                           paste('N', 
                                 1:length(km_model_object$stats$n), 
                                 sep = ''), 
                           'P_value', 
                           'Cutoff')
        
      } else {
        
        output <- output %>% 
          as_tibble
        
        names(output) <- c('ChiSq', 
                           'Df', 
                           paste('N', 
                                 1:length(km_model_object$stats$n), 
                                 sep = ''), 
                           'P_value')
        
      }
      
      return(output %>% 
               as_tibble)
      
    }
    
    plot_km <- function(km_model_object, p_value = NULL, plot_title = NULL, x_lab = 'Time', ...) {
      
      # plots the km model of interest with the ggsurvplot function
      # ... are additional arguments for the ggsurvplot function
      # handles multiple km_model_objects by recursion
      
      if(class(km_model_object) == 'list') {
        
        res_lst <- km_model_object %>% 
          map(plot_km, ...) %>% 
          set_names(names(km_model_object))
        
        return(res_lst)
        
      }
      
      surv_object <- km_model_object$surv_object

      fit <- surv_fit(surv_object ~ var_strata, data = km_model_object$test_table)
      
      if(is.null(p_value)) {
        
        p_value <- signif(km_model_object$p_value, 2)
        
      }
      
      if(is.null(km_model_object$cutpoint)) {
        
        plot_subtitle <- NULL
        
      } else {
        
        plot_subtitle <- paste('Cutoff = ', signif(km_model_object$cutpoint, 2))
        
      }

      
      output <- ggsurvplot(fit = fit, 
                           ggtheme = globals$common_theme, 
                           pval = p_value, 
                           pval.size = 2.75, ...)  + 
        labs(title = plot_title, 
             subtitle = plot_subtitle, 
             tag = paste('\nLow: n = ', 
                         km_model_object$stats$n[1], 
                         ', High: n = ', 
                         km_model_object$stats$n[2], 
                         sep = ''), 
             x = x_lab)
        
      
      return(output$plot)
      
    }
    
    save_km_plot <- function(km_plot, file_name, width, height) {
      
      ## saves a km plot on the disc together with the risk table
      
      cairo_pdf(file_name, 
                width = cm_to_inch(width), 
                height = cm_to_inch(height))
      
      print(km_plot, 
            newpage = FALSE)
      
      dev.off()
      
    }

# varia ----
    
    cm_to_inch <- function(inp_cm){
      
      return(0.393701 * inp_cm)
      
    }
    
# END ---

