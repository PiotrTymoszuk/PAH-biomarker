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
      ddply(.(ID), summarise, death_study = if(sum(dead_alive, na.rm = TRUE) > 0) 1 else 0) %>% 
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
  
  translate_vars <- function(variable_vector, value = 'label', lexicon = pah_study$legend) {
    
    ## translates the variable vector into a vector of corresponding labels
    
    transl_vec <- lexicon[[value]] %>% 
      set_names(lexicon[['variable']])
    
    transl_vec <- transl_vec[variable_vector]
    
    transl_vec[is.na(transl_vec)] <- ''
    
    return(transl_vec)
    
  }
  
# wrappers for Cox modeling and model validation -----
  
  model_spline_cox <- function(data, event_variable, time_variable, indep_variable, sec_order = TRUE, ...) {
    
    ## model formula
    
    surv_chunk <- paste0('Surv(', time_variable, ', ', event_variable, ')')
    
    if(is.numeric(data[[indep_variable]])) {
      
      if(sec_order) {
        
        indep_chunk <- paste0(indep_variable, ' + ', 'I(', indep_variable, '^2)')
        
      } else {
        
        indep_chunk <- indep_variable
        
      }
      
      num_index <- TRUE
      
    } else {
      
      indep_chunk <- indep_variable
      
      num_index <- FALSE
      
    }
    
    mod_form <- paste(surv_chunk, indep_chunk, sep = '~') %>% 
      as.formula
    
    ## model
    
    cox_call <- call2('coxph', formula = mod_form, data = data)
    
    cox_mod <- eval(cox_call)
    
    ## model summary
    
    cox_summary <- get_cox_results(cox_mod, 
                                   exponentiate = TRUE, 
                                   r_stats = TRUE)
    
    c_obj <- concordance(cox_mod)
    
    cox_summary <- cox_summary %>% 
      mutate(variable = indep_variable, 
             level = if(!num_index) stri_replace(parameter, fixed = variable, replacement = '') else NA, 
             order = if(num_index & sec_order) c(1, 2) else NA, 
             aic = AIC(cox_mod), 
             c_index = c_obj$concordance, 
             c_lower_ci = c_index + sqrt(c_obj$var) * qnorm(0.025), 
             c_upper_ci = c_index + sqrt(c_obj$var) * qnorm(0.975))
    
    cox_summary <- cox_summary[c('parameter', 
                                 'variable', 
                                 'level', 
                                 'order', 
                                 'estimate', 
                                 'se', 
                                 'z', 
                                 'lower_ci', 
                                 'upper_ci', 
                                 'p_value', 
                                 'n_complete', 
                                 'rsq_mev', 
                                 'aic', 
                                 'c_index', 
                                 'c_lower_ci', 
                                 'c_upper_ci')]
    
    ## assumptions
    
    cox_assum <- cox.zph(cox_mod)
    
    ## output
    
    list(model = cox_mod, 
         summary = as_tibble(cox_summary), 
         assumptions = cox_assum)
    
  }
  
  
# displaying the modeling results -----
  
  plot_summ_forest <- function(inp_tbl, 
                               variable = 'variable', 
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
                               color_lab = 'Risk\ncorrelation', 
                               cutpoint = 1, 
                               signif_digits = 2, 
                               facet = T) {
    
    ## plots a two-cohort Forest plot
    
    ## meta
    
    grid_formula <- expr(!!ensym(variable) ~ !!ensym(cohort_var)) %>% 
      as.formula
    
    ## plot
    
    inp_tbl <- inp_tbl %>% 
      mutate(axis_lab = ifelse(.data[[level]] %in% c('no', 'yes', '') | is.na(.data[[level]]), 
                               translate_vars(.data[[variable]]), 
                               paste(translate_vars(.data[[variable]]), .data[[level]], sep = ': ')), 
             axis_lab = ifelse(order == 1 | is.na(order), axis_lab, paste0(axis_lab, '\u00B2')), 
             axis_lab = stri_replace(axis_lab, fixed = 'cm2', replacement = 'cm\u00B2'), 
             axis_lab = stri_replace(axis_lab, fixed = '=', replacement = '\u2265'), 
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

# clustering graphics -----
  
  
  plot_vio_panel <- function(clust_object, 
                             plot_title = NULL, 
                             plot_subtitle = NULL, 
                             x_lab = 'Median-normalized value') {
    
    ## plots a violin panel of the clustering features
    
    ## plotting table in the long format
    
    plotting_tbl <- clust_object %>% 
      extract('data')
    
    plotting_ft <- names(plotting_tbl)
    
    plotting_tbl <- plotting_tbl %>% 
      rownames_to_column('ID') %>% 
      gather(key = 'feature', 
             value = 'value', 
             all_of(plotting_ft))
    
    plotting_tbl <- clust_object %>% 
      extract('assignment') %>% 
      set_names(c('ID', 'clust_id')) %>% 
      left_join(plotting_tbl, ., by = 'ID') %>% 
      mutate(feature = translate_vars(feature))
    
    ## n numbers, presented in the plot tag
    
    plot_tag <- ngroups(clust_object)
    
    plot_tag <- map2_chr(plot_tag[[1]], 
                         plot_tag[[2]], 
                         ~paste0(.x, ': n = ', .y)) %>% 
      paste(collapse = ', ') %>% 
      paste0('\n', .)
    
    ## violin plot panel
    
    plotting_tbl %>% 
      ggplot(aes(x = value, 
                 y = feature, 
                 fill = clust_id)) + 
      geom_violin(scale = 'width', 
                  alpha = 0.25, 
                  show.legend = FALSE) + 
      geom_point(aes(color = clust_id), 
                 size = 2, 
                 shape = 16, 
                 position = position_jitterdodge(jitter.width = 0.1, 
                                                 dodge.width = 0.9), 
                 alpha = 0.5) + 
      scale_fill_manual(values = globals$cluster_colors, 
                        name = 'Cluster') + 
      scale_color_manual(values = globals$cluster_colors, 
                         name = 'Cluster') +
      globals$common_theme + 
      theme(axis.title.y = element_blank(), 
            strip.background = element_blank(), 
            strip.text = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = x_lab)
    
  }
  
# varia ------
  
  set_common_y <- function(plot_IBK, plot_LZ) {
    
    ## sets the common y numeric scale for two plots
    ## with the clustering features
    
    cmm_limits <- map(list(plot_IBK, plot_LZ), 
                      ~.x$data$variable) %>% 
      reduce(c) %>% 
      range

    plot_IBK <- plot_IBK + 
      labs(y = paste0('IBK, ', plot_IBK$labels$y))
    
    plot_LZ <- plot_LZ + 
      labs(y = paste0('LZ/W, ', plot_LZ$labels$y), 
           title = '')
    
    list(plot_IBK, plot_LZ) %>% 
      map(~.x + scale_y_continuous(limits = cmm_limits))
    
  }
  
  mm_inch <- function(input_mm) {
    
    return(0.0393700787 * input_mm)
    
  }
  
# END ----