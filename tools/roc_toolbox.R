# A set of useful functions for ROC plotting and analysis

# libraries -----

  require(tidyverse)
  require(stringi)

  library(OptimalCutpoints)
  library(plotROC)

# globals ----

  ## a common theme for ROC
  
    roc_theme <-  theme(legend.text = element_text(size = 8, face = 'plain', color = 'black'), 
                        legend.title = element_text(size = 8, face = 'plain', color = 'black'), 
                        axis.text = element_text(size = 8, face = 'plain', color = 'black'), 
                        axis.title = element_text(size = 8, face = 'plain', color = 'black'), 
                        plot.title = element_text(size = 10, face = 'bold'), 
                        plot.subtitle = element_text(size = 8, face = 'plain', color = 'black'), 
                        plot.tag = element_text(size = 8, face = 'plain', color = 'black', hjust = 0), 
                        plot.tag.position = 'bottom')
    
    
# ROC analysis and plotting toolbox ----
    
    find_optimal_cutoff <- function(inp_table, status_variable, marker, tag.healthy = 0, method = 'Youden', ...) {
      
      ## returns an output of the optimal.cutpoints from the package OptimalCutpoints
      
      analysis_table <- inp_table %>% 
        data.frame
      
      output <- optimal.cutpoints(status = status_variable, X = marker, data = analysis_table, 
                                  tag.healthy = tag.healthy, method = 'Youden', ...)
      
      return(output)
      
    }
    
    get_optimal_cutpoint_coefs <- function(optimal_cutpoint_object) {
      
      ## extracts coefficient values from the optimal_cutpoint_object generated e.g. by 
      ## the find_optimal_cutoff function
      
      output <- optimal_cutpoint_object %>% 
        summary
      
      output <- output$p.table$Global$Youden
      
      return(output)
      
    }
    
    get_auc <- function(optimal_cutpoint_object) {
      
      require(stringi)
      
      ## extracts AUC with its confidence intetervals from the optimal_cutpoint_object generated e.g. by 
      ## the find_optimal_cutoff function
      
      output <- optimal_cutpoint_object %>% 
        summary
      
      output <- output$p.table$Global$AUC_CI
      
      output <- output %>% 
        stri_replace(., fixed = ",", replacement = '') %>% 
        stri_replace(., fixed = '(', replacement = '') %>% 
        stri_replace(., fixed = ')', replacement = '') %>%
        stri_split(., fixed = ' ', simplify = T) %>% 
        as.numeric %>% 
        c
      names(output) <- c('AUC', 'lowerCI', 'upperCI')
      
      return(output)
      
    }
    
    plot_roc <- function(inp_table, m_variable, d_variable, marker_variable = NULL, ...) {
      
      ## draws a plot with ROC curves. ROCs for diverse markers/parameter defined by the marker_variable
      ## are plotted with separate colors. M variable defines the value of the marker, d codes for the disease
      ## status
      
      if(is.null(marker_variable)) {
        
        output <- inp_table %>% 
          ggplot(aes(m = .data[[m_variable]], d = .data[[d_variable]]))
        
      } else {
        
        output <- inp_table %>% 
          ggplot(aes(m = .data[[m_variable]], d = .data[[d_variable]], color = .data[[marker_variable]]))
        
      }
      
      output <- output +  
        geom_roc(...) + 
        style_roc() + 
        roc_theme
      
      return(output)
      
    }
    
# END -----