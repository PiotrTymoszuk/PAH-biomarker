# This script contains tools for calculation of discrete and contonuous scores from a given set of variables

# libraries ----

  library(plyr)
  library(dplyr)
  library(purrr)

# discrete score calculation ----

  calculate_score <- function(inp_table, cutoff_list, score_list, stratify_only = F, tibble_output = F, show_var_cuts = F) {
    
    
    ## calculates a discrete-variable score fo the data in the given table. Variables
    ## to be used for score calculation are specified as names of cutoff_list and score_list.
    ## cutoff list contains values of the breaks the variable is cut into, 
    ## the score list contains the numbers to be assigned to the cut intervals of the variable
    ## The returned score variable/vector contains the sum of scores
    
    new_variable_names <- paste(names(cutoff_list), '_new', sep = '') %>% 
      set_names(names(cutoff_list))
    
    output_table <- inp_table
    
    for(var in names(cutoff_list)) {
      
      output_table[[new_variable_names[var]]] <- output_table[[var]] %>% 
        cut(breaks = cutoff_list[[var]], 
            labels = score_list[[var]])
      
    }
    
    if(stratify_only) {
      
      return(output_table)
      
    }
    
    output_table[['score']] <- 0
    
    for(var in new_variable_names) {
      
      output_table[['score']] <- output_table[['score']] + output_table[[var]]
      
    }
    
    if(tibble_output) {
      
      if(show_var_cuts) {
        
        return(output_table)
        
      } else {
        
        return(output_table[, !names(output_table) %in% new_variable_names])
        
      }
      
    } else {
      
      return(output_table$score)
      
    }
  
  }
  
  calculate_var_sum <- function(inp_table, var_vector, tibble_output = F) {
    
    ## calculates a simple sum of given variables
    
    output_table <- inp_table
    
    output_table[['score']] <- 0
    
    for(var in var_vector) {
      
      output_table[['score']] <- output_table[['score']] + output_table[[var]]
      
    }
    
    if(tibble_output){
      
      return(output_table)
      
    } else {
      
      return(output_table[['score']])
      
    }
    
  }


