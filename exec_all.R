# This program executes the scripts of the PAH biomarker project

# libraries -----

  require(plyr)
  require(tidyverse)

# listing figure scripts ----

  executable_scripts <- c('data_import.R', 
                          'univariate_modeling.R', 
                          'multivariate_modeling.R', 
                          'participant_clustering.R', 
                          'render_paper.R')
  
# error container ----
  
  exec_errors <- list()
  
# executing the script list ----
  
  for(file in executable_scripts) {
    
    diagn_output <- try(source(file, encoding = 'UTF-8'), silent = T) ## UTF-8 encoding cause of greek letters 
                                                                      ##in plots and tables
    
    if(class(diagn_output) == 'try-error') {
      
      exec_errors <- c(exec_errors, 'execution failed')
      
    } else {
      
      exec_errors <- c(exec_errors, 'execution successful')
      
    }

  }
  
  exec_errors <- exec_errors %>% 
    reduce(c) %>% 
    tibble(execution_result = .)
  
  exec_errors <- exec_errors %>% 
    mutate(file = executable_scripts) %>% 
    select(file, execution_result)
  
  rm(diagn_output, executable_scripts)
  
  message(rep('-', 60))
  print(exec_errors)
  message(rep('-', 60))
  
# saving execution results on the disc ----

  exec_errors %>% 
    write_tsv('exec_log.log')
  
  save.image()
  
# END ----