# This program executes the scripts of the PAH biomarker project

# libraries -----

  library(soucer)

# listing figure scripts ----

  print(source_all(c('import.R', 
                     'exploration.R', 
                     'analysis.R', 
                     'paper.R'), 
                   message = TRUE, crash = FALSE))

# END -----