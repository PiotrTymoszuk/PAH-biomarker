# Generates paper tables and figures

# tools --------

  library(knitr)
  library(rmarkdown)
  library(bookdown)
  library(flextable)
  library(writexl)
  
  library(soucer)
  library(figur)
  
  source_all(c('./tools/project_tools.R'), 
             message = TRUE, crash = TRUE)
  
  insert_head()
  
# paper scripts -----
  
  insert_msg('Sourcing the paper scripts')
  
  source_all(c('./paper scripts/tables.R', 
               './paper scripts/figures.R',
               './paper scripts/biblio.R', 
               './paper scripts/links.R', 
               './paper scripts/render.R'), 
             message = TRUE, crash = TRUE)
  
# END -----
  
  insert_tail()