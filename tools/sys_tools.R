# This script contains the toolbox for directory navigation, reading and importing flow-jo tables generated in course
# of the low-density granulocyte RA study

# libraries -----

  library(plyr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringi)
  library(scriptName)

# file system toolbox ----

  enter_directory <- function(path, relative = T) {
  
  ## creates a directory and sets the working directory respectively
  ## if relative is True, a subfolder of the current working directory is entered
  
  if (relative) {
    
    curr_path <- getwd()
    
    dir.create(paste(curr_path, '/', path, sep = ''), showWarnings = F)
    
  } else {
    
    dir.create(path, showWarnings = F)
    
  }
  
  
  setwd(path)
  
}

  go_proj_directory <- function() {
  
  ## sets the working directory back to the project directory
  
  setwd(rprojroot::find_rstudio_root_file())
  
  }
  
  go_up_directory <- function(){
    
    ## sets the working directory to the direct upward directory
    
    require(stringi)
    
    splitted_dir <- stri_split_fixed(getwd(), pattern = '/') %>% 
      unlist
    
    up_dir <- splitted_dir[1:length(splitted_dir) - 1] %>% 
      reduce(function(x, y) paste(x, '/', y, sep = ''))
    
    setwd(up_dir)
    
  }
  
  list_dirs <- function(){
    
    ## lists all the subfolders of the given directory in a non-recursive manner
    
    return(list.dirs(full.names = F, recursive = F))
    
  }
  
  list_file_type <- function(directory = '.', extension = '.txt') {
    
    ## looks for potential tables (e.g. = .txt) files in a given directory
    
    return(list.files(path = directory, pattern = paste('\\', extension, '$', sep = '')))
    
  }

# messages ------
  
  insert_msg <- function(text = NULL, separator = '>', sep_len = 60) {
    
    message(rep(separator, sep_len))
    
    if(!is.null(text)) {
      
      message(text)
      
    }
    
  }
  
  insert_head <- function(def_prefix = 'Executing:') {
    
    insert_msg(paste(def_prefix, current_filename()), separator = '<>')
    
  }
  
  insert_tail <- function(def_suffix = 'succesfully sourced') {
    
    insert_msg(paste(current_filename(), def_suffix), separator = '<>')
    
  }