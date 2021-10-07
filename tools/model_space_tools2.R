# This toolbox script extends contains general tools for model space generation with model space generators
# enabling instant model filtering

# libraries ----

  library(plyr)
  library(dplyr)
  library(purrr)
  library(gtools)
  library(furrr)

  source('./tools/model_space_tools.R')

# model filtering wrapper -----
  
  filter_model <- function(inp_model, filter_fun) {
    
    ## with a provided filtering function (returning TRUE or FALSE), model is kept or rejected
    
    if(filter_fun(inp_model)) {
      
      return(inp_model)
      
    } else {
      
      return(NULL)
      
    }
    
  }

# serial modeling function ----
  
  mdl_lst_ft_ <- function(mod_fun, 
                          var_comb_list, 
                          response, 
                          inp_data, 
                          filter_fun = NULL, 
                          parallel = T, 
                          globals = globals, ...) {
    
    ## for a given list of variable combinations, a list of models created by mod_fun function 
    ## with the given response is returned. ... are additional arguments passed to the modeling function
    ## filter_fun enables instant filtering out of unwanted models
    
    message(paste('Models to calculate:', length(var_comb_list)))
    start_time <- Sys.time()
    
    if(parallel) {
      
      plan('multiprocess')
      
      model_list <- var_comb_list %>% 
        future_map(function(x) make_model(mod_fun = mod_fun, 
                                          formula = make_formula(response, x), 
                                          inp_data = inp_data, ...), 
                   .options = future_options(globals = globals))
      
      plan('sequential')
      
    } else {

      model_list <- var_comb_list %>% 
        map(function(x) make_model(mod_fun = mod_fun, 
                                   formula = make_formula(response, x), 
                                   inp_data = inp_data, ...))
      
    }
   
    
    if(is.null(filter_fun)) {
      
      return(list(vars = var_comb_list, 
                  models = model_list))
      
    }
    
    message('Filtering model list')
    
    model_list <- model_list %>% 
      map(filter_model, filter_fun = filter_fun) %>% 
      compact
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    return(list(vars = var_comb_list[names(model_list)], 
                models = model_list))
    
  }

  make_model_list_ft <- function(mod_fun, 
                                 var_comb_list, 
                                 response, 
                                 inp_data, 
                                 filter_fun, 
                                 chunks = 1, 
                                 parallel = T, 
                                 globals = globals, ...) {
    
    ## for a given list of variable combinations, a list of models created by mod_fun function 
    ## with the given response is returned. ... are additional arguments passed to the modeling function
    ## filter_fun enables instant filtering out of unwanted models. Chunks define number of serial
    ## iterations run on model list chunks of the equal size
    
    sep = rep('=', 60)
    master_time = Sys.time()
    chunk_no = 1
    
    message(sep)
    message('Generating a model space for a total of ', length(var_comb_list), 'models')
    message(sep)
    
    if(chunks == 1) {
      
      ## by default, a 100% parallel run
      
      return(mdl_lst_ft_(mod_fun = mod_fun, 
                         var_comb_list = var_comb_list, 
                         response = response, 
                         inp_data = inp_data, 
                         filter_fun = filter_fun, 
                         parallel = parallel, 
                         globals = globals, ...))
      
      
    } else {
      
      chunk_size = ceiling(length(var_comb_list)/chunks)
      
      var_list_split <- split_vec(var_comb_list, chunk_size = chunk_size)
      
      result_vars <- list()
      result_models <- list()
      
      for(chunk in var_list_split) {
        
        message('Working on list chunk no.', chunk_no, ' of ', chunks)
        chunk_no = chunk_no + 1
        
        curr_output <- mdl_lst_ft_(mod_fun = mod_fun, 
                                   var_comb_list = chunk, 
                                   response = response, 
                                   inp_data = inp_data, 
                                   filter_fun = filter_fun, 
                                   globals = globals, ...)
        
        result_vars <- c(result_vars, curr_output$vars)
        result_models <- c(result_models, curr_output$models)
        
      }
      
      message(sep)
      message('Total time elapsed:', Sys.time() - master_time)
      message(sep)
      
      return(list(vars = result_vars, 
                  models = result_models))
      
    }

  }
  
# model space generating function ----

  gen_model_space_ft <- function(mod_fun, 
                              variable_vect, 
                              response, 
                              inp_data, 
                              filter_fun, 
                              min_var = 1, 
                              max_var = length(variable_vect), 
                              add_variables = NULL,  
                              globals = T, ...) {
    
    ## a wrapper for the make_model_list and make_comb_tree functions. Generates a model space for a given reponse, 
    ## variable list (from which combinations are generated) and modeling function. min_param and max_param enable
    ## specifying the minimal and maximal number of dependent variables in the model, respectively.
    ##  ... are additional arguments for the modeling function.
    ## The add_variables vector enables inclusion of some non-changing confoundes as age or sex
    
    message('Model space generator')
    start_time <- Sys.time()
    
    message(paste('Generating combination tree for n =', min_var, '..', max_var, 'fold combinations'))
    
    var_comb_list <- make_comb_tree(inp_vector = variable_vect, 
                                    min_length = min_var, 
                                    max_length = max_var) %>% 
      set_names(paste('test_', 1:length(.), sep = ''))
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    if(!is.null(add_variables)) {
      
      var_comb_list <- var_comb_list %>% 
        map(function(x) c(x, add_variables))
      
    }
    
    return(make_model_list_ft(mod_fun = mod_fun, 
                              var_comb_list = var_comb_list, 
                              response = response, 
                              inp_data = inp_data, 
                              filter_fun = filter_fun, 
                              chunks = 1, 
                              globals = globals))

  }
  


# END ----