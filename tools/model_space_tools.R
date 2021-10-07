# This toolbox script contains general tools for model space generation

# libraries ----

  library(plyr)
  library(dplyr)
  library(purrr)
  library(gtools)
  library(furrr)

# formula creation by sting concatenation and the 'as.formula' function -----

  make_formula <- function(response, variable_vector, ...) {
    
    ## generates a formulas of the following type:
    ## response ~ variable1 + variable2 + ... + variable_n
    
    variable_vector <- paste0(variable_vector, collapse = '+')
    
    return(paste(response, '~', variable_vector) %>% 
             as.formula(object = ., ...))
    
  }

# combination tree generating function -----

  make_comb_tree <- function(inp_vector, min_length = 1, max_length = length(inp_vector)) {
    
    ## for a given inp_vector a list of all n-long combinations is returned, where min_length =< n =< max_length
    
    comb_list <- list()
    
    for(i in min_length: max_length) {
      
      res_list <- combinations(length(inp_vector), i, inp_vector) %>% 
        asplit(., 1)
      
      comb_list <- c(comb_list, res_list)
      
    } 
    
    return(comb_list)
    
  }

# error-resistant model generator -----

  make_model <- function(mod_fun, formula, inp_data, ...) {
    
    ## returns a model defined by a given modeling function, formula, input data and optional arguments
    ## or error message
    
    out_model <- try(mod_fun(formula = formula, data = inp_data, ...), silent = T)
    
    if(any(class(out_model) == 'try-error')) {
      
      return(out_model[1])
      
    } else {
      
      return(out_model)
    }
    
  }

# serial modeling function ----

  make_model_list <- function(mod_fun, 
                              var_comb_list, 
                              response, 
                              inp_data, 
                              parallel = F, 
                              globals = T, ...) {
    
    ## for a given list of variable combinations, a list of models created by mod_fun function 
    ## with the given response is returned. ... are additional arguments passed to the modeling function
    ## can be parallelized. The globals argument enables export of objects to the workers
    
    message(paste('Models to calculate:', length(var_comb_list)))
    start_time <- Sys.time()
    
    if(parallel) {

      require(furrr)

      plan('multiprocess')
      
      model_list <- future_map(var_comb_list, function(x) make_model(mod_fun = mod_fun, 
                                                                     formula = make_formula(response, x), 
                                                                     inp_data = inp_data, ...), 
                               .options = future_options(globals = globals, lazy = F))
      

    } else {
      
      model_list <- map(var_comb_list, function(x) make_model(mod_fun = mod_fun, 
                                                              formula = make_formula(response, x), 
                                                              inp_data = inp_data, ...))
      
    }
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    return(model_list)
    
  }

# model space generating function ----

  gen_model_space <- function(mod_fun, 
                              variable_vect, 
                              response, 
                              inp_data, 
                              parallel = T, 
                              min_var = 1, 
                              max_var = length(variable_vect), 
                              globals = T,  ...) {
    
    ## a wrapper for the make_model_list and make_comb_tree functions. Generates a model space for a given reponse, 
    ## variable list (from which combinations are generated) and modeling function. min_param and max_param enable
    ## specifying the minimal and maximal number of dependent variables in the model, respectively.
    ## Can be run in parallel. ... are additional arguments for the modeling function.

    message('Model space generator')
    start_time <- Sys.time()
    
    message(paste('Generating combination tree for n =', min_var, '..', max_var, 'fold combinations'))
    
    var_comb_list <- make_comb_tree(inp_vector = variable_vect, 
                                    min_length = min_var, 
                                    max_length = max_var)
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    model_list <- make_model_list(mod_fun = mod_fun, 
                                  var_comb_list = var_comb_list, 
                                  response = response, 
                                  inp_data = inp_data, 
                                  parallel = parallel, 
                                  globals = globals)
    
    return(list(vars = var_comb_list, 
                models = model_list))
    
  }
  
# LRT model space scanning tools ----
  
  test_lrt <- function(inp_model, null_model_formula = NULL, inp_data = NULL, 
                       mod_fun = NULL, test = 'Chi') {
    
    ## an error-insensitive LRT caclulating function for moset fixed effect model classes
    ## test specifies which test ha to be used to compare deviance differences (see: anova)
    ## memory saver reduces model frame to a reponse-only table and will work only with 'real' null models
    
    null_model <- make_null_model(inp_model = inp_model, 
                                  null_model_formula = null_model_formula, 
                                  inp_data = inp_data, 
                                  mod_fun = mod_fun, ...)
    
    res_table <- anova(inp_model, null_model, test = test)
    
    return(res_table[2, 2:ncol(res_table)])
    
  }
  
  make_null_model <- function(inp_model, null_model_formula = NULL, inp_data = NULL, 
                              mod_fun = NULL, simple_return = T...) {
    
    ## generates a null model for a given input model. The memory saver function is likely
    ## to work with 'real' null models

    null_model <- try(mod_fun(formula = null_model_formula, 
                              data = model.frame(inp_model, data = inp_data, ...)), 
                      silent = T)

    return(null_model)
    
  }
  
  scan_mod_space_lrt <- function(inp_model_list, null_model_formula = NULL, 
                                 mod_fun = NULL, inp_data = NULL, test = 'Chi', 
                                 extr_fun = NULL, parallel = c(T, T, T), tibble_output = T, 
                                 globals = T, simple_lrt = T, ...) {
    
    ## a list wrapper for scanning a space of models (inp_model_list) with a series of LRT tests
    ## ... are additional arguments passed to mod_fun modeling function
    ## extr_fun defines a function extracting additional data from the model and concatenating them with the result
    ## The globals argument enables export to the workers
    ## the simple_lrt option will work with Cox models (it's fast), where LRT results are extracted from
    ## the model summary directly
    
    ## generating null model list
    
    start_time <- Sys.time()
    
    if(!simple_lrt) {
      
      message(paste('Generating null model list for n =', length(inp_model_list), 'input models'))
      
      if(parallel[1]) {
        
        plan('multiprocess')
        
        null_model_list <- inp_model_list %>% 
          future_map(function(x) make_null_model(inp_model = x, 
                                                 null_model_formula = null_model_formula, 
                                                 inp_data = inp_data, 
                                                 mod_fun = mod_fun, ...), 
                     .options = future_options(globals = globals, lazy = F))
        
      } else {
        
        null_model_list <- inp_model_list %>% 
          map(function(x) make_null_model(inp_model = x, 
                                          null_model_formula = null_model_formula, 
                                          inp_data = inp_data, 
                                          mod_fun = mod_fun, ...))
        
      }
      
      message(paste('Time elapsed:', Sys.time() - start_time))
      
    }
    
    ## generating list of anova results
    
    start_time <- Sys.time()
    
    message(paste('LRT testing n =', length(inp_model_list), 'input models'))
    
    if(!simple_lrt) {
      
      if(parallel[2]) {
        
        plan('multiprocess')
        
        lrt_results <- future_map2(inp_model_list, 
                                   null_model_list, 
                                   function(x, y) try(anova(x, y, test = test), silent = T))
        
      } else {
        
        lrt_results <- map2(inp_model_list, 
                            null_model_list, 
                            function(x, y) try(anova(x, y, test = test), silent = T))
        
      }
      
      lrt_results <- lrt_results %>% 
        map(function(x) try(x[2, 2:ncol(x)], silent = T))
      
    } else {
      
      lrt_results <- inp_model_list %>% 
        map(summary) %>% 
        map(function(x) x$logtest) %>% 
        map(set_names, 
            c('ChiSq', 'Df', 'P(>|Chi|)'))
      
    }
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    
    if(!is.null(extr_fun)) {
      
      if(parallel[3]) {
        
        message('Extracting additional model information')
        start_time <- Sys.time()
        
        plan('multiprocess')
        
        extr_results <- future_map(inp_model_list, function(x) extr_fun(x), 
                                   .options = future_options(globals = globals, lazy = F))
        
        message(paste('Time elapsed:', Sys.time() - start_time))
        
        
      } else {
        
        
        message('Extracting additional model information')
        start_time <- Sys.time()
        
        extr_results <- map(inp_model_list, function(x) extr_fun(x))
        
        message(paste('Time elapsed:', Sys.time() - start_time))
        
        
      }
      
      
    }
    
    if(tibble_output){
      
      lrt_results <- lrt_results %>% 
        do.call('rbind', .) %>% 
        as_tibble %>% 
        mutate(model_ID = names(inp_model_list))
      
      if(!is.null(extr_fun)) {
        
        extr_results <- extr_results %>% 
          do.call('rbind', .) %>% 
          as_tibble()
        
        return(cbind(lrt_results, extr_results) %>%  as_tibble())
        
      } else {
        
        return(lrt_results)
        
      }
      
    } else {
      
      if(!is.null(extr_fun)) {
        
        return(list(lrt_results = lrt_results, extr_results = extr_results))
        
      } else {
        
        return(lrt_results)
        
      }
      
    }
    
  }
  
  large_mod_space_lrt <- function(inp_model_list, null_model_formula = NULL, 
                                  mod_fun = NULL, inp_data = NULL, test = 'Chi', 
                                  extr_fun = NULL, globals = T, chunk_size = 1000, 
                                  parallel = c(T, F, F), main_parallel = F, ...) {
    
    ## tackles large model lists in a seria-parallel manner, chunk size controls
    ## the size of input list chunk and hence number of internal iterations
    
    start_time = Sys.time()
    
    message(paste(rep('=', 60), collapse = ''))
    message('Scanning LRT and extracting model data for total of n = ', length(inp_model_list), ' models')
    
    model_chunks <- split_vec(inp_vector = inp_model_list, chunk_size = chunk_size)
    
    if(!main_parallel) {
      
      lrt_result <- list()
      
      for(i in model_chunks) {
        
        lrt_result <- scan_mod_space_lrt(inp_model_list = i, 
                                         null_model_formula = null_model_formula, 
                                         mod_fun = mod_fun, 
                                         inp_data = inp_data, 
                                         test = test, 
                                         extr_fun = extr_fun, 
                                         parallel = parallel, 
                                         tibble_output = T, 
                                         globals = globals, ...) %>% 
          c(lrt_result, .)
        
      }
      
    } else {
      
      plan('multiprocess')
      
      lrt_result <- model_chunks %>% 
        future_map(function(x) scan_mod_space_lrt(inp_model_list = x, 
                                                  null_model_formula = null_model_formula, 
                                                  mod_fun = mod_fun, 
                                                  inp_data = inp_data, 
                                                  test = test, 
                                                  extr_fun = extr_fun, 
                                                  parallel = c(F, F, F), 
                                                  tibble_output = T, ...), 
                   .options = future_options(globals = globals))
    
    }
    
    lrt_result <- lrt_result %>% 
      do.call('rbind', .)
    
    message('Time elapsed ', Sys.time() - start_time)
    message(paste(rep('=', 60), collapse = ''))
    
    return(lrt_result)
    
  }

# C-index calculation ----
  
  calculate_c_space <- function(inp_model_list, result_table = NULL, parallel = T) {
    
    ## calculates concordance indexes for a list of models (for details, see: cocnordance function in survival)
    ## if an output table wiht preexisting results and a model_ID column is provided, C-indexes are merged into it
    
    require(survival)
    
    calc_concordance <- function(inp_model) {
      
      ## error-resistant concordance calculating function
      
      conco_output <- try(concordance(inp_model)$concordance, silent = T)
      
      if(any(class(conco_output) == 'try-error')) {
        
        return('calculation error')
        
      } else {
        
        return(conco_output)
        
      }
      
    }
    
    if(parallel) {
      
      require(doParallel)
      require(furrr)
      
      cl <- makeCluster(detectCores() - 1)
      plan(future::cluster, workers = cl, homogeneous = F)
      
      message(paste('C-index calculation. Models to test:', length(inp_model_list)))
      start_time <- Sys.time()
      
     conco_table <- inp_model_list %>% 
        future_map(calc_concordance) %>% 
        do.call('rbind', .) %>% 
        as_tibble %>% 
        set_names(c('c_index')) %>% 
        mutate(model_ID = names(inp_model_list))
     
     stopCluster(cl)
      
    } else {
      
      message(paste('C-index calculation. Models to test:', length(inp_model_list)))
      start_time <- Sys.time()
      
      conco_table <- inp_model_list %>% 
        map(calc_concordance) %>% 
        do.call('rbind', .) %>% 
        as_tibble %>% 
        set_names(c('c_index')) %>% 
        mutate(model_ID = names(inp_model_list))
      
    }
    
    if(!is.null(result_table)) {
      
      message('Updating the result table with C-indexes')
      start_time <- Sys.time()
      
      conco_table <- left_join(result_table, conco_table, by = 'model_ID')
      
      message(paste('Time elapsed:', Sys.time() - start_time))
      
    } 
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    return(conco_table)
    
  }
  
# Covariance penalty calculation ----
  
  calculate_cp_space <- function(inp_model_list, result_table = NULL, pen_fun = AIC, parallel = T, ...) {
    
    ## calculates covariance penalty for a list of models (for details, see: AIC())
    ## if an output table wiht preexisting results and a model_ID column is provided, Cp values are merged into it
    
    calc_cp <- function(inp_model, pen_fun, ...) {
      
      ## error-resistant Cp calculating function
      
      cp_output <- try(pen_fun(inp_model, ...), silent = T)
      
      if(any(class(cp_output) == 'try-error')) {
        
        return('calculation error')
        
      } else {
        
        return(cp_output)
        
      }
      
    }
    
    if(parallel) {
      
      require(doParallel)
      require(furrr)
      
      cl <- makeCluster(detectCores() - 1)
      plan(future::cluster, workers = cl, homogeneous = F)
      
      message(paste('Cp calculation. Models to test:', length(inp_model_list)))
      start_time <- Sys.time()
      
      cp_table <- inp_model_list %>% 
        future_map(calc_cp, pen_fun = pen_fun, ...) %>% 
        do.call('rbind', .) %>% 
        as_tibble %>% 
        set_names(c('Cp')) %>% 
        mutate(model_ID = names(inp_model_list))
      
      stopCluster(cl)
      
    } else {
      
      message(paste('Cp calculation. Models to test:', length(inp_model_list)))
      start_time <- Sys.time()
      
      cp_table <- inp_model_list %>% 
        map(calc_cp, pen_fun = pen_fun, ...) %>% 
        do.call('rbind', .) %>% 
        as_tibble %>% 
        set_names(c('Cp')) %>% 
        mutate(model_ID = names(inp_model_list))
      
    }
    
    if(!is.null(result_table)) {
      
      message('Updating the result table with Cp-stats')
      start_time <- Sys.time()
      
      cp_table <- left_join(result_table, cp_table, by = 'model_ID')
      
    } 
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    return(cp_table)
    
  }
  
# A result cleaning wrapper for stage objects, i.e. lists containers for vars, models and results ----
  
  clean_results <- function(stage_object) {
    
    ## renames the result table, adjusts the LRT pvalues (Benjamini-Hochberg)
    ## adds variable names and checks significance of relevant (non-intercept)
    ## model estimates
    
    start_time = Sys.time()
    
    message('Cleaning results')
    
    result_table <- stage_object$results
    
    if(('P(>|Chi|)' %in% names(result_table)) | ('test' %in% names(result_table))) {
      
      if(!all(c('r_cod', 
             'r_mer', 
             'r_mev') %in% names(result_table))) {
        
        result_table <- result_table %>% 
          set_names(c('chisg', 'df', 'p_lrt', 'model_id', 'p_est'))
        
      } else {
        
        result_table <- result_table %>% 
          set_names(c('chisg', 'df', 'p_lrt', 'model_id', 'p_est', 'rsq_cod', 'rsq_mer', 'rsq_mev'))
        
      }
      
      
    } else if('Pr(>F)' %in% names(result_table)) {
      
      result_table <- result_table %>% 
        set_names(c('rss', 'df', 'ssq', 'f', 'p_lrt', 'model_id', 'p_est'))
      
    } else {
      
      message('The stage object has no valis result table')
      return('Error: no valid result table detected')
      
    }
    
    ## BH correction, extraction of relavant estimate p values
    
    result_table <- result_table %>% 
      mutate(p_lrt_adj = p.adjust(p_lrt, method = 'BH')) %>% ## Benjamini-Hochberg correction
      mutate(p_est_relevant = map(result_table$p_est, function(x) x[names(x) != '(Intercept)'])) ## relevant estimates' p
    
    ## getting variable number
    
    len_table <- stage_object$var %>% 
      map(function(x) length(x)) %>% 
      reduce(c) %>% 
      tibble(model_id = names(stage_object$vars), 
             var_number = .)
    
    ## getting the AIC and C-index data
    
    message(paste('Calculating AIC for n = ', length(stage_object$vars), 'models'))
    
    aic_table <- calculate_cp_space(stage_object$models, parallel = F) %>% 
      set_names(c('aic', 'model_id'))
    
    message(paste('Calculating C-indexes for n = ', length(stage_object$vars), 'models'))
    
    c_table <- calculate_c_space(stage_object$models, parallel = F) %>% 
      set_names(c('c_index', 'model_id'))
    
    # updating the result table
    
    result_table <- add_var_names(result_table, stage_object$vars)
    
    result_table <- list(result_table, len_table, aic_table, c_table) %>% 
      reduce(function(x, y) left_join(x, y, by = 'model_id'))
    
    ## identifying significant models: p_lrt_adj < 0.05 and all estimate p values < 0.05
    
    result_table <- result_table %>% 
      mutate(significant = map(result_table$p_est_relevant, function(x) all(x < 0.05)) %>% 
               reduce(c)) %>% 
      mutate(significant = significant & p_lrt_adj < 0.05)
    
    
    return(list(vars = stage_object$vars, 
                models = stage_object$models,
                results = result_table))
    
  }
  
# Other result cleaning functions ------
  
  ## a function identifying significant var combinations and models
  
    id_signif <- function(stage_object) {
      
      result_table <- stage_object$results %>% 
        filter(significant)
      
      signif_ids <- result_table$model_id
      
      return(list(vars = stage_object$vars[signif_ids], 
                  models = stage_object$models[signif_ids], 
                  results = result_table))
      
    }
  
  ## a function for indetifying and removing nested variable sets
  
    filter_out_nested <- function(stage_object) {
      
      ## filters out variable sets, models and summaries for models
      ## which are nested in other var sets/models
      
      nested_var_vec <- find_nested_var_sets(stage_object$vars)
      
      non_nested_ids <- names(stage_object$vars[!nested_var_vec])
      
      new_stage_object <- list()
      
      new_stage_object$vars <- stage_object$vars[non_nested_ids]
      new_stage_object$models <- stage_object$models[non_nested_ids]
      new_stage_object$results <- stage_object$results %>% 
        filter(model_id %in% non_nested_ids)
      
      return(new_stage_object)
      
    }
  
  ## a function identifying top n best-performing models
  
    identify_top <- function(stage_object, n_best = 10, sel_variable = 'c_index', new_var_name = 'best_c') {
      
      ## indentifies top n best-predicting models: adds a variable to the result table
      
      new_result_tbl <- stage_object$results %>% 
        mutate(sel_variable = .[[sel_variable]]) %>% 
        top_n(n_best, sel_variable) %>% 
        mutate(best_model = TRUE) %>% 
        select(model_id, best_model)
      
      new_result_tbl <- left_join(stage_object$results, new_result_tbl, by = 'model_id') %>% 
        mutate(best_model = ifelse(is.na(best_model), FALSE, best_model))
      
      new_result_tbl[[new_var_name]] <- new_result_tbl$best_model
      new_result_tbl$best_model <- NULL
      
      return(list(vars = stage_object$vars, 
                  models = stage_object$models, 
                  results = new_result_tbl))
      
    }
  
  ## identifying best class model (i.e. top AIC and top C-index)
  
    identify_best_class <- function(stage_object, sel_vars = c('best_c', 'best_aic')) {
    
    new_result_tbl <- stage_object$results %>% 
      mutate(best_class = .[[sel_vars[1]]] & .[[sel_vars[2]]])
    
    return(list(vars = stage_object$vars, 
                models = stage_object$models, 
                results = new_result_tbl))
    
  }
  
  
# Result saving function ----

  result_editor <- function(stage_object, file_suffix) {
      
      ## changes the format of the given result tables by adding p_values for the relevant estimates
      ## and saves in on a disc with the provided file suffix
      
    stage_object %>% map(function(x) x[['results']]) %>% 
        map(function(x) select(x, - p_est)) %>% 
        map(function(x) mutate(x, p_est_relevant = map_chr(x$p_est_relevant, zip_vector))) %>% 
        walk2(., paste(names(.), file_suffix, sep = ''), write_tsv)
      
  }    
   
# novel parameter number and fraction for a paramater set -----
  
  determine_novel <- function(inp_var_list, old_vars, result_table = NULL) {
    
    ## for each element of the given variable set list
    ## number of novel paramaters is determined
    ## old paramater names are provided by the old_vars vector
    ## updates a result table if provided
    
    message('Determining the number and fraction of novel paramater per parameter set')
    message(paste('Sets to calculate:', length(inp_var_list)))
    
    start_time = Sys.time()
    
    novelty_res <- map(inp_var_list, function(x) c(no_novel_var = sum(!x %in% old_vars), 
                                                   frac_novel_var = sum(!x %in% old_vars)/length(x))) %>% 
      do.call('rbind', .) %>% 
      as_tibble %>% 
      mutate(model_ID = names(inp_var_list))
    
    if(!is.null(result_table)) {
      
      novelty_res <- left_join(result_table, novelty_res, by = 'model_ID')
      
    }
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    return(novelty_res)
    
  }
  
# testing for nesting -----
  
  test_nested <- function(inp_vector, other_vector) {
    
    # checks whether a vector is nested within the other given vector
    # i.e. whether all elements of inp_vector are present in the other_vector
    
    return(all(inp_vector %in% other_vector))
    
  }
  
  identify_nested <- function(inp_vector, other_vector_list) {
    
    ## checks whteher the given vector is nested in at least one 
    ## vectors in the vector list
    
    for(i in other_vector_list) {
      
      if(test_nested(inp_vector, i)) {
        
        return(TRUE)
        
      }
      
    }
    
    return(FALSE)
    
  }
  
  find_nested_var_sets <- function(var_set_list) {
    
    ## In a named list containing multiple variable sets, the sets nested within other variables sets in the list
    ## are identified
    
    nested_sets <- c()
    
    for(i in names(var_set_list)) {
      
      nested_sets <- c(nested_sets, identify_nested(var_set_list[[i]], var_set_list[names(var_set_list) != i]))
      
    }
    
    return(set_names(nested_sets, names(var_set_list)))
    
  }
  
# result table annotation -----
  
  add_var_names <- function(inp_table, var_list) {
    
    # adds a column with names of variables used for modeling.
    # The annotation needs a model_ID column in the inp_table and a list named with model_IDs
    
    message(paste('Adding variable names. Names to add:', length(var_list)))
    start_time = Sys.time()
    
    annotated_table <- var_list %>%
      map(function(x) paste(x, collapse = ', ')) %>% 
      do.call('rbind', .) %>% 
      as_tibble %>% 
      set_names(c('variables')) %>% 
      mutate(model_id = names(var_list)) %>% 
      left_join(inp_table, ., by = 'model_id')
    
    message(paste('Time elapsed:', Sys.time() - start_time))
    
    return(annotated_table)
    
  }
  
# varia ----
  
  fast_rbind <- function(inp_list, retain_class = F) {
    
    new_colnames <- names(inp_list[[1]])
    
    output_tibble <- matrix(unlist(inp_list), nrow = length(inp_list), byrow = T) %>% 
      as_tibble %>% 
      set_names(new_colnames)
    
    if(retain_class){
      
      output_tibble <- map_dfc(output_tibble, parse_guess)
      
    }
    
    return(output_tibble)
    
  }
  
  split_vec <- function(inp_vector, chunk_size) {
    
    return(split(inp_vector, ceiling(seq_along(inp_vector)/chunk_size)))
    
  }
  
  zip_vector <- function(inp_vector) {
    
    ## a name zipping function
    
    out_vector <- map2_chr(names(inp_vector), 
                           inp_vector, 
                           function(x, y) paste('p', x, ' = ', signif(y, 2), sep = '')) %>% 
      paste(., collapse = ', ')
    
    return(out_vector)
    
  }

# END ----