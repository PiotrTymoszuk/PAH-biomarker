# libraries ----

  require(tidyverse)
  require(cowplot)

# varia -----
    
    my_wrap <- function(inp_vector, segment_length) {
      
      ## divides a vector into segments of a given length
      ## elements within the segment are collapsed with comma, segments collapsed with a new line
      
      split_vec <- split(inp_vector, ceiling(seq_along(inp_vector)/segment_length))
      
      wrap_vec <- split_vec %>% 
        map(function(x) paste(x, collapse = ', ')) %>% 
        paste(collapse = '\n')
      
      return(wrap_vec)
      
    }
    
    
# Figure object and saving functions -----
    
    add_title <- function(figure_plot, title_text, 
                          x = 0.05, y = 0.5, size = 12, hjust = 0, vjust = 0.5, rel_heights = c(0.05, 0.95), 
                          face = 'bold', ...) {
      
      new_figure <- plot_grid(ggdraw() + 
                                draw_text(title_text, 
                                          x = x, 
                                          y = y, 
                                          size = size, 
                                          hjust = hjust, 
                                          vjust = vjust, 
                                          fontface = face, 
                                          ...), 
                              figure_plot, 
                              nrow = 2, 
                              rel_heights = rel_heights)
      
      return(new_figure)
      
    }
    
    save_figure <- function(figure_plot, figure_name, h, w = 180, target_folder = 'figures', ...) {
      
      ## an error-resistant figure-saving function
      
      insert_msg(paste('Saving:', figure_name))
      
      enter_directory(target_folder)
      
      tryCatch(ggsave(filename = figure_name, 
                      plot = figure_plot, 
                      width = w, 
                      height = h, 
                      units = 'mm', 
                      ...), 
               error = simpleError('Saving failed'), 
               finally = go_proj_directory())
      
    }
    
    as_figure_object <- function(figure_plot, figure_label, h, w = 180) {
      
      ## creates a ready-to-save list with the figure and metadata
      
      figure_obj <- list(plot = figure_plot, 
                         figure_label = figure_label, 
                         h = h, 
                         w = w)
      
      return(figure_obj)
      
    }
    
    save_figure_object <- function(figure_obj, format = 'pdf', target_folder = 'figures', ...) {
      
      ## saves a figure object
      
      save_figure(figure_plot = figure_obj$plot, 
                  figure_name = paste(figure_obj$figure_label, format, sep = '.'), 
                  h = figure_obj$h, 
                  w = figure_obj$w, 
                  target_folder = target_folder, ...)
      
    }
    
    