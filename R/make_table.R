#' @title Helper function: Make a table
#' 
#' @description Make table for simulation results.
#'
#'
#' @param df data frame
#' @param x character specifying the variable of interest
#' @param rnames vector containing row names
#' @param cnames vector containing column names
#' @param panel_names vector containing the panel names
#' @param title title of the table
#' @param filename write table to the filename
#'
#' @return saves a LaTeX table under the filename
#' 
#' 


make_table <- function(df, x, rnames, cnames = c("$H = 1$" , "$H = 10 $",  "$H = 25$"), panel_names = NULL, title, filename){
  
  df_tab <- df %>% 
    pivot_wider(names_from = N, 
                values_from = x, 
                names_prefix = paste0(x, "_N")) %>% 
    pivot_wider(names_from = H, 
                values_from = c(paste0(x, "_N50"), paste0(x, "_N200"), paste0(x, "_N800"))) %>% 
    dplyr::select(starts_with(x))
  
  
  # Into data frame
  my_tab <- as.data.frame(df_tab)
  
  # Number of rows in each panel
  no_row <- length(rnames)
  
  # Add rownames
  my_tab <- cbind(stat = rep(rnames, nrow(my_tab)/no_row), my_tab)
  
  # If panel_names is given
  if(is.null(panel_names) == FALSE){
    
    # Index of the panel
    panel_index <- rep(no_row, length(panel_names))
    names(panel_index) <- panel_names
    
    # Export to LaTeX
    kable(my_tab, col.names = c("", rep(cnames, 3)), booktabs = TRUE, format = "latex", escape = FALSE, digits=3, 
          caption = title) %>%
      add_header_above(c(" " = 1, "$T = 50$" = 3, "$T = 200$" = 3,  "$T = 800$" = 3), escape = FALSE) %>%
      pack_rows(
        index = panel_index,
        latex_align = "c",
        hline_before = TRUE,
        hline_after = TRUE,
        escape = FALSE
      ) %>%
      column_spec (1:10, latex_column_spec = "l|ccc|ccc|ccc") %>%
      kable_styling(latex_options="scale_down") %>% 
      save_kable(filename)
  }
  
  # No panel_names
  if(is.null(panel_names)){
    
    # Export to LaTeX
    kable(my_tab, col.names = c("", rep(cnames, 3)), booktabs = TRUE, format = "latex", escape = FALSE, digits=3, 
          caption = title) %>%
      add_header_above(c(" " = 1, "$T = 50$" = 3, "$T = 200$" = 3,  "$T = 800$" = 3), escape = FALSE) %>%
      column_spec (1:10, latex_column_spec = "l|ccc|ccc|ccc") %>%
      kable_styling(latex_options="scale_down") %>% 
      save_kable(filename)
  }
  
  
}