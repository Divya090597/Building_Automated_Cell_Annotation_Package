#' Generate Models and Initialize New Columns List
#'
#' This function creates a list of models and initializes an empty list to store
#' new columns that can be added or processed later.
#'
#' @param model_names A character string representing the name of the model.
#' @param model_objects The object representing the model (e.g., gene list or matrix).
#' @return A list containing the models and an initialized empty list for new columns.
#' @export
generate_model_and_columns <- function(model_names = character(),
                                       model_objects = list()) {

  models <- list()

  for (i in seq_along(model_names)) {
    models[[model_names[i]]] <- model_objects[[i]]
  }


  new_columns_list <- list()


  return(list(models = models, new_columns_list = new_columns_list))
}
