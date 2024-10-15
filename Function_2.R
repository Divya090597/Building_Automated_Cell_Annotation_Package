#' Load Macrophage Model from a URL
#'
#' This function downloads a gene list model from a provided URL and loads it using
#'
#' @param Macrophage A character string representing the URL of the gene list model file. The file
#'   should be accessible online and in the proper format for `scGate::load_scGate_model()`.
#' @return The loaded macrophage model.
#' @importFrom utils download.file
#' @export

load_model <- function(Macrophage) {

  temp_file <- tempfile(fileext = ".tsv")

  download.file(Macrophage, destfile = temp_file, mode = "wb", method = "curl", extra = "-k")

  mac_model <- scGate::load_scGate_model(temp_file)

  unlink(temp_file)

  return(mac_model)
}


