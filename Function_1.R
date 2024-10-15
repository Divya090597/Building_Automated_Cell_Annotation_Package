#'
#' Load Macrophage Gene List Model
#'
#' This function loads a gene list model from a file, specifically for myeloid cells.
#' The model is typically used for cell classification or gating in single-cell transcriptomics data.
#'
#' @param file_path A string specifying the file path to the TSV file containing the myeloid cell gene list.
#'
#' @return An object representing the loaded gene list model. Typically, this will be a list containing gene signatures, levels, thresholds and annotations for cell classification.
#'
#'
#' @export
load_mac_model <- function(file_path) {

  gene_list_model <- scGate::load_scGate_model(file_path)

  return(gene_list_model)
}



