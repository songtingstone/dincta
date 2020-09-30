#' Dincta: Data INtegration and Cell Type Annotation of Single Cell Transcriptomes
#' 
#' Algorithm for single cell integration and cell type annotation. 
#'
#' @section Usage: 
#' 
#' \enumerate{
#' \item ?DinctaMatrix to run Dincta on gene expression or PCA 
#' embeddings matrix.
#' \item ?DinctaClassification to run Dincta Classification on on gene expression or PCA 
#' embeddings matrix.
#' }
#' @section Useful links:
#' 
#' \enumerate{
#' \item Report bugs at \url{https://github.com/songtingstone/dincta/issues}
#' \item Read the manuscript
#' \href{https://www.biorxiv.org/content/tobeappear}{online}.
#' }
#' 
#'
#' @name dincta
#' @docType package
#' @useDynLib dincta
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp loadModule
#' @importFrom methods new
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
loadModule("dincta_full_module", TRUE)
loadModule("dincta_marginal_module", TRUE)
loadModule("dincta_classification_module", TRUE)
NULL
