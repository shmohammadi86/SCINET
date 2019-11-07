#' SCINET (Single-Cell Imputation and NETwork construction)
#'
#' An analytical framework to construct context-specific interactomes at a single-cell resolution
#'
#' @section Usage: 
#' 
#' \enumerate{
#' \item ?compute_gene_activities_full to compute activity scores for cell states computed by ACTIONet
#' \item ?pair.datasets to match rows of expression matrix and input graph
#' \item ?construct_cell_networks to construct networks
#' }
#' @section Useful links:
#' 
#' \enumerate{
#' \item Report bugs at \url{https://github.com/shmohammadi86/SCINET/issues}
#' \item Read the manuscript: 
#' \href{https://www.biorxiv.org/content/10.1101/586859v1}{Single-cell interactomes of the human brain reveal cell-type specific convergence of brain disorders}.
#' }
#' 
#'
#' @name SCINET
#' @docType package
#' @importFrom Rcpp evalCpp
#' @exportPattern ^[[:alpha:]]+ 
#' @useDynLib SCINET, .registration=TRUE
NULL

