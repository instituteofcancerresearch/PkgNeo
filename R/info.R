#' Returns the genes that have pre-saved dna
#'
#' @description This function is for the neoantigens project.
#' We call it with a mutation, gene and sequence.
#'
#' @return a list of genes with pre-saved dna sequences
#' @export
#' @examples
#' get_pre_saved_genes_list()
#' @references
#' \url{https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html}
#' @author Rachel Alcraft, \url{https://github.com/instituteofcancerresearch/PkgStopGain}
#'


get_pre_saved_genes_list <- function() {
  return(c("BRCA1", "BRCA2", "RAD51C", "RAD51D", "PALB2"))
}