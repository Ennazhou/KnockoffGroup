#' Unified KnockoffGroup Function
#'
#' Unified interface for knockoff group analysis, supporting both trio and population designs.
#'
#' @param design "trio" or "population".
#' @param ... Arguments passed to the design-specific function.
#' @return Results from the design-specific function.
#' @export
#'
#' @examples
#' \dontrun{
#' library(KnockoffGroup)
#'
#' # Trio example
#' data(KnockoffGroup.example)
#' dat.ko.trio <- create_knockoff(KnockoffGroup.example$window_matrix,
#'                                KnockoffGroup.example$gene.clusters,
#'                                KnockoffGroup.example$dat.hap,
#'                                KnockoffGroup.example$pos,
#'                                design = "trio", M = 10)
#'
#' window.result.trio <- KnockoffGroup(design = "trio",
#'                                     KnockoffGroup.example$geneinfo,
#'                                     KnockoffGroup.example$dat,
#'                                     dat.ko.trio,
#'                                     KnockoffGroup.example$pos,
#'                                     chr = as.character(1))
#'
#' # Population example
#' dat.ko.pop <- create_knockoff(KnockoffGroup.example$window_matrix,
#'                               KnockoffGroup.example$gene.clusters,
#'                               KnockoffGroup.example$dat.pop,
#'                               KnockoffGroup.example$pos,
#'                               design = "population",
#'                               M = 10)
#'
#' y.b <- KnockoffGroup.example$y.b
#' X <- as.matrix(rnorm(nrow(KnockoffGroup.example$dat.pop), 0, 1))
#' result.prelim <- KS.prelim(y.b, X = X, out_type = "D")
#'
#' window.results.pop <- KnockoffGroup(
#'   design        = "population",
#'   result.prelim = result.prelim,
#'   chr           = 1,
#'   dat_ko        = dat.ko.pop,
#'   dat           = KnockoffGroup.example$dat.pop,
#'   pos           = KnockoffGroup.example$pos,
#'   window.bed    = KnockoffGroup.example$window.bed,
#'   M             = 10,
#'   thres.single  = 0.01,
#'   maxcor        = 0,
#'   maxbp         = 5e5,
#'   thres.ultrarare = 0
#' )
#' }
#'
KnockoffGroup <- function(design = "trio", ...) {
  if (design == "trio") {
    return(KnockoffGroup_trio(...))
  } else if (design == "population") {
    return(KS.chr_ko(...))
  } else {
    stop("Invalid design: must be 'trio' or 'population'")
  }
}
