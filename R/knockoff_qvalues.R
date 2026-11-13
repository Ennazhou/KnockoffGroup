#' Calculate Multiple-Knockoff q-values
#'
#' Apply the multiple-knockoff threshold across the complete set of tested genes.
#' Tied scores enter the threshold together. Original winners have kappa = 0.
#' @param kappa Integer winner indices, from 0 to M.
#' @param tau Non-negative importance gaps, one per gene.
#' @param M Number of knockoff copies.
#' @return A numeric vector of q-values in input order. Missing statistics,
#'   zero gaps, and knockoff winners receive q = 1.
#' @export
#' @examples
#' knockoff_qvalues(c(0, 0, 1), c(3, 2, 1), M = 10)
knockoff_qvalues <- function(kappa, tau, M) {
    if (length(M) != 1L || is.na(M) || M < 1 || M != as.integer(M))
        stop("M must be a positive integer.")
    if (length(kappa) != length(tau)) stop("kappa and tau must have the same length.")
    valid <- !is.na(kappa) & !is.na(tau)
    if (any(kappa[valid] < 0 | kappa[valid] > M | kappa[valid] != floor(kappa[valid])) ||
        any(tau[valid] < 0)) stop("Invalid kappa or tau.")
    q <- rep(1, length(tau))
    order_index <- which(valid & tau > 0)
    order_index <- order_index[order(tau[order_index], decreasing = TRUE)]
    if (!length(order_index)) return(q)
    winners <- kappa[order_index] == 0
    positive <- cumsum(winners)
    negative <- cumsum(!winners)
    threshold_end <- which(!duplicated(tau[order_index], fromLast = TRUE))
    ratio <- (1 + negative[threshold_end]) / (M * pmax(positive[threshold_end], 1))
    threshold_q <- pmin(1, rev(cummin(rev(ratio))))
    group <- match(tau[order_index], tau[order_index][threshold_end])
    q[order_index[winners]] <- threshold_q[group[winners]]
    q
}
