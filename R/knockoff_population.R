#' Knockoff Group Analysis for Population Design
#'
#' This file contains functions for performing knockoff group analysis in population designs.
#' The main function \code{KS.chr_ko} computes knockoff statistics for population data.
#' Auxiliary functions include preliminary model fitting, p-value calculations, and more.
#'
#' @name knockoff_population_functions
#' @import Matrix
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp glm gaussian binomial var pchisq rbinom dbeta na.omit
#' @importFrom SPAtest ScoreTest_SPA
NULL

#' Preliminary Model for Analysis
#'
#' Fit a null model for the outcome (continuous or dichotomous).
#'
#' @param Y Outcome vector.
#' @param X Covariate matrix (optional).
#' @param id ID vector (optional).
#' @param out_type Outcome type: "C" for continuous, "D" for dichotomous.
#' @return A list containing fitted model components (residuals, inverse matrices, etc.).
#' @export
KS.prelim <- function(Y, X = NULL, id = NULL, out_type = "C") {
  Y <- as.matrix(Y)
  n <- nrow(Y)

  if (length(X) != 0) {
    X0 <- svd(as.matrix(X))$u
  } else {
    X0 <- NULL
  }
  X0 <- cbind(rep(1, n), X0)

  if (out_type == "C") {
    glm.null <- glm(Y ~ 0 + X0, family = gaussian)
  } else if (out_type == "D") {
    glm.null <- glm(Y ~ 0 + X0, family = binomial)
  } else {
    stop("out_type must be 'C' or 'D'")
  }

  mu <- glm.null$fitted.values
  Y.res <- Y - mu

  # Calculate variance
  if (out_type == "D") {
    v <- mu * (1 - mu)
  } else {
    v <- rep(var(Y.res), n)
  }

  # Pre-compute inverse matrices for score tests
  inv.X0  <- solve(t(X0) %*% (v * X0))
  inv.vX0 <- solve(t(X0) %*% (v * X0))

  if (length(id) == 0) id <- seq_len(n)

  list(Y = Y, id = id, n = n, X0 = X0, nullglm = glm.null, out_type = out_type,
       inv.X0 = inv.X0, inv.vX0 = inv.vX0)
}

#' Compute Single Variant P-Values
#'
#' Compute p-values for single variants using SPA (for binary) or linear score test (for continuous).
#'
#' @param X Genotype matrix (n_samples x n_variants).
#' @param result.prelim Preliminary model result from \code{KS.prelim}.
#' @return Vector of p-values.
#' @noRd
Get.p <- function(X, result.prelim) {
  X <- as.matrix(X)

  mu  <- result.prelim$nullglm$fitted.values
  Y.res <- result.prelim$Y - mu

  if (result.prelim$out_type == "D") {
    ScoreTest_SPA(t(X), result.prelim$Y, result.prelim$X0,
                  method = "fastSPA", minmac = -Inf)$p.value
  } else {
    v <- rep(var(Y.res), nrow(X))
    A <- (t(X) %*% Y.res)^2
    B <- colSums(v * X^2)
    C <- t(X) %*% (v * result.prelim$X0) %*% result.prelim$inv.X0
    D <- t(t(result.prelim$X0) %*% (v * X))

    # Calculate Score Statistic
    stat <- A / (B - rowSums(C * D))
    stat[stat < 0] <- 0
    pchisq(stat, df = 1, lower.tail = FALSE)
  }
}

#' Impute Missing Values
#'
#' Impute missing genotypes based on MAF.
#'
#' @param Z Genotype matrix.
#' @param method Imputation method: "fixed", "random", "bestguess".
#' @return Imputed matrix.
#' @noRd
Impute <- function(Z, method = "fixed") {
  Z <- as.matrix(Z)
  if (sum(is.na(Z)) == 0) return(Z)

  for (j in seq_len(ncol(Z))) {
    na.idx <- is.na(Z[, j])
    if (!any(na.idx)) next

    maf <- mean(Z[!na.idx, j]) / 2

    if (method == "fixed")      Z[na.idx, j] <- 2 * maf
    if (method == "random")     Z[na.idx, j] <- rbinom(sum(na.idx), 2, maf)
    if (method == "bestguess")  Z[na.idx, j] <- round(2 * maf)
  }
  Z
}

#' Cauchy Combination
#'
#' Combine p-values using the Cauchy combination method (CCT).
#'
#' @param p Vector of p-values.
#' @return Combined p-value.
#' @noRd
Get.cauchy <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) == 0) return(NA)

  p[p > 0.99] <- 0.99

  tmp <- rep(NA, length(p))
  small <- p < 1e-16
  reg   <- p >= 1e-16

  tmp[small] <- 1 / p[small] / pi
  tmp[reg]   <- tan((0.5 - p[reg]) * pi)

  cct <- mean(tmp)

  if (is.na(cct)) return(NA)
  if (cct > 1e15) return(1 / cct / pi)

  1 - pcauchy(cct)
}

#' Cauchy Scan for Windows
#'
#' Compute Cauchy p-values aggregated over windows.
#'
#' @param p Vector of p-values.
#' @param window.matrix Sparse matrix indicating window membership.
#' @return Vector of window-based p-values.
#' @noRd
Get.cauchy.scan <- function(p, window.matrix) {
  p[p > 0.99] <- 0.99

  tmp <- numeric(length(p))
  small <- p < 1e-16 & !is.na(p)

  tmp[small] <- 1 / p[small] / pi
  tmp[!small] <- tan((0.5 - p[!small]) * pi)

  # Aggregation via matrix multiplication
  cct <- as.numeric(t(tmp) %*% window.matrix / colSums(window.matrix))

  pval <- numeric(length(cct))
  large <- cct > 1e15 & !is.na(cct)

  pval[large] <- 1 / cct[large] / pi
  pval[!large] <- 1 - pcauchy(cct[!large])

  pval
}

#' SKAT/Dispersion P-Values
#'
#' Compute SKAT (Variance component) p-values.
#'
#' @param score Score vector.
#' @param K Kernel matrix part.
#' @param window.matrix Window membership matrix.
#' @param weight Weight vector.
#' @param result.prelim Preliminary model.
#' @return Matrix of p-values.
#' @noRd
Get.p.SKAT <- function(score, K, window.matrix, weight, result.prelim) {
  Q <- as.vector(t(score^2) %*% (weight * window.matrix)^2)
  K.temp <- weight * t(weight * K)
  p <- rep(NA, length(Q))

  for (i in seq_along(Q)) {
    win.idx <- window.matrix[, i] != 0
    if (sum(win.idx) == 0) next

    tmp <- K.temp[win.idx, win.idx, drop = FALSE]
    if (sum(tmp^2) == 0) next

    lambda <- eigen(tmp, symmetric = TRUE, only.values = TRUE)$values
    if (anyNA(lambda)) next

    # Match the dispersion test used by the manuscript analysis.
    probability <- CompQuadForm::davies(Q[i], lambda, acc = 1e-3)$Qq
    if (probability > 1 || probability <= 0) {
      probability <- Get_Liu_PVal.MOD.Lambda(Q[i], lambda)
    }
    p[i] <- probability
  }
  as.matrix(p)
}

#' Liu Parameters for Lambda
#' @noRd
Get_Liu_Params_Mod_Lambda <- function(lambda) {
  c1 <- sapply(1:4, function(i) sum(lambda^i))
  muQ <- c1[1]
  sigmaQ <- sqrt(2 * c1[2])
  s1 <- c1[3] / c1[2]^(3/2)
  s2 <- c1[4] / c1[2]^2

  if (s1^2 > s2) {
    a <- 1 / (s1 - sqrt(s1^2 - s2))
    d <- s1 * a^3 - a^2
    l <- a^2 - 2 * d
  } else {
    l <- 1 / s2; a <- sqrt(l); d <- 0
  }
  list(l = l, d = d, muQ = muQ, muX = l + d, sigmaQ = sigmaQ, sigmaX = sqrt(2) * a)
}

#' Liu P-Value
#' @noRd
Get_Liu_PVal.MOD.Lambda <- function(Q, lambda) {
  param <- Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm <- (Q - param$muQ) / param$sigmaQ
  Q.Norm1 <- Q.Norm * param$sigmaX + param$muX
  pchisq(Q.Norm1, df = param$l, ncp = param$d, lower.tail = FALSE)
}

#' Core KS Test
#' @noRd
KS.test <- function(G, G_k, p.single, p.single_k, result.prelim,
                    window.matrix, weight.matrix) {
  mu  <- result.prelim$nullglm$fitted.values
  Y.res <- result.prelim$Y - mu
  outcome <- result.prelim$out_type
  M <- length(G_k)

  v <- if (outcome == "D") mu * (1 - mu) else rep(var(Y.res), nrow(G))

  # Ensure G is dense for crossprod to avoid S4 dispatch issues if not loaded
  G_dense <- as.matrix(G)

  A <- crossprod(G_dense, v * G_dense)
  B <- crossprod(G_dense, v * result.prelim$X0)
  C <- result.prelim$inv.vX0
  K <- A - B %*% C %*% t(B)

  score    <- t(G_dense) %*% Y.res
  score_k <- matrix(vapply(G_k, function(gk) as.numeric(crossprod(as.matrix(gk), Y.res)),
                           numeric(ncol(G))), nrow = ncol(G), ncol = M)

  # Burden Test
  wm <- weight.matrix[, 1] * window.matrix
  # Use dense matrix multiplication
  p.burden      <- Get.p(as.matrix(G %*% wm), result.prelim)
  p.burden_k    <- matrix(vapply(G_k, function(gk)
    as.numeric(Get.p(as.matrix(gk %*% wm), result.prelim)), numeric(ncol(window.matrix))),
    nrow = ncol(window.matrix), ncol = M)

  # Dispersion Test
  p.dispersion    <- Get.p.SKAT(score, K, window.matrix, weight.matrix[, 1], result.prelim)
  p.dispersion_k <- matrix(vapply(seq_len(M), function(m)
    as.numeric(Get.p.SKAT(score_k[, m], K, window.matrix, weight.matrix[, 1], result.prelim)),
    numeric(ncol(window.matrix))), nrow = ncol(window.matrix), ncol = M)

  # Single-variant Cauchy
  p.V1    <- Get.cauchy.scan(p.single, window.matrix)
  p.V1_k <- matrix(apply(p.single_k, 2, Get.cauchy.scan, window.matrix = window.matrix),
                   nrow = ncol(window.matrix), ncol = M)

  # Combine
  p.individual <- cbind(burden = p.burden, dispersion = p.dispersion, singleCauchy = p.V1)
  p.KS <- apply(p.individual, 1, Get.cauchy)

  p.KS_k <- matrix(NA, nrow(p.individual), M)
  for (m in 1:M) {
    tmp <- cbind(burden = p.burden_k[, m],
                 dispersion = p.dispersion_k[, m],
                 singleCauchy = p.V1_k[, m])
    p.KS_k[, m] <- apply(tmp, 1, Get.cauchy)
  }

  list(p.KS = p.KS, p.KS_k = p.KS_k, p.individual = p.individual)
}

#' Main KS Function for Chromosome
#'
#' Compute knockoff statistics for a chromosome in population design.
#'
#' @param result.prelim Preliminary model fitted by \code{KS.prelim}.
#' @param chr Chromosome identifier.
#' @param dat_ko Knockoff genotype data (Array: n x p x M).
#' @param dat Original genotype data (Matrix: n x p).
#' @param pos Position vector.
#' @param window.bed Matrix defining windows (chr, start, end).
#' @param excluded.pos Optional excluded positions; retained for compatibility.
#' @param M Number of knockoffs (default 10).
#' @param thres.single Threshold for single variant output.
#' @param maxcor Max correlation (unused).
#' @param maxbp Max base pairs (unused).
#' @param thres.ultrarare Threshold for ultrarare variants.
#' @param thres.missing Missingness threshold; retained for compatibility.
#' @param midout.dir Optional output directory; retained for compatibility.
#' @param temp.dir Optional temporary directory; retained for compatibility.
#' @param jobtitle Optional job title; retained for compatibility.
#' @param Gsub.id Optional sample subset ID; retained for compatibility.
#' @param impute.method Imputation method ("fixed", "random", "bestguess").
#' @param bigmemory Logical flag retained for compatibility.
#' @param leveraging Logical flag retained for compatibility.
#' @param LD.filter Optional LD filter retained for compatibility.
#' @return A list containing results for single variants and windows.
#' @export
KS.chr_ko <- function(result.prelim, chr, dat_ko, dat, pos, window.bed,
                      excluded.pos = NULL, M = 10, thres.single = 0.01,
                      maxcor = 0, maxbp = 5e5, thres.ultrarare = 0,
                      thres.missing = 0.10, midout.dir = NULL, temp.dir = NULL,
                      jobtitle = NULL, Gsub.id = NULL, impute.method = "fixed",
                      bigmemory = TRUE, leveraging = TRUE, LD.filter = NULL) {

  # 1. Preprocessing
  G <- dat
  G[G < 0 | G > 2] <- NA
  G <- Impute(G, impute.method)

  MAF <- colMeans(G) / 2
  MAF[is.na(MAF)] <- 0
  flip <- MAF > 0.5
  MAF[flip] <- 1 - MAF[flip]
  if (any(flip)) G[, flip] <- 2 - G[, flip]

  # Convert to sparse for storage efficiency
  G <- Matrix(G, sparse = TRUE)

  # 2. Single Variant P-values (Original)
  # G is sparse here, but Get.p now handles it via as.matrix() internally
  p.single <- Get.p(G, result.prelim)

  # 3. Pre-process Knockoff Data
  G_k <- lapply(seq_len(M), function(m)
    matrix(dat_ko[, , m, drop = FALSE], nrow = nrow(G), ncol = ncol(G)))

  # 4. Knockoff Single Variant P-values (Chunked)
  p.single_k <- matrix(NA, ncol(G), M)

  for (m in 1:M) {
    tmp <- matrix(0, ncol(G), 1)
    chunk_size <- 1000
    for (i in seq(1, ncol(G), by = chunk_size)) {
      j <- min(i + chunk_size - 1, ncol(G))
      # Use drop=FALSE to keep matrix structure
      subset_G_k <- G_k[[m]][, i:j, drop = FALSE]
      # Get.p will convert this chunk to dense matrix
      tmp[i:j, ] <- Get.p(subset_G_k, result.prelim)
    }
    p.single_k[, m] <- tmp
  }

  # ---------- Output: Single Variant ----------
  result.summary.single <- NULL
  comm.idx <- which(MAF >= thres.single)

  if (length(comm.idx) > 0) {
    pc    <- p.single[comm.idx]
    pc_k  <- p.single_k[comm.idx, , drop = FALSE]

    stats_neg_log <- -log10(pc)
    stats_neg_log_k <- -log10(pc_k)

    W <- (stats_neg_log - apply(stats_neg_log_k, 1, median)) *
      (stats_neg_log >= apply(stats_neg_log_k, 1, max))
    W[is.na(W)] <- 0

    mk <- MK.statistic(stats_neg_log, stats_neg_log_k)
    result.summary.single <- cbind(chr, pos[comm.idx], pos[comm.idx],
                                   pos[comm.idx], pos[comm.idx],
                                   mk, W, pc, pc_k, MAF[comm.idx])
    colnames(result.summary.single) <- c("chr", "start", "end", "actual_start", "actual_end",
                                         "kappa", "tau", "W_KS", "P_KS",
                                         paste0("P_KS_k", 1:M), "MAF")
  }

  # ---------- Output: Window ----------
  weight.beta <- 1 / sqrt(nrow(G) * MAF * (1 - MAF))
  weight.beta[is.infinite(weight.beta)] <- 0
  weight.matrix <- matrix(weight.beta, ncol = 1)
  colnames(weight.matrix) <- "w1"

  win.tmp <- window.bed[window.bed[, 2] <= max(pos) & window.bed[, 3] >= min(pos), , drop = FALSE]
  if (nrow(win.tmp) == 0) stop("No overlapping windows found.")

  win.mat0 <- matrix(vapply(seq_len(nrow(win.tmp)), function(i)
    as.numeric(pos >= win.tmp[i, 2] & pos < win.tmp[i, 3]), numeric(length(pos))),
    nrow = length(pos), ncol = nrow(win.tmp))

  # Force conversion to dense for colSums to avoid ambiguity
  win.MAC  <- colSums(colSums(as.matrix(G)) * win.mat0)
  membership <- apply(win.mat0, 2, paste0, collapse = "")
  keep <- which(!duplicated(membership) & colSums(win.mat0) > 0 & win.MAC >= 10)
  if (length(keep) == 0L) {
    return(list(result.single = result.summary.single, result.window = data.frame()))
  }

  window.matrix  <- Matrix(win.mat0[, keep, drop = FALSE], sparse = TRUE)
  actual_ranges <- t(apply(win.mat0[, keep, drop = FALSE], 2, function(x) {
    idx <- which(x == 1)
    if(length(idx) > 0) c(min(pos[idx]), max(pos[idx])) else c(NA, NA)
  }))
  window.summary <- cbind(win.tmp[keep, 2:3, drop = FALSE], actual_ranges)

  # Run KS Test
  KS.fit <- KS.test(G, G_k, p.single, p.single_k, result.prelim,
                    window.matrix = window.matrix, weight.matrix = weight.matrix)

  nlog_p <- -log10(KS.fit$p.KS)
  nlog_p_k <- -log10(KS.fit$p.KS_k)
  W <- (nlog_p - apply(nlog_p_k, 1, median)) * (nlog_p >= apply(nlog_p_k, 1, max))
  W[is.na(W)] <- 0
  mk <- MK.statistic(nlog_p, nlog_p_k)

  result.summary.window <- cbind(
    chr, window.summary, mk, W_KS = W, P_KS = KS.fit$p.KS,
    KS.fit$p.KS_k, KS.fit$p.individual
  )
  colnames(result.summary.window) <- c(
    "chr", "start", "end", "actual_start", "actual_end",
    "kappa", "tau", "W_KS", "P_KS",
    paste0("P_KS_k", 1:M), "burden", "dispersion", "singleCauchy"
  )

  list(result.single = result.summary.single, result.window = result.summary.window)
}
