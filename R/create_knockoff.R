#' Create Knockoff Genotypes
#'
#' This file contains functions for creating integrated knockoff genotype data for trio or population designs.
#' The main function `create_knockoff` generates knockoffs using phased haplotype data or unphased genotype data.
#'
#' @name create_knockoff_functions
NULL

#' Create integrated knockoff genotype data
#'
#' Create knockoff genotype data using phased haplotype data for trio or population design.
#' This function integrates the trio-based and population-based designs.
#'
#' @param window_matrix A matrix defining window boundaries for groups.
#' @param gene.clusters A vector of cluster assignments for genes/variants.
#' @param dat Data input: For 'trio' design, a 6n*p haplotype matrix (n trios, p variants, 2 haplotypes per individual in trio order: father, mother, offspring). For 'population' design, an n*p genotype or haplotype matrix (n individuals, p variants).
#' @param pos A numeric vector of length p for variant positions.
#' @param design Character specifying the design: "trio" (default) or "population".
#' @param M Positive integer for the number of knockoffs (default: 10).
#' @param cormin Real number in [0,1] for the minimum correlation threshold (default: 0.01).
#' @param phasing.dad Numeric vector for paternal phasing (required for trio if not NULL; auto-calculated if NULL).
#' @param phasing.mom Numeric vector for maternal phasing (required for trio if not NULL; auto-calculated if NULL).
#' @param maxbp Positive integer for the size of neighboring base pairs (default: 500000).
#' @param maxN.neighbor Maximum number of neighboring variants to use in each
#'   conditional regression. The effective cap is also bounded by
#'   \code{floor(nrow(X)^(1/3))}, matching the runtime fast implementation.
#'
#' @return For trio design: A 3n*p*M array for the knockoff genotype data. For population design: An n*p*M array for the knockoff data.
#'
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp dbeta na.omit
#' @export
#'
#' @examples
#' # Trio example
#' data(KnockoffGroup.example)
#' dat.ko.trio <- create_knockoff(KnockoffGroup.example$window_matrix,
#'                                KnockoffGroup.example$gene.clusters,
#'                                KnockoffGroup.example$dat.hap,
#'                                KnockoffGroup.example$pos, design="trio", M=10)
#' # Population example
#' dat.ko.pop <- create_knockoff(KnockoffGroup.example$window_matrix,
#'                               KnockoffGroup.example$gene.clusters,
#'                               KnockoffGroup.example$dat.pop,
#'                               KnockoffGroup.example$pos,
#'                               design="population", M=10)
create_knockoff <- function(window_matrix, gene.clusters, dat, pos, design = "trio", M=10, cormin=0.01,
                            phasing.dad=NULL, phasing.mom=NULL, maxbp=500000,
                            maxN.neighbor=60) {
  if (design == "trio") {
    if (nrow(dat) %% 6 != 0) stop("The number of rows of the input haplotype matrix must be a multiple of six")
    if (!all(dat %in% c(0,1))) stop("The input haplotype matrix can only contain 0 or 1.")
    n <- nrow(dat)/6 # number of trios
    nsnp <- ncol(dat) # number of variants

    ind.hap.off <- c((1:n)*6-1, (1:n)*6)
    if (is.null(phasing.dad) | is.null(phasing.mom)) {
      phasing <- get_phasing(dat)
      phasing.dad <- phasing$phasing.dad
      phasing.mom <- phasing$phasing.mom
    }
    phasing.dad <- phasing.dad %% 2
    phasing.mom <- phasing.mom %% 2
    # remove offspring
    dat.hap <- dat[-ind.hap.off,, drop=FALSE]

    dad_hap <- dat.hap[sort(c(seq(1,(4*n),4), seq(2,(4*n),4))),, drop=FALSE]
    ind_dad <- 2*(1:n) + sample(c(0,-1), n, replace = T)
    ind_dad2 <- (1:(2*n))[-ind_dad]
    dad1 <- dad_hap[ind_dad,, drop=FALSE]
    dad2 <- dad_hap[ind_dad2,, drop=FALSE]

    mom_hap <- dat.hap[sort(c(seq(3,(4*n),4), seq(4,(4*n),4))),, drop=FALSE]
    ind_mom <- 2*(1:n) + sample(c(0,-1), n, replace = T)
    ind_mom2 <- (1:(2*n))[-ind_mom]
    mom1 <- mom_hap[ind_mom,, drop=FALSE]
    mom2 <- mom_hap[ind_mom2,, drop=FALSE]

    dadk1 <- create_knockoff_core(window_matrix, gene.clusters, dad1, pos, M, cormin, maxbp, maxN.neighbor)
    dadk2 <- create_knockoff_core(window_matrix, gene.clusters, dad2, pos, M, cormin, maxbp, maxN.neighbor)
    momk1 <- create_knockoff_core(window_matrix, gene.clusters, mom1, pos, M, cormin, maxbp, maxN.neighbor)
    momk2 <- create_knockoff_core(window_matrix, gene.clusters, mom2, pos, M, cormin, maxbp, maxN.neighbor)
    dadk <- dadk1 + dadk2
    momk <- momk1 + momk2
    kidk <- array(dim=dim(dadk))

    for (i in 1:n) {
      if (ind_dad[i] %% 2 == phasing.dad[i]) hap1 <- dadk1[i,,, drop=FALSE] else hap1 <- dadk2[i,,, drop=FALSE]
      if (ind_mom[i] %% 2 == phasing.mom[i]) hap2 <- momk1[i,,, drop=FALSE] else hap2 <- momk2[i,,, drop=FALSE]
      kidk[i,,] <- hap1 + hap2
    }

    dat1 <- array(dim=c(3*n, nsnp, M)) # knockoff for trio, dim=(3*n,nsnp,M)
    for (i in 1:n) {
      dat1[i*3-2,,] <- dadk[i,,]
      dat1[i*3-1,,] <- momk[i,,]
      dat1[i*3,,] <- kidk[i,,]
    }
    return(dat1)
  } else if (design == "population") {
    X <- as.matrix(dat)
    return(create_knockoff_core(window_matrix, gene.clusters, X, pos, M, cormin, maxbp, maxN.neighbor))
  } else {
    stop("Invalid design: must be 'trio' or 'population'")
  }
}

#' @title Core Knockoff Creation Function (Internal)
#' @description Unified core function for creating knockoffs, used in both designs.
#' @noRd
create_knockoff_core <- function(window_matrix, gene.clusters, X, pos, M=6, cormin=0.01,
                                 maxbp=500000, maxN.neighbor=60) {
  gene.clusters[is.na(gene.clusters)] <- max(gene.clusters, na.rm = TRUE) + seq_len(sum(is.na(gene.clusters)))

  X <- as.matrix(X)
  sparse.fit <- sparse.cor(X)
  cor.X <- sparse.fit$cor
  cor.X[is.na(cor.X)] <- 0
  cor.X[is.infinite(cor.X)] <- 0
  cov.X <- sparse.fit$cov
  cov.X[is.na(cov.X)] <- 0
  cov.X[is.infinite(cov.X)] <- 0

  clusters <- na.omit(unique(gene.clusters))

  X_k <- array(0, dim=c(nrow(X), ncol(X), M)); index.exist <- c()
  hard_cap <- min(maxN.neighbor, floor(nrow(X)^(1 / 3)))
  for (k in unique(gene.clusters)) {
    cluster.fitted <- cluster.residuals <- matrix(NA, nrow(X), sum(gene.clusters==k, na.rm=TRUE))
    indk <- which(gene.clusters==k)
    ind_gene <- which(pos >= window_matrix$start[k] & pos <= window_matrix$end[k])

    for (i in indk) {
      temp <- abs(cor.X[i,]); temp[indk] <- 0; temp[ind_gene] <- 0
      index.pos <- which(pos >= max(pos[i]-maxbp, pos[1]) & pos <= min(pos[i]+maxbp, pos[length(pos)]))
      temp[-index.pos] <- 0

      index <- order(temp, decreasing=T)
      max_take <- min(length(index), sum(temp>cormin), hard_cap)
      index <- setdiff(index[1:max_take], i)

      y <- X[,i]
      if (length(index)==0) { fitted.values <- mean(y) } else {
        x <- X[,index, drop=F]; temp.xy <- rbind(mean(y), crossprod(x,y)/length(y) - colMeans(x)*mean(y))
        x.exist <- c()
        for (j in 1:M) {
          x.exist <- cbind(x.exist, X_k[,intersect(index, index.exist), j])
        }
        temp.xy <- rbind(temp.xy, crossprod(x.exist,y)/length(y) - colMeans(x.exist)*mean(y))

        temp.cov.cross <- sparse.cov.cross(x, x.exist)$cov
        temp.cov <- sparse.cor(x.exist)$cov
        temp.xx <- cov.X[index, index]
        temp.xx <- rbind(cbind(temp.xx, temp.cov.cross), cbind(t(temp.cov.cross), temp.cov))

        temp.xx <- cbind(0, temp.xx)
        temp.xx <- rbind(c(1, rep(0, ncol(temp.xx)-1)), temp.xx)

        pca.fit <- princomp(covmat=temp.xx)
        v <- pca.fit$loadings
        cump <- cumsum(pca.fit$sdev^2)/sum(pca.fit$sdev^2)
        n.pc <- which(cump >= 0.999)[1]
        pca.index <- intersect(1:n.pc, which(pca.fit$sdev != 0))

        temp.inv <- v[,pca.index, drop=F] %*% (pca.fit$sdev[pca.index]^(-2) * t(v[,pca.index, drop=F]))
        temp.beta <- temp.inv %*% temp.xy

        temp.j <- 1
        fitted.values <- temp.beta[1] + crossprod(t(x), temp.beta[(temp.j+1):(temp.j+ncol(x)),, drop=F]) - sum(colMeans(x)*temp.beta[(temp.j+1):(temp.j+ncol(x)),, drop=F])
        temp.j <- temp.j + ncol(x)
        for (j in 1:M) {
          temp.x <- as.matrix(X_k[,intersect(index, index.exist), j])
          if (ncol(temp.x) >= 1) {
            fitted.values <- fitted.values + crossprod(t(temp.x), temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),, drop=F]) - sum(colMeans(temp.x)*temp.beta[(temp.j+1):(temp.j+ncol(temp.x)),, drop=F])
          }
          temp.j <- temp.j + ncol(temp.x)
        }
      }
      residuals <- y - fitted.values
      cluster.fitted[, match(i, indk)] <- as.vector(fitted.values)
      cluster.residuals[, match(i, indk)] <- as.vector(residuals)

      index.exist <- c(index.exist, i)
    }
    cluster.sample.index <- sapply(1:M, function(x) sample(1:nrow(X)))
    for (j in 1:M) {
      X_k[, indk, j] <- cluster.fitted + cluster.residuals[cluster.sample.index[,j],, drop=F]
    }
  }
  return(X_k)
}

#' @title Get Phasing Information
#' @description Calculate phasing for haplotypes (used in trio design).
#' @param dat.hap Haplotype matrix.
#' @param xchr Logical for X chromosome.
#' @return List with phasing.dad and phasing.mom.
#' @noRd
get_phasing <- function(dat.hap, xchr=FALSE) {
  nsnp <- ncol(dat.hap)
  ntrio <- nrow(dat.hap)/6
  phasing.dad <- phasing.mom <- rep(1, ntrio)
  if (!xchr) {
    for (i in 1:ntrio) {
      index <- ((i-1)*6+1):(i*6)
      if (length(which(dat.hap[index[1],]==dat.hap[index[5],]))/nsnp >= length(which(dat.hap[index[2],]==dat.hap[index[5],]))/nsnp) phasing.dad[i] <- 1 else phasing.dad[i] <- 2
      if (length(which(dat.hap[index[3],]==dat.hap[index[6],]))/nsnp >= length(which(dat.hap[index[4],]==dat.hap[index[6],]))/nsnp) phasing.mom[i] <- 1 else phasing.mom[i] <- 2
    }
  } else {
    for (i in 1:ntrio) {
      index <- ((i-1)*6+1):(i*6)
      momhap <- 6
      if (length(which(dat.hap[index[1],]==dat.hap[index[5],]))/nsnp < length(which(dat.hap[index[1],]==dat.hap[index[6],]))/nsnp) momhap <- 5
      if (length(which(dat.hap[index[3],]==dat.hap[index[momhap],]))/nsnp >= length(which(dat.hap[index[4],]==dat.hap[index[momhap],]))/nsnp) phasing.mom[i] <- 1 else phasing.mom[i] <- 2
    }
  }
  out <- list()
  out$phasing.dad <- phasing.dad
  out$phasing.mom <- phasing.mom
  return(out)
}

#' @title Sparse Correlation
#' @description Compute sparse correlation and covariance.
#' @param x Matrix.
#' @return List with cov and cor.
#' @noRd
sparse.cor <- function(x) {
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat, cor=cormat)
}

#' @title Sparse Cross Covariance
#' @description Compute sparse cross covariance.
#' @param x Matrix.
#' @param y Matrix.
#' @return List with cov.
#' @noRd
sparse.cov.cross <- function(x, y) {
  n <- nrow(x)
  cMeans.x <- colMeans(x); cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x, cMeans.y))/(n-1)
  list(cov=covmat)
}
