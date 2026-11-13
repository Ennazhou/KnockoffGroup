#' Knockoff Group Analysis for Trio Design
#'
#' This file contains functions for performing knockoff group analysis in trio designs.
#' The main function `KnockoffGroup_trio` computes feature statistics for trio data.
#' Auxiliary functions include FBAT calculations, ACAT, and more.
#'
#' @name knockoff_trio_functions
NULL

#' Calculate KnockoffTrio's Feature Statistics
#'
#' Compute statistics using original and knockoff genotype data for trios.
#'
#' @param geneinfo Gene information matrix.
#' @param dat Original genotype matrix.
#' @param dat.ko Knockoff genotype array.
#' @param pos Variant positions.
#' @param start Start position (optional).
#' @param end End position (optional).
#' @param size Window sizes retained for compatibility with the runtime scripts.
#' @param p_value_only If TRUE, only compute p-values.
#' @param adjust_for_cov If TRUE, adjust for covariates.
#' @param y Residual vector.
#' @param chr Chromosome.
#' @param xchr If TRUE, X chromosome analysis.
#' @param sex Sex vector for offspring.
#' @return Data frame with results.
#' @importFrom stats as.dist cutree hclust median pcauchy pnorm princomp dbeta na.omit
#' @importFrom CompQuadForm davies
#' @export
KnockoffGroup_trio <- function(geneinfo, dat, dat.ko=NA, pos, start=NA, end=NA,
                               size=c(1,1000,5000,10000,20000,50000),
                               p_value_only=FALSE, adjust_for_cov=FALSE,
                               y=NA, chr="1", xchr=FALSE, sex=NA){
  if (nrow(geneinfo) > 1L) {
    named <- all(c("start", "end", "nsnp", "gene_id") %in% colnames(geneinfo))
    results <- lapply(seq_len(nrow(geneinfo)), function(i) {
      gene <- geneinfo[i, , drop = FALSE]
      left <- if (named) gene$start else gene[1, 5]
      right <- if (named) gene$end else gene[1, 6]
      columns <- which(pos >= left & pos <= right)
      if (length(columns) == 0L) return(NULL)
      if (named) gene$nsnp <- length(columns) else gene[1, 4] <- length(columns)
      KnockoffGroup_trio(gene, dat[, columns, drop = FALSE],
        if (p_value_only) NA else dat.ko[, columns, , drop = FALSE],
        pos[columns], size = size, p_value_only = p_value_only,
        adjust_for_cov = adjust_for_cov, y = y, chr = chr, xchr = xchr, sex = sex)
    })
    return(do.call(rbind, results))
  }
  if (nrow(dat) %% 3!=0) stop("The number of rows of the original trio matrix must be a multiple of three.")
  if (xchr & is.na(sex)) {
    stop("Gender information is required if KnockoffGroup_trio is applied to the X chromosome")
  }
  #if (!all(dat %in% c(0,1,2))) stop("The original trio matrix can only contain 0, 1, or 2.")
  if (is.na(start)) start<-min(pos)
  if (is.na(end)) end<-max(pos)
  n<-nrow(dat)/3 #number of trios
  nsnp<-ncol(dat)
  out_fbat<-fbat(dat,adjust_for_cov,y,xchr=xchr,sex=sex)
  z<-out_fbat$additive
  p1<-2*pnorm(-abs(z)) #p-values for SNP-based FBAT in original trios
  p1_sign<-sign(z)
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  mac<-apply(dat[-index_off,,drop=FALSE],2,sum) #minor allele count for parents
  if (xchr==FALSE) maf<-mac/(4*n) else maf<-(apply(dat[index_dad,,drop=FALSE],2,sum)/n+apply(dat[index_mom,,drop=FALSE],2,sum)/(2*n))/2
  geno<-dat[index_off,,drop=FALSE]-(dat[index_dad,,drop=FALSE]+dat[index_mom,,drop=FALSE])/2
  Y=rep(1,n)
  Ymean <- sum(Y)
  #
  #
  #
  #Prepare Z matrix for set-based fbat
  weight<-1/(sqrt(n*maf*(1-maf)))
  weight[which(weight==Inf)]<-0
  #original
  Z<-sweep(out_fbat$Z, 2, weight, '*')
  Zsq<-t(Z)%*%Z ###var
  Zsum<-apply(Z,2,sum)
  #prepare windows gene based
  window<-data.frame()
  if (all(c("start", "end", "nsnp", "gene_id") %in% colnames(geneinfo))) {
    window<-data.frame(start=geneinfo$start,end=geneinfo$end,ind3=1,ind4=geneinfo$nsnp,n=geneinfo$nsnp, gene_no=geneinfo$gene_id)
  } else if (ncol(geneinfo) >= 6) {
    window<-data.frame(start=geneinfo[,5],end=geneinfo[,6],ind3=1,ind4=geneinfo[,4],n=geneinfo[,4], gene_no=geneinfo[,3])
  } else {
    stop("geneinfo must contain columns start, end, nsnp, and gene_id, or at least six runtime columns.")
  }
  window<-window[which(window$n>0),]
  window<-window[!duplicated(window[,c(3,4)]),]
  nwindow<-nrow(window)
  q1<-p.fbat<-z1<-rep(NA,nwindow) #window p-values for original trios q1.cv<-q1.rv<-q1.pa<-q1.ma<-
  #q1.burden.cv<-q1.burden.rv<-q1.burden.urv<-rep(NA,nwindow) #p.denovo
  #n.denovo<-rep(0,nwindow)
  direction<-rep(NA,nwindow) #direction.pa<-direction.ma<-
  if (!p_value_only){
    #The following three warnings only work when xchr==FALSE
    #if (dim(dat.ko)[1] %% 3!=0) stop("The number of rows of the knockoff trio matrix must be a multiple of three.")
    #if (nrow(dat)!=dim(dat.ko)[1]) stop("The number of rows of the two matrices must be equal.")
    #if (ncol(dat)!=dim(dat.ko)[2]) stop("The number of columns of the two matrices must be equal.")
    M<-dim(dat.ko)[3]
    out_fbat<-fbat(dat,adjust_for_cov,y,dosage=TRUE,dat.ko,xchr,sex) #SNP-based FBAT statistics for knockoff trios
    z.ko<-out_fbat$additive
    Z2<-out_fbat$Z
    p2<-2*pnorm(-abs(z.ko)) #p-values for SNP-based FBAT in knockoff trios, dim=c(nsnp,M)
    #
    Z2sq<-array(dim=c(nsnp,nsnp,M))
    Z2sum<-array(dim=c(nsnp,M))
    if (dim(Z2)[2]>1){
      for (m in 1:M){
        Z2[,,m]<-sweep(Z2[,,m], 2, weight, '*')
        Z2sq[,,m]<-t(Z2[,,m])%*%Z2[,,m]
        Z2sum[,m]<-apply(Z2[,,m],2,sum)
      }
    } else if (dim(Z2)[2]==1){
      for (m in 1:M){
        Z2[,,m]<-Z2[,,m]*weight
        Z2sq[,,m]<-t(Z2[,,m])%*%Z2[,,m]
        Z2sum[,m]<-sum(Z2[,,m])
      }
    }
    q2<-z2<-array(dim=c(M,nwindow))
    p.burden.ko<-array(dim=c(M,nwindow))
  }
  for (i in 1:nwindow) {
    ind<-window$ind3[i]:window$ind4[i]
    if (length(ind)==1){
      q1[i]<-p.fbat[i]<-p1[ind] #q1.cv[i]
      z1[i]<-z[ind]
      direction[i]<-p1_sign[ind]
    }else{
      z1[i]<-fbat_set(Zsum,Zsq,ind)
      p.fbat[i]<-2*pnorm(-abs(z1[i]))
      #single-variant for all
      ind.single<-ind #[mac[ind]>=5]
      p.single<-p1[ind.single]
      q1[i]<-ACAT(c(p.single,p.fbat[i]))### p.weight.ACAT2
      #direction of q1[i]
      stat<-c(z1[i],z[ind.single])
      tmp<-which.max(abs(stat))
      if (length(tmp)>0) direction[i]<-sign(stat[tmp])
    }
    if (!p_value_only){
      if (length(ind)==1){
        q2[,i]<-p2[ind,] #q2.cv[,i]<-
        z2[,i]<-z.ko[ind,]
        p.burden.ko[,i]<-2*pnorm(-abs(z2[,i]))
      }else{
        for (m in 1:M){
          p.single.ko<-p2[ind.single,m]
          z2[m,i]<-fbat_set(Z2sum[,m],Z2sq[,,m],ind)
          p.burden.ko.temp<-2*pnorm(-abs(z2[m,i]))
          q2[m,i]<-ACAT(c(p.single.ko,p.burden.ko.temp))
        }
      }
    }
  }
  window$ind3<-pos[window$ind3]
  window$ind4<-pos[window$ind4]
  colnames(window)[3]<-"actual_start"
  colnames(window)[4]<-"actual_end"
  if (!p_value_only){
    out<-calculate_w_kappatau(q1,q2)
    w<-out$w
    kappatau<-out$kappatau
    rownames(q2)<-paste0("p_",1:M)
    q2<-t(q2)
    rownames(z2)<-paste0("z_",1:M)
    z2<-t(z2)
    rownames(p.burden.ko)<-paste0("p.burden_",1:M)
    p.burden.ko<-t(p.burden.ko)
    #p: acat p, z: z score for single-variant or set FBAT (depends on if a window contains 1 or more variants), p.burden: p value for single-variant or set FBAT
    window<-cbind(chr,window,dir=direction,w,p=q1,z=z1,p.burden=p.fbat,kappatau,q2,z2)
  } else window<-cbind(chr,window,dir=direction,p.WA=q1,z=z1,p.burden=p.fbat)
  return(window)
}

#' FBAT Calculation
#'
#' Compute FBAT statistics.
#'
#' @param dat Genotype matrix.
#' @param adjust_for_covariates Adjust for covariates.
#' @param y Residuals.
#' @param dosage Dosage mode.
#' @param dat1 Dosage data.
#' @param xchr X chromosome.
#' @param sex Sex vector.
#' @return List with statistics.
#' @noRd
fbat<-function(dat,adjust_for_covariates=FALSE,y=NA,dosage=FALSE,dat1=NA,xchr=FALSE,sex=NA){
  #if (xchr & is.na(sex)) stop("Gender information is required if FBAT is applied to the X chromosome")
  if (is.null(ncol(dat))) dat<-as.matrix(dat)
  n<-nrow(dat)/3
  index_dad<-seq(1,(3*n),3) #index for parent1
  index_mom<-seq(2,(3*n),3) #index for parent2
  index_off<-seq(3,(3*n),3) #index for offspring
  if (!xchr){
    # difference of genotype between offspring and the average of parents
    if (dosage==FALSE) Z<-dat[index_off,,drop=FALSE]-(dat[index_dad,,drop=FALSE]+dat[index_mom,,drop=FALSE])/2 else Z<-dat1[index_off,,,drop=FALSE]-(dat1[index_dad,,,drop=FALSE]+dat1[index_mom,,,drop=FALSE])/2
  } else Z<-xcontribution(dat,dat1,dosage,sex)
  tat<-pcontribution(dat,dat1,dosage,xchr,sex)
  Z.pa<-tat$Z.pa
  Z.ma<-tat$Z.ma
  if (adjust_for_covariates){
    Z<-sweep(Z,1,y,'*')
    Z.pa<-sweep(Z.pa,1,y,'*')
    Z.ma<-sweep(Z.ma,1,y,'*')
  }
  additive<-colSums(Z)/sqrt(colSums(Z^2)) #if dosage=T, dim=c(nsnp,M) ###variance
  paternal<-colSums(Z.pa)/sqrt(colSums(Z.pa^2)) #if dosage=T, dim=c(nsnp,M)
  maternal<-colSums(Z.ma)/sqrt(colSums(Z.ma^2)) #if dosage=T, dim=c(nsnp,M)
  additive[is.na(additive)]<-0
  paternal[is.na(paternal)]<-0
  maternal[is.na(maternal)]<-0
  out<-list()
  out$additive<-additive
  out$paternal<-paternal
  out$maternal<-maternal
  out$Z<-Z
  out$Z.pa<-Z.pa
  out$Z.ma<-Z.ma
  return (out)
}

#' Set-Based FBAT
#'
#' Compute set-based FBAT.
#'
#' @param W Sum vector.
#' @param V Variance matrix.
#' @param ind Indices.
#' @return Statistic.
#' @noRd
fbat_set<-function(W,V,ind){
  if (length(ind)>0){
    s1<-sum(W[ind])
    s2<-sum(V[ind,ind]) ###V: var-cov
    if (s2==0) s2<-1
    return (unname(s1/sqrt(s2)))
  } else return(NA)
}

#' Parent Contribution
#'
#' Compute parent contributions.
#'
#' @param dat Genotype.
#' @param dat1 Dosage data.
#' @param dosage Dosage mode.
#' @param xchr X chromosome.
#' @param sex Sex vector.
#' @return List with Z.pa and Z.ma.
#' @noRd
pcontribution <- function(dat,dat1=NA,dosage=FALSE,xchr=FALSE,sex=NA){ #dat: genotype matrix, dat1: related to dosage
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]
  if (dosage==FALSE){
    Z.pa<-Z.ma<-array(0,dim=c(n.row/3,ncol(dat)))
    het <- (mom == 1) # mom geno == 1
    hethom <- het & (dad == 2) # mom geno == 1 & dad geno == 2
    # ind212<-hethom & (kid == 2)
    # ind211<-hethom & (kid == 1)
    # n212 <- colSums(ind212, na.rm=TRUE)
    # n211 <- colSums(ind211, na.rm=TRUE)
    Z.ma[hethom & (kid == 2)]<-1 # dad, mom, off == 2,1,2; Z.ma = 1
    Z.ma[hethom & (kid == 1)]<-(-1) # dad, mom, off == 2,1,1; Z.ma = -1
    hethom <- het & (dad == 0) # mom geno == 1 & dad geno == 0
    # ind011<-hethom & (kid == 1)
    # ind010<-hethom & (kid == 0)
    # n011 <- colSums(ind011, na.rm=TRUE)
    # n010 <- colSums(ind010, na.rm=TRUE)
    Z.ma[hethom & (kid == 1)]<-1 # dad, mom, off == 2,1,1; Z.ma = 1 # Entries set to be 1 or -1?
    Z.ma[hethom & (kid == 0)]<-(-1) # dad, mom, off == 2,1,0; Z.ma = -1
    het <- (dad == 1)
    hethom <- het & (mom == 2)
    # ind122<-hethom & (kid == 2)
    # ind121<-hethom & (kid == 1)
    # n122 <- colSums(ind122, na.rm=TRUE)
    # n121 <- colSums(ind121, na.rm=TRUE)
    Z.pa[hethom & (kid == 2)]<-1
    Z.pa[hethom & (kid == 1)]<-(-1)
    hethom <- het & (mom == 0)
    # ind101<-hethom & (kid == 1)
    # ind100<-hethom & (kid == 0)
    # n101 <- colSums(ind101, na.rm=TRUE)
    # n100 <- colSums(ind100, na.rm=TRUE)
    Z.pa[hethom & (kid == 1)]<-1
    Z.pa[hethom & (kid == 0)]<-(-1)
    het <- (mom == 1) & (dad == 1)
    ind112<-het & (kid == 2)
    ind110<-het & (kid == 0)
    # n112 <- colSums(ind112, na.rm=TRUE)
    # n110 <- colSums(ind110, na.rm=TRUE)
    Z.pa[ind112]<-Z.ma[ind112]<-1
    Z.pa[ind110]<-Z.ma[ind110]<-(-1)
    #mat<-cbind(n100,n110,n121,n101,n112,n122,n010,n110,n211,n011,n112,n212)
    #colnames(mat)<-c("p-1","p-1","p-1","p+1","p+1","p+1","m-1","m-1","m-1","m+1","m+1","m+1") #contributions
    if (xchr){
      ind <- (dad == 1) & (mom == 1) & (kid==1)
      ind0<-sweep(ind,1,sex==0,'&')
      Z.pa[ind0]<-1
      Z.ma[ind0]<-(-1)
      ind0<-sweep(ind,1,sex==1,'&')
      Z.pa[ind0]<-(-1)
      Z.ma[ind0]<-1
    }
    out<-list()
    out$Z.pa<-Z.pa
    out$Z.ma<-Z.ma
    return(out)
  }else{
    dims<-dim(dat1)
    dims[1]<-dims[1]/3
    Z.pa<-Z.ma<-array(0,dim=dims)
    dadk <- dat1[seq.int(1, n.row, 3),,,drop=FALSE]
    momk <- dat1[seq.int(2, n.row, 3),,,drop=FALSE]
    kidk <- dat1[seq.int(3, n.row, 3),,,drop=FALSE]
    het <- (mom == 1)
    hethom <- het & (dad == 2)
    ind212<-hethom & (kid == 2)
    ind211<-hethom & (kid == 1)
    # n212 <- colSums(ind212, na.rm=TRUE)
    # n211 <- colSums(ind211, na.rm=TRUE)
    Z.pa[ind212]<-kidk[ind212]-dadk[ind212]
    Z.ma[ind212]<-kidk[ind212]-momk[ind212]
    Z.pa[ind211]<-kidk[ind211]-momk[ind211]
    Z.ma[ind211]<-kidk[ind211]-dadk[ind211]
    hethom <- het & (dad == 0)
    ind011<-hethom & (kid == 1)
    ind010<-hethom & (kid == 0)
    # n011 <- colSums(ind011, na.rm=TRUE)
    # n010 <- colSums(ind010, na.rm=TRUE)
    Z.pa[ind011]<-kidk[ind011]-momk[ind011]
    Z.ma[ind011]<-kidk[ind011]-dadk[ind011]
    Z.pa[ind010]<-kidk[ind010]-dadk[ind010]
    Z.ma[ind010]<-kidk[ind010]-momk[ind010]
    het <- (dad == 1)
    hethom <- het & (mom == 2)
    ind122<-hethom & (kid == 2)
    ind121<-hethom & (kid == 1)
    # n122 <- colSums(ind122, na.rm=TRUE)
    # n121 <- colSums(ind121, na.rm=TRUE)
    Z.pa[ind122]<-kidk[ind122]-dadk[ind122]
    Z.ma[ind122]<-kidk[ind122]-momk[ind122]
    Z.pa[ind121]<-kidk[ind121]-momk[ind121]
    Z.ma[ind121]<-kidk[ind121]-dadk[ind121]
    hethom <- het & (mom == 0)
    ind101<-hethom & (kid == 1)
    ind100<-hethom & (kid == 0)
    # n101 <- colSums(ind101, na.rm=TRUE)
    # n100 <- colSums(ind100, na.rm=TRUE)
    Z.pa[ind101]<-kidk[ind101]-momk[ind101]
    Z.ma[ind101]<-kidk[ind101]-dadk[ind101]
    Z.pa[ind100]<-kidk[ind100]-dadk[ind100]
    Z.ma[ind100]<-kidk[ind100]-momk[ind100]
    het <- (mom == 1) & (dad == 1)
    ind112<-het & (kid == 2)
    ind110<-het & (kid == 0)
    # n112 <- colSums(ind112, na.rm=TRUE)
    # n110 <- colSums(ind110, na.rm=TRUE)
    Z.pa[ind112]<-kidk[ind112]-momk[ind112]
    Z.ma[ind112]<-kidk[ind112]-dadk[ind112]
    Z.pa[ind110]<-kidk[ind110]-dadk[ind110]
    Z.ma[ind110]<-kidk[ind110]-momk[ind110]
    if (xchr){
      ind <- (dad == 1) & (mom == 1) & (kid==1)
      ind0<-sweep(ind,1,sex==0,'&')
      Z.pa[ind0]<-dadk[ind0]
      Z.ma[ind0]<-(-1)*momk[ind0]
      ind0<-sweep(ind,1,sex==1,'&')
      Z.pa[ind0]<-(-1)*dadk[ind0]
      Z.ma[ind0]<-momk[ind0]
    }
    #mat<-cbind(n100,n110,n121,n101,n112,n122,n010,n110,n211,n011,n112,n212)
    #colnames(mat)<-c("p-1","p-1","p-1","p+1","p+1","p+1","m-1","m-1","m-1","m+1","m+1","m+1") #contributions
    out<-list()
    out$Z.pa<-Z.pa
    out$Z.ma<-Z.ma
    return(out)
  }
}

#' X Contribution
#'
#' Compute X chromosome contributions.
#'
#' @param dat Genotype.
#' @param dat1 Dosage data.
#' @param dosage Dosage mode.
#' @param sex Sex vector.
#' @return Z matrix.
#' @noRd
xcontribution <- function(dat,dat1=NA,dosage=FALSE,sex=NA){
  #sex: gender info for the offspring, 0 for females and 1 for males
  #Only knockoffs of mothers are required
  n.row <- nrow(dat)
  dad <- dat[seq.int(1, n.row, 3),, drop=FALSE]
  mom <- dat[seq.int(2, n.row, 3),, drop=FALSE]
  kid <- dat[seq.int(3, n.row, 3),, drop=FALSE]
  if (dosage==FALSE){
    Z<-array(0,dim=c(n.row/3,ncol(dat)))
    het <- (mom == 1)
    hethom <- het & (dad == 0)
    # dad geno == 0, mom geno == 1
    Z[hethom & (kid == 0)]<-(-0.5) # 010
    Z[hethom & (kid == 1)]<-0.5 # 011
    hethom <- het & (dad == 1)
    # dad geno == 1, mom geno == 1
    ind111<-hethom & (kid == 1) # 111
    Z[hethom & (kid == 0)]<-(-0.5) #110 trios on chrX can only have male offspring
    Z[sweep(ind111,1,sex==0,'&')]<-(-0.5)
    # sweep(ind111,1,sex==0,'&'): TRUE -> NA
    Z[sweep(ind111,1,sex==1,'&')]<-0.5
    # sweep(ind111,1,sex==1,'&'): TRUE -> NA. It is the same index as the above one?
    Z[hethom & (kid == 2)]<-0.5 #112 trios on chrX can only have female offspring
  }else{
    dims<-dim(dat1)
    dims[1]<-dims[1]/3
    Z<-array(0,dim=dims)
    dadk <- dat1[seq.int(1, n.row, 3),,,drop=FALSE]
    momk <- dat1[seq.int(2, n.row, 3),,,drop=FALSE]
    kidk <- dat1[seq.int(3, n.row, 3),,,drop=FALSE]
    het<-dad==0 & mom==1
    Z[het]<-kidk[het]-0.5*(dadk[het]+momk[het])
    het<-dad==1 & mom==1
    het0<-sweep(het,1,sex==0,'&')
    Z[het0]<-kidk[het0]-0.75*(dadk[het0]+momk[het0])
    het0<-sweep(het,1,sex==1,'&')
    Z[het0]<-kidk[het0]-0.25*(dadk[het0]+momk[het0])
  }
  return(Z)
}

#' Calculate W, Kappa, Tau
#'
#' Compute W, kappa, tau from p-values.
#'
#' @param q1 Original p-values.
#' @param q2 Knockoff p-values.
#' @return List with w and kappatau.
#' @noRd
calculate_w_kappatau<-function(q1,q2){
  out<-list()
  t1<--log10(q1)
  t2<--log10(q2)
  t2_med<-apply(t2,2,median)
  t2_max<-apply(t2,2,max)
  out$w<-(t1-t2_med)*(t1>=t2_max)
  out$w.raw<-t1-t2_med
  out$kappatau<-MK.statistic(t1,t(t2),method="median")
  #out$q<-MK.q.byStat(out$kappatau[,1],out$kappatau[,2],M=M)
  return(out)
}

#' MK Statistic
#'
#' Compute MK statistic.
#'
#' @param T_0 Original.
#' @param T_k Knockoff.
#' @param method Method.
#' @return Matrix with kappa and tau.
#' @noRd
MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  kappa<-apply(T.temp,1,which.max)-1
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}

#' ACAT Combination
#'
#' ACAT p-value combination.
#'
#' @param p P-values.
#' @return Combined p-value.
#' @noRd
ACAT<-function(p){
  p[p>0.99]<-0.99 #p[p>1-1e-2]<-1-1e-2
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))
  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}
