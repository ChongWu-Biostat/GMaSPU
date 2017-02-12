#' Sum of Powered Score (SPU) tests and adaptive SPU (aSPU) test for single trait - SNP set association.
#'
#' It gives p-values of the SPU tests and aSPU test.
#'
#' @param U Score vector for the marker set we are interested in.
#'
#' @param V Corresponding covariance matrix for the score vector.
#'
#' @param pow power used in SPU test. A vector of the powers.
#'
#' @param n.perm number of permutations or bootstraps.
#'
#' @return A list object, Ts : test statistics for the SPU tests (in the order of the specified pow) and finally for the aSPU test.
#'         pvs : p-values for the SPU and aSPU tests.
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014)
#' A powerful and adaptive association test for rare variants,
#' Genetics, 197(4), 1081-95
#'

aSPU <- function(U,V,pow=c(1:8,Inf), n.perm = 1000){
    
    Ts <- rep(0, length(pow))
    for(j in 1:length(pow)){
        if (pow[j] < Inf)
        Ts[j] = sum(U^pow[j]) else Ts[j] = max(abs(U))
    }
    
    eV<-eigen(V)
    CovSsqrt<- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] = 0 # pass 0 as infitiy
    T0sC <- calcT0sim(as.matrix(CovSsqrt), as.matrix(pow), n.perm)
    T0s <- T0sC$T0
    pPerm0 = rep(NA,length(pow))
    
    for (j in 1:length(pow)) {
        pPerm0[j] = sum(abs(Ts[j])<=abs(T0s[,j])) / n.perm
        P0s = ( ( n.perm - rank( abs(T0s[,j]) ) ) + 1 ) / (n.perm )
        if (j == 1 ) minp0  = P0s else minp0[which(minp0>P0s)] = P0s[which(minp0>P0s)]
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm+1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    names(Ts) <- c(paste("SPU", pow, sep=""), "aSPU")
    names(pvs) = names(Ts)
    
    list(Ts = Ts, pvs = pvs)
    
}
