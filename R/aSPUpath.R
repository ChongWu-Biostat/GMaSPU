#' Pathway based Sum of Powered Score tests (SPUsPath) and adaptive SPUpath (aSPUsPath) test for single trait - pathway association starting with score vector and its covariance matrix.
#'
#' It gives p-values of the SPUsPath tests and aSPUsPath test starting with score vector and its covariance matrix.
#'
#' @param U Score vector for the marker set we are interested in.
#'
#' @param V Corresponding covariance matrix for the score vector.
#'
#' @param info A table for mapping each marker to the specific genes and chromosome.
#'
#' @param pow SNP or single CpG sites specific power(gamma values) used in SPUpath test.
#'
#' @param pow2 GENE specific power(gamma values) used in SPUpath test.
#'
#' @param n.perm number of permutations.
#'
#' @export
#' @return P-values for SPUMpath tests and aSPUMpath test.
#'
#' @author Chong Wu and Wei Pan
#'
#' @references
#' Il-Youp Kwak, Wei Pan (2015)
#' Adaptive Gene- and Pathway-Trait Association Testing with GWAS Summary Statistics,
#' Bioinformatics, 32(8):1178-1184
#'
aSPUPath <- function (U, V, info, pow = c(1, 2, 4, 8, Inf), pow2 = c(1, 2, 4, 8),
 n.perm = 1000) {
    
    info = methy.info
    
    chrs <- unique(methy.info[,"chrom.id"])
    CH <- list(0)
    CH.CovSsqrt <- list(0)
    for (i in 1:length(chrs)) {
        c = chrs[i]
        CH[[i]] = unname(U[methy.info[methy.info[,"chrom.id"]==c,1],1])
        
        Covtemp = unname(V[methy.info[methy.info[,"chrom.id"]==c,1],methy.info[methy.info[,"chrom.id"]==c,1]])
        eV <- eigen(Covtemp, symmetric = TRUE)
        covsqrt = t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
        rownames(covsqrt) = rownames(Covtemp)
        colnames(covsqrt) = colnames(Covtemp)
        CH.CovSsqrt[[i]] <- covsqrt
    }
    
    nSNPs0 = table(methy.info[,"gene.id"])
    gene.index = unique(methy.info[,"gene.id"])
    nGenes = length(gene.index)
    StdTs <- rep(0, length(pow) * nGenes)
    for (j in 1:length(pow)) {
        for (iGene in 1:nGenes) {
            
            indx = methy.info[methy.info[,"gene.id"] == iGene,1]
            if (pow[j] < Inf) {
                a = (sum(U[indx,1]^pow[j]))
                StdTs[(j - 1) * nGenes + iGene] = sign(a) * ((abs(a)/nSNPs0[iGene])^(1/pow[j]))
            } else {
                StdTs[(j - 1) * nGenes + iGene] = max(abs(U[indx,1]))
            }
        }
    }
    
    Ts2 <- rep(0, length(pow) * length(pow2))
    for (j2 in 1:length(pow2)) {
        for (j in 1:length(pow)) {
            if (pow2[j2] < Inf) {
                Ts2[(j2 - 1) * length(pow) + j] = sum(StdTs[((j -
                1) * nGenes + 1):(j * nGenes)]^pow2[j2])
            }
            else {
                Ts2[(j2 - 1) * length(pow) + j] = max(StdTs[((j -
                1) * nGenes + 1):(j * nGenes)])
            }
        }
    }
    
    if(max(pow) == Inf) {
        pow[which(pow ==Inf)] = -1
    }
    
    if(max(pow2) == Inf) {
        pow2[which(pow2 ==Inf)] = -1
    }
    
    s <- sample(1:10^7, 1)
    

    k = dim(U)[1]
    nSNPs0 =  as.vector(nSNPs0)
    Ts2 = as.vector(Ts2)
    
    nChrom0 = as.vector(unname(table(methy.info[,"chrom.id"])))

    Results = aSPUsPathEngine2(CH, CH.CovSsqrt, pow, pow2, nGenes, n.perm , k, nSNPs0,nChrom0, Ts2, s)
    minp0 <- Results$minp0
    pPerm0 <- Results$pPerm0
    P0s <- Results$P0
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    Ts <- c(Ts2, min(pPerm0))
    
    if(min(pow) == -1) {
        pow[which(pow == -1 )] = Inf
    }
    
    if(min(pow2) == -1) {
        pow2[which(pow2 == -1)] = Inf
    }
    
    nmvec <- NULL
    for (ii in pow2) {
        for (jj in pow) {
            nmvec <- c(nmvec, paste("SPUsPath", jj, ",", ii,
            sep = ""))
        }
    }
    nmvec <- c(nmvec, "aSPUsPath")
    names(Ts) = nmvec
    names(pvs) = nmvec
    list(Ts = Ts, pvs = pvs)
}
