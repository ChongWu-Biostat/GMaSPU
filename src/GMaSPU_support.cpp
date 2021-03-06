#include <armadillo>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List calcT0sim (arma::mat CvSqrt, arma::mat powV, int nperm) {
    //const int n = X1.n_rows;
    const int npow = powV.n_rows;
    const int k = CvSqrt.n_rows;
    // containers
    arma::mat T0s1(nperm,npow);
    T0s1.fill(0);
    
    for (int i = 0; i < nperm; i++) {
        arma::mat U00 = arma::randn(k,1);
        arma::mat U0 = CvSqrt * U00;
        
        for (int j = 0; j < npow; j++) {
            if (powV(j,0) == 0) {
                arma::mat tmpU01 = abs(U0);
                T0s1(i,j) = tmpU01.max();
            } else {
                T0s1(i,j) = accu(pow(U0,powV(j,0)));
            }
        }
    }
    Rcpp::List res;
    res["T0"] =T0s1;
    
    return(res);
}


// [[Rcpp::export]]
void set_seed(unsigned int seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);
}

// This function is taken from http://stackoverflow.com/questions/39153082/rcpp-rank-function-that-does-average-ties
class Comparator {
    private:
    const Rcpp::NumericVector& ref;
    
    bool is_na(double x) const
    {
        return Rcpp::traits::is_na<REALSXP>(x);
    }
    
    public:
    Comparator(const Rcpp::NumericVector& ref_)
    : ref(ref_)
    {}
    
    bool operator()(const int ilhs, const int irhs) const
    {
        double lhs = ref[ilhs], rhs = ref[irhs];
        if (is_na(lhs)) return false;
        if (is_na(rhs)) return true;
        return lhs < rhs;
    }
};


// [[Rcpp::export]]
Rcpp::NumericVector avg_rank(Rcpp::NumericVector x)
{
    R_xlen_t sz = x.size();
    Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
    std::sort(w.begin(), w.end(), Comparator(x));
    
    Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
    for (R_xlen_t n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
        for (R_xlen_t k = 0; k < n; k++) {
            r[w[i + k]] = i + (n + 1) / 2.;
        }
    }
    
    return r;
}





// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List aSPUsPathEngine2(Rcpp::List CH, Rcpp::List CHcovSq, arma::vec pow1, arma::vec pow2, int nGenes, int n_perm, int k, arma::vec nSNPs0, arma::vec nChrom0, arma::vec Ts2, int s) {
    const int n_ch = CH.size();
    const int n_pow1 = pow1.size();
    const int n_pow2 = pow2.size();
    
    arma::vec T0s(n_perm);
    arma::vec T0st(nGenes*n_pow1);
    arma::vec Ts2t(n_pow1*n_pow2);
    arma::vec pPerm0(n_pow1*n_pow2);
    arma::vec minp0(n_perm);
    arma::vec P0s(n_perm);
    
    // iterate for pow2
    for(int j2=0 ; j2 < n_pow2; j2++) {
        // iterate for pow1
        for(int j=0 ; j < n_pow1; j++) {
            
            // set seed to use same random numbers for same b's
            // This is necessary to use efficient memory
            set_seed(s);
            for(int b=0; b < n_perm; b++) {
                
                // Generate Score from null distribution
                Rcpp::NumericVector U00 = rnorm(k,0,1);
                arma::vec u0 = as<arma::vec>(U00);
                arma::vec U0(1);
                
                int SNPstart = 0;
                // iterate for chromosome
                for(int b2 = 0; b2 < n_ch ; b2++) {
                    
                    
                    // itreate for chromosome
                    if( b2 != 0) {
                        SNPstart = sum( nChrom0.subvec(0,b2-1) ) ;
                    }
                    int idx1 = SNPstart;
                    int idx2 = SNPstart+nChrom0(b2)-1;
                    
                    arma::vec TT = as<arma::mat>(CHcovSq[b2])*u0.subvec(idx1, idx2) ;
                    U0 = join_cols(U0,TT) ;
                }
                
                arma::vec UU2 = U0.subvec(1,U0.size()-1);
                
                
                // iterate for genes
                SNPstart = 0;
                for(int iGene = 0 ; iGene < nGenes ; iGene++ ) {
                    
                    // calculate starting and ending position of each gene from nSNPs0 vector
                    if( iGene != 0) {
                        SNPstart = sum( nSNPs0.subvec(0,iGene-1) ) ;
                    }
                    int idx1 = SNPstart;
                    int idx2 = SNPstart+nSNPs0(iGene)-1;
                    
                    // calculate 1st level test statistic
                    if( pow1[j] > 0) {
                        arma::vec tmp1 = pow(UU2.subvec(idx1, idx2),pow1[j]);
                        double tmp2 = sum(tmp1);
                        
                        if( tmp2 > 0 ) {
                            T0st(j*nGenes+iGene) = pow(std::abs(tmp2)/nSNPs0[iGene] , 1/pow1[j]);
                        } else {
                            T0st(j*nGenes+iGene) = -pow(std::abs(tmp2)/nSNPs0[iGene] , 1/pow1[j]);
                        }
                        
                    } else {
                        arma::vec T0tp = abs(UU2.subvec(idx1, idx2));
                        T0st(j*nGenes+iGene) = max(abs(T0tp));
                    }
                }
                
                // calculate 2nd level test statistics
                if( pow2[j2] > 0) {
                    arma::vec tmp3 = pow(T0st.subvec(j*nGenes,(j+1)*nGenes-1), pow2[j2]);
                    double tmp4 = sum(tmp3);
                    T0s(b) = tmp4;
                } else {
                    T0s(b) = max( arma::abs(T0st.subvec(j*nGenes,(j+1)*nGenes-1)) ) ;
                }
            }
            
            
            // Calculate P-values
            int tmp3 = 0;
            arma::vec T0sabs = arma::abs(T0s);
            Rcpp::NumericVector a( T0sabs.begin(), T0sabs.end() );
            Rcpp::NumericVector ranka = avg_rank(a);
            arma::vec rankarma = as<arma::vec>(ranka);
            
            for( int tt=0 ; tt < n_perm ; tt++) {
                if( std::abs(Ts2(j2*n_pow1 + j) ) <= std::abs(T0s(tt))) {
                    tmp3++;
                }
                
                P0s(tt) = (double) (n_perm - rankarma(tt) + 1) / (double) n_perm;
            }
            
            minp0 = P0s;
            
            if(j == 1) {
                for( int ii=0; ii < n_perm; ii++) {
                    minp0(ii) = P0s(ii);
                }
            } else {
                for( int ii=0; ii < n_perm; ii++) {
                    if( minp0(ii) > P0s(ii) ) {
                        minp0(ii) = P0s(ii);
                    }
                }
            }
            
            pPerm0(j2*n_pow1 + j) = (double) tmp3 / (double) n_perm;
            
        }
    }
    return Rcpp::List::create(Rcpp::Named("minp0") = minp0,
                              Rcpp::Named("pPerm0") = pPerm0,
                              Rcpp::Named("P0s") = P0s);
}

