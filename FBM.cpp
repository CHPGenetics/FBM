#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <stdlib.h>
// [[Rcpp::depends(BH)]]
//#include <boost/math/distributions/normal.hpp>
//#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/math/distributions.hpp>
#include <boost/regex.hpp>
#include <boost/random/variate_generator.hpp>
#include <algorithm>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
//So...alright no comment is allowed between export line and the function

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

double sigma0 = 1000;
double tau0 = 0.05;
double tau0Omega= 0.1;
double tauAlpha=0.1;
double tauBeta=0.1; // #once set differently, it'll be informative prior on tau
double smallNumber=0.000001;
double sigmam0 = 1;
double ae = 200; // Anti explosion parameter
int itrCheckConvergence=10000;
int itrStart=800; //if maxItr>itrCheckConvergence, from 800 start to check convergence.
boost::random::mt19937 rng;

// [[Rcpp::export]]
double rInvGamma(const double shape, const double scale){
    boost::gamma_distribution<> gd( shape );
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > var_gamma( rng, gd );
    double out = 0;
    out = scale/(var_gamma());
    //we don't use randg because it doesn't take in double shape and scale.
    return out;
}


// [[Rcpp::export]]
double dNorm(const double x, const double mean, const double sd){
    boost::math::normal_distribution<> my_normal(mean,sd);
    double out = boost::math::pdf(my_normal,x);
    return out;
}

// [[Rcpp::export]]
double rBeta(const double shape1, const double shape2){
    double randFromUnif = randu();
    //parameters and the random value on (0,1) you drew
    boost::math::beta_distribution<> dist(shape1, shape2);
    double randFromDist = quantile(dist, randFromUnif);
    return randFromDist;
}

// [[Rcpp::export]]
bool naiveGeweke(const vec a, const vec b){
    bool converged=false;
    double geweke = abs(arma::median(a) - arma::median(b))/(std::sqrt(arma::var(a)/20.0 + arma::var(b)/100.0));
    if (geweke < 1.96){
        converged=true;
    }
    return converged;
}


class parmUpdate{
public:
    int nItr;
    int nItrActual;
    int nBurninActual;
    vec Y; // N x 1
    vec mapMethy; // J
    mat C ; //NobsGene x L
    vec missGeneIdx; // NmissGene
    vec obsGeneIdx; // NobsGene
    vec missMethyIdx; // NmissMethy
    vec obsMethyIdx; // NobsMethy
    int nMissGene ;
    int nMissMethy;
    int nSample;
    int nGene;
    int nMethy;
    bool noCFlag;
    int nC;
    field<vec> JkAll;

    mat methy;
    mat gene;
    mat geneM;
    mat geneMbar;
    vec omega;
    double sigma ;
    vec sigmak ;

    double tauM ;
    double tauMbar;
    double tauOmega;

    double piM;
    double piMbar;
    double piOmega;

    vec IM;
    vec IMbar;
    vec IOmega;

    vec gammaM;
    vec gammaMbar;
    vec gammaC;

    mat omegaAll; // cube is more effecient than field when usable
    mat gammaMAll;
    mat gammaMbarAll;
    mat sigmakAll;
    mat IMAll;
    mat IMbarAll;
    mat IOmegaAll;
    mat gammaCAll;
    vec tauMAll;
    vec tauMbarAll;
    vec tauOmegaAll;
    vec sigmaAll;
    vec piMAll;
    vec piMbarAll;
    vec piOmegaAll;

    void updateOmega(int k, int itr){

        int JkSize = JkAll(k).size();

        if(JkSize==1){
            int j0 = JkAll(k)(0);
            double tmpMethy2 = sum(methy.col(j0) % methy.col(j0));
            double tmpTau = tau0Omega;
            if(IOmega(j0)==1)
            tmpTau = tauOmega;

            vec gammaM_k(nGene-1);
            mat geneM_k(nSample, nGene-1);
            int tmpIdx =0;
            for (int l=0; l<nGene; l++){
                if(l!=k){
                    gammaM_k(tmpIdx) = gammaM(l);
                    geneM_k.col(tmpIdx) = geneM.col(l);
                    tmpIdx++;
                }
            }

            vec A = Y - geneM_k*gammaM_k - geneMbar*gammaMbar - C*gammaC;
            vec B = gene.col(k);
            double factor1 = 1.0/(pow(tmpTau,-2) + tmpMethy2*( pow(sigmak(k),-2) + pow(sigma,-2)*pow(as_scalar(gammaM(k)),2)));
            double factor2 = sum( (A*as_scalar(gammaM(k))*pow(sigma,-2) + pow(sigmak(k),-2)*B) % methy.col(j0) ) ;

            omega(j0) = factor1*factor2 + as_scalar(randn(1))*std::sqrt(factor1);
            if(abs(omega(j0)) > sigma0 || !isfinite(omega(j0)) ){
                omega(j0)=0;
                IOmega(j0)=0;
            }
            omegaAll(itr,j0)=omega(j0);
        }
        else{
            for (int j0=JkAll(k)(0); j0<=JkAll(k)(JkSize-1); j0++){

                vec idxJk_j0(JkSize-1);
                int tmpIdx =0;
                for (int j=JkAll(k)(0); j<=JkAll(k)(JkSize-1); j++){
                    if(j!=j0){
                        idxJk_j0(tmpIdx) = j;
                        tmpIdx++;
                    }
                }
                double tmpMethy2 = as_scalar(methy.col(j0).t() * methy.col(j0));
                vec methyOmegaJk_j0 = methy.cols(idxJk_j0(0), idxJk_j0(idxJk_j0.size()-1)) * omega.subvec(idxJk_j0(0), idxJk_j0(idxJk_j0.size()-1));
                double tmpTau = tau0Omega;
                if(IOmega(j0)==1)
                tmpTau = tauOmega;

                vec gammaM_k(nGene-1);
                mat geneM_k(nSample, nGene-1);
                tmpIdx =0;
                for (int l=0; l<nGene; l++){
                    if(l!=k){
                        gammaM_k(tmpIdx) = gammaM(l);
                        geneM_k.col(tmpIdx) = geneM.col(l);
                        tmpIdx++;
                    }
                }
                vec A = Y - geneM_k*gammaM_k - geneMbar*gammaMbar - C*gammaC - methyOmegaJk_j0 * as_scalar(gammaM(k));
                vec B = gene.col(k) - methyOmegaJk_j0;
                double factor1 = 1.0/(pow(tmpTau,-2) + tmpMethy2*( pow(sigmak(k),-2) + pow(sigma,-2)*pow(as_scalar(gammaM(k)),2)));
                double factor2 = as_scalar( (A*as_scalar(gammaM(k))*pow(sigma,-2) + pow(sigmak(k),-2)*B).t()*methy.col(j0) ) ;

                omega(j0) = factor1*factor2 + as_scalar(randn(1))*std::sqrt(factor1);
                if(abs(omega(j0)) > ae || !isfinite(omega(j0)) ){
                    omega(j0)=0;
                }
                omegaAll(itr,j0)=omega(j0);

            }
        }
    }

    void updateGeneM(int k){
        if(JkAll(k).size()==1)
        geneM.col(k) = methy.col(JkAll(k)(0)) * omega(JkAll(k)(0));
        else
        geneM.col(k) = methy.cols(JkAll(k)(0), JkAll(k)(JkAll(k).size()-1)) * omega.subvec(JkAll(k)(0), JkAll(k)(JkAll(k).size()-1));
    }
    void updateGeneMbar(int k){
        geneMbar.col(k) = gene.col(k) - geneM.col(k);
        if(nMissGene>0){
            vec Ymissed(nMissGene);
            mat Cmissed(nMissGene,nC);
            mat geneMmissed(nMissGene, nGene);
            mat geneMbar_kmiss(nMissGene, nGene-1);
            int tmpIdx =0;
            vec gammaMbar_k(nGene-1);
            mat geneMbar_k(nSample, nGene-1);
            for (int l=0; l<nGene; l++){
                if(l!=k){
                    gammaMbar_k(tmpIdx) = gammaMbar(l);
                    geneMbar_k.col(tmpIdx) = geneMbar.col(l);
                    tmpIdx++;
                }
            }
            for(int l=0; l<nMissGene; l++){
                Ymissed(l)=Y(missGeneIdx(l));
                Cmissed.row(l) = C.row(missGeneIdx(l));
                geneMmissed.row(l) = geneM.row(missGeneIdx(l));
                geneMbar_kmiss.row(l) = geneMbar_k.row(missGeneIdx(l));
            }
            vec A = Ymissed - Cmissed*gammaC - geneMmissed * gammaM - geneMbar_kmiss*gammaMbar_k;
            double factor1 = 1.0/(pow(gammaMbar(k),2) * pow(sigma,-2) + pow(sigmak(k),-2));
            vec factor1_2 = A * gammaMbar(k) * pow(sigma,-2) * factor1;
            for(int l=0; l<nMissGene; l++)
            geneMbar(missGeneIdx(l),k) = factor1_2(l) + as_scalar(randn(1)) * std::sqrt(factor1);
        }
    }

    void updateGeneMiss(int k){
        for(int l=0; l<nMissGene; l++)
        gene(missGeneIdx(l),k) = geneM(missGeneIdx(l),k) + geneMbar(missGeneIdx(l),k);
    }
    void updateMethyMiss(int k){
        int JkSize = JkAll(k).size();

        vec Ymissed(nMissMethy);
        mat Cmissed(nMissMethy,nC);
        mat genemissed(nMissMethy, nGene);
        mat geneM_k(nSample, nGene-1);
        mat geneM_kmissed(nMissMethy, nGene-1);
        mat methymissed(nMissMethy, nMethy);

        vec gammaMbar_k(nGene-1);
        vec gammaM_k(nGene-1);

        int tmpIdx =0;
        for (int l=0; l<nGene; l++){
            if(l!=k){
                gammaM_k(tmpIdx) = gammaM(l);
                gammaMbar_k(tmpIdx) = gammaMbar(l);
                geneM_k.col(tmpIdx) = geneM.col(l);
                tmpIdx++;
            }
        }
        for(int l=0; l<nMissMethy; l++){
            Ymissed(l)=Y(missMethyIdx(l));
            Cmissed.row(l) = C.row(missMethyIdx(l));
            genemissed.row(l) = gene.row(missMethyIdx(l));
            geneM_kmissed.row(l) = geneM_k.row(missMethyIdx(l));
            methymissed.row(l) = methy.row(missMethyIdx(l));
        }

        if(JkSize==1){
            int j0=JkAll(k)(0);
            double tmpGammaM_Mbar = gammaM(k) - gammaMbar(k);
            double tmpOmega2 = pow(omega(j0),2);

            vec A = Ymissed - genemissed*gammaMbar - geneM_kmissed*(gammaM_k - gammaMbar_k) - Cmissed*gammaC;
            vec B = genemissed.col(k);

            double factor1 = 1.0/(pow(sigmam0,-2) + tmpOmega2* (pow(sigmak(k), -2) + pow(sigma, -2)*pow(tmpGammaM_Mbar,2)));
            vec factor2 =  (A*tmpGammaM_Mbar*pow(sigma,-2) + pow(sigmak(k),-2)*B)*omega(j0);

            for(int l=0; l<nMissMethy; l++)
            methy(missMethyIdx(l),j0) = as_scalar( factor1*factor2(l) + randn(1)*std::sqrt(factor1));
        }
        else{
            for (int j0=JkAll(k)(0); j0<=JkAll(k)(JkSize-1); j0++){

                vec idxJk_j0(JkSize-1);
                int tmpIdx =0;
                for (int j=JkAll(k)(0); j<=JkAll(k)(JkSize-1); j++){
                    if(j!=j0){
                        idxJk_j0(tmpIdx) = j;
                        tmpIdx++;
                    }
                }

                double tmpGammaM_Mbar = gammaM(k) - gammaMbar(k);
                double tmpOmega2 = pow(omega(j0),2);

                vec methyOmegaJk_j0 = methymissed.cols(idxJk_j0(0), idxJk_j0(idxJk_j0.size()-1)) * omega.subvec(idxJk_j0(0), idxJk_j0(idxJk_j0.size()-1));
                vec A = Ymissed - genemissed*gammaMbar - geneM_kmissed*(gammaM_k - gammaMbar_k) - methyOmegaJk_j0*(gammaM(k)-gammaMbar(k)) - Cmissed*gammaC;
                vec B = genemissed.col(k) - methyOmegaJk_j0;

                double factor1 = 1.0/(pow(sigmam0,-2) + tmpOmega2* (pow(sigmak(k), -2) + pow(sigma, -2)*pow(tmpGammaM_Mbar,2)));
                vec factor2 =  (A*tmpGammaM_Mbar*pow(sigma,-2) + pow(sigmak(k),-2)*B)*omega(j0) ;

                for(int l=0; l<nMissMethy; l++)
                methy(missMethyIdx(l),j0) = as_scalar( factor1*factor2(l) + randn(1)*std::sqrt(factor1));
            }
        }
    }

    void updateGammaC(int l, int itr){
        if(noCFlag){
            gammaC(l)=0;
        }
        else{
            double tmpC2 = sum(C.col(l) % C.col(l));
            double factor1;
            if(nC>1){
                mat C_l(nSample, nC-1);
                vec gammaC_l(nC-1);
                int tmpIdx =0;
                for (int m=0; m<nC; m++){
                    if(m!=l){
                        gammaC_l(tmpIdx) = gammaC(m);
                        C_l.col(tmpIdx) = C.col(m);
                        tmpIdx++;
                    }
                }
                factor1 = sum(C.col(l) % (Y-geneM*gammaM - geneMbar*gammaMbar - C_l*gammaC_l ));
            }
            else{
                factor1 = sum(C.col(l) % (Y-geneM*gammaM - geneMbar*gammaMbar) );
            }
            double factor2 = 1.0/(pow(sigma0,-2) + pow(sigma,-2)*tmpC2);
            gammaC(l) =  factor1*factor2*pow(sigma,-2) + as_scalar(randn(1)) * std::sqrt(factor2);
            gammaCAll(itr,l) = gammaC(l);
        }
    }


    void updateGammaM(int k, int itr){
        double tauTmp = tau0;
        if(IM(k)==1){
            tauTmp = tauM;
        }
        double tmpGene2 = sum(geneM.col(k) % geneM.col(k));

        int tmpIdx =0;
        vec gammaM_k(nGene-1);
        mat geneM_k(nSample, nGene-1);
        for (int l=0; l<nGene; l++){
            if(l!=k){
                gammaM_k(tmpIdx) = gammaM(l);
                geneM_k.col(tmpIdx) = geneM.col(l);
                tmpIdx++;
            }
        }
        double factor1 = sum(geneM.col(k) % (Y-geneMbar*gammaMbar - geneM_k*gammaM_k - C*gammaC ));
        double factor2 = 1.0/(pow(sigma,-2)*tmpGene2 + pow(tauTmp,-2));
        gammaM(k) =  factor1*factor2*pow(sigma,-2) + as_scalar(randn(1)) * std::sqrt(factor2);
        if(abs(gammaM(k)) > ae || !isfinite(gammaM(k)) ){
            gammaM(k)=0;
            IM(k)=0;
        }
        gammaMAll(itr, k) = gammaM(k);
    }

    void updateGammaMbar(int k, int itr){
        double tauTmp = tau0;
        if(IMbar(k)==1){
            tauTmp = tauMbar;
        }
        double tmpGene2 = sum(geneMbar.col(k) % geneMbar.col(k));

        int tmpIdx =0;
        vec gammaMbar_k(nGene-1);
        mat geneMbar_k(nSample, nGene-1);
        for (int l=0; l<nGene; l++){
            if(l!=k){
                gammaMbar_k(tmpIdx) = gammaMbar(l);
                geneMbar_k.col(tmpIdx) = geneMbar.col(l);
                tmpIdx++;
            }
        }

        double factor1 = sum(geneMbar.col(k) % (Y-geneM*gammaM - geneMbar_k*gammaMbar_k - C*gammaC ));
        double factor2 = 1.0/(pow(sigma,-2)*tmpGene2 + pow(tauTmp,-2));
        gammaMbar(k) =  factor1*factor2*pow(sigma,-2) + as_scalar(randn(1)) * std::sqrt(factor2);
        if(abs(gammaMbar(k)) > ae || !isfinite(gammaMbar(k)) ){
            gammaMbar(k)=0;
            IMbar(k)=0;
        }
        gammaMbarAll(itr, k) = gammaMbar(k);
    }

    void updateIM (int k, int itr){
        double f1Tmp = dNorm(gammaM(k), 0, tauM);
        double f0Tmp = dNorm(gammaM(k), 0, tau0);
        double pTmp = piM*f1Tmp / ((1-piM)*f0Tmp + piM*f1Tmp);
        if(f1Tmp == 0 && f0Tmp ==0){
            pTmp = 1;
        }
        double tmpU = as_scalar(randu());
        IM(k)=0;
        IMAll(itr,k)=0;
        if(tmpU < pTmp){
            IM(k) = 1;
            IMAll(itr,k)=1;

            double tmpSum=0;
            for(int j=JkAll(k)(0); j<=JkAll(k)(JkAll(k).size()-1); j++){
                tmpSum += IOmega(j);
            }
            if(tmpSum==0){
                IM(k)=0;
                IMAll(itr,k)=0;
            }
        }
    }

    void updateIMbar (int k, int itr){
        double f1Tmp = dNorm(gammaMbar(k), 0, tauMbar);
        double f0Tmp = dNorm(gammaMbar(k), 0, tau0);
        double pTmp = piMbar*f1Tmp / ((1-piMbar)*f0Tmp + piMbar*f1Tmp);
        if(f1Tmp == 0 && f0Tmp ==0){
            pTmp = 1;
        }
        double tmpU = as_scalar(randu());
        IMbar(k)=0;
        IMbarAll(itr,k)=0;
        if(tmpU < pTmp){
            IMbar(k) = 1;
            IMbarAll(itr,k)=1;
        }
    }

    void updateIOmega(int j, int itr){
        int k = mapMethy(j);
        if(!isfinite(omega(j)))
        cout << "omega " << k <<endl;
        double f1Tmp = dNorm(omega(j), 0, tauOmega);
        double f0Tmp = dNorm(omega(j), 0, tau0Omega);
        double pTmp = piOmega*f1Tmp / ((1-piOmega)*f0Tmp + piOmega*f1Tmp);
        if(f1Tmp == 0 && f0Tmp ==0){
            pTmp = 1;
        }
        double tmpU = as_scalar(randu());
        IOmega(j)=0;
        IOmegaAll(itr,j)=0;
        if(tmpU < pTmp){
            IOmega(j) = 1;
            IOmegaAll(itr,j)=1;
        }
    }

    void updatePiM (int itr){
        double piTmp = rBeta(1+sum(IM), 1+nGene-sum(IM));
	piM = piTmp;
        piMAll(itr) = piM;
    }
    void updatePiMbar (int itr){
        double piTmp = rBeta(1+sum(IMbar), 1+nGene-sum(IMbar));
	        piMbar = piTmp;
	  piMbarAll(itr) = piMbar;
    }
    void updatePiOmega (int itr){
        double piTmp = rBeta(1+sum(IOmega), 1+nMethy-sum(IOmega));
	        piOmega = piTmp;
        piOmegaAll(itr) = piOmega;
    }

    void updateSigma(int itr){
        vec factor1 = Y - geneMbar*gammaMbar - geneM*gammaM - C*gammaC;
        double factor1sq = sum(factor1%factor1);
        sigma = std::sqrt(rInvGamma(nSample/2.0+smallNumber, factor1sq/2.0+smallNumber));
        sigmaAll(itr) = sigma;
    }

    void updateSigmak(int k, int itr){
        double factor1 = sum(geneMbar.col(k) % geneMbar.col(k));
        sigmak(k)=std::sqrt(rInvGamma(nSample/2.0 + smallNumber, factor1/2.0 + smallNumber));
        sigmakAll(itr,k) = sigmak(k);
    }


    void updateTauM(int itr){
      double a=0;
        double b = 0;
        for (int k=0; k<nGene; k++){
            b += pow(gammaM(k), 2)*IM(k);
	    a += IM(k);
        }
	tauM = std::sqrt(rInvGamma(a/2.0+tauAlpha, b/2.0+tauBeta));
	if(abs(tauM) > std::sqrt(ae) || !isfinite(tauM) ){
	  tauM =std::sqrt(ae);
	}
        tauMAll(itr) = tauM;
    }

    void updateTauMbar(int itr){
        double b = 0;
	double a=0;
        for (int k=0; k<nGene; k++){
            b += pow(gammaMbar(k), 2)*IMbar(k);
	    a += IMbar(k);
        }
	tauMbar = std::sqrt(rInvGamma(a/2.0+tauAlpha, b/2.0+tauBeta));
	if(abs(tauMbar) > std::sqrt(ae) || !isfinite(tauMbar) ){
	  tauMbar =std::sqrt(ae);
	}
        tauMbarAll(itr) = tauMbar;
    }

    void updateTauOmega(int itr){
        double b = 0;
	double a =0;
        for (int j=0; j<nMethy; j++){
            int k = mapMethy(j);
            b += pow(omega(j), 2)*IOmega(k);
	    a+= IOmega(k);
        }
	tauOmega = std::sqrt(rInvGamma(a/2.0+tauAlpha, b/2.0+tauBeta));
	if(abs(tauOmega) > std::sqrt(ae) || !isfinite(tauOmega) ){
	  tauOmega =std::sqrt(ae);
	}
        tauOmegaAll(itr) = tauOmega;
    }

    void Gibbs(){

        if(nItr > itrCheckConvergence){ // check convergence after itrStart=800
            for (int itr=0; itr<itrStart; itr++){
                for(int k=0; k<nGene; k++){
                    updateOmega(k,itr);
                }

                for(int j=0; j<nMethy; j++){
                    updateIOmega(j, itr);
                }
                for(int k=0; k<nGene; k++){
                    updateIM(k, itr);
                }

                for(int k=0; k<nGene; k++){
                    updateIMbar(k, itr);
                }

                updatePiM(itr);
                updatePiMbar(itr);
                updatePiOmega(itr);

                for(int k=0; k<nGene; k++){
                    updateGeneM(k);
                }
                for(int k=0; k<nGene; k++){
                    updateGeneMbar(k);
                }

                if(nMissGene>0){
                    for(int k=0; k<nGene; k++){
                        updateGeneMiss(k);
                    }
                }
                if(nMissMethy>0){
                    for(int k=0; k<nGene; k++){
                        updateMethyMiss(k);
                    }
                }
                for(int k=0; k<nGene; k++){
                    updateSigmak(k,itr);
                }
                for(int k=0; k<nGene; k++){
                    updateGammaM(k,itr);
                }
                for(int k=0; k<nGene; k++){
                    updateGammaMbar(k,itr);
                }
                for(int l=0; l<nC; l++){
                    updateGammaC(l,itr);
                }

                updateSigma(itr);
                updateTauOmega(itr);
                updateTauM(itr);
                updateTauMbar(itr);
            }
            int itr=itrStart;
            vec windowHead = sigmaAll.subvec(itrStart-201, itrStart-181);
            vec windowTail = sigmaAll.subvec(itrStart-101, itrStart-1);
            while(itr < nItr && !naiveGeweke(windowHead,windowTail)){
                for(int k=0; k<nGene; k++){
                    updateOmega(k,itr);
                }

                for(int j=0; j<nMethy; j++){
                    updateIOmega(j, itr);
                }
                for(int k=0; k<nGene; k++){
                    updateIM(k, itr);
                }

                for(int k=0; k<nGene; k++){
                    updateIMbar(k, itr);
                }

                updatePiM(itr);
                updatePiMbar(itr);
                updatePiOmega(itr);

                for(int k=0; k<nGene; k++){
                    updateGeneM(k);
                }
                for(int k=0; k<nGene; k++){
                    updateGeneMbar(k);
                }

                if(nMissGene>0){
                    for(int k=0; k<nGene; k++){
                        updateGeneMiss(k);
                    }
                }
                if(nMissMethy>0){
                    for(int k=0; k<nGene; k++){
                        updateMethyMiss(k);
                    }
                }
                for(int k=0; k<nGene; k++){
                    updateSigmak(k,itr);
                }
                for(int k=0; k<nGene; k++){
                    updateGammaM(k,itr);
                }
                for(int k=0; k<nGene; k++){
                    updateGammaMbar(k,itr);
                }
                for(int l=0; l<nC; l++){
                    updateGammaC(l,itr);
                }

                updateSigma(itr);
                updateTauOmega(itr);
                updateTauM(itr);
                updateTauMbar(itr);
                if( itr%20==0 ){
                    windowHead = sigmaAll.subvec(itr-200, itr-180);
                    windowTail = sigmaAll.subvec(itr-100, itr);
                }
                itr++;
            }
            itr--;
            if(itr < nItr){
                nItrActual = itr;
                nBurninActual = itr-200;
            }
            else{
                nItrActual = nItr;
                nBurninActual = std::round(nItr/2.0);
            }
        }
        else{
            for(int itr=0; itr<nItr; itr++){
                for(int k=0; k<nGene; k++){
                    updateOmega(k,itr);
                }

                for(int j=0; j<nMethy; j++){
                    updateIOmega(j, itr);
                }
                for(int k=0; k<nGene; k++){
                    updateIM(k, itr);
                }

                for(int k=0; k<nGene; k++){
                    updateIMbar(k, itr);
                }

                updatePiM(itr);
                updatePiMbar(itr);
                updatePiOmega(itr);

                for(int k=0; k<nGene; k++){
                    updateGeneM(k);
                }
                for(int k=0; k<nGene; k++){
                    updateGeneMbar(k);
                }

                if(nMissGene>0){
                    for(int k=0; k<nGene; k++){
                        updateGeneMiss(k);
                    }
                }
                if(nMissMethy>0){
                    for(int k=0; k<nGene; k++){
                        updateMethyMiss(k);
                    }
                }
                for(int k=0; k<nGene; k++){
                    updateSigmak(k,itr);
                }
                for(int k=0; k<nGene; k++){
                    updateGammaM(k,itr);
                }
                for(int k=0; k<nGene; k++){
                    updateGammaMbar(k,itr);
                }
                for(int l=0; l<nC; l++){
                    updateGammaC(l,itr);
                }

                updateSigma(itr);
                updateTauOmega(itr);
                updateTauM(itr);
                updateTauMbar(itr);
            }
            nItrActual = nItr;
            nBurninActual = std::round(nItr/2.0);
        }
    }
};


// [[Rcpp::export]]
Rcpp::List FBMcpp(const Rcpp::List& data, int nItr, const int seed, const string method){
    //extract matrices from list
    parmUpdate parm;

    parm.nItr = nItr;

    vec Y = Rcpp::as<vec>(data(0)); // N x 1
    parm.Y = Y;
    mat geneObs= Rcpp::as<mat>(data(1)); // NobsGene x K
    mat methyObs = Rcpp::as<mat>(data(2)); // NobsMethy x J
    parm.mapMethy = Rcpp::as<vec>(data(3))-1; // J
    mat C = Rcpp::as<mat>(data(4)); //NobsGene x L
    parm.C=C;
    parm.missGeneIdx = Rcpp::as<vec>(data(5))-1; // NmissGene
    parm.obsGeneIdx = Rcpp::as<vec>(data(6))-1; // NobsGene
    parm.missMethyIdx = Rcpp::as<vec>(data(7))-1; // NmissMethy
    parm.obsMethyIdx = Rcpp::as<vec>(data(8))-1; // NobsMethy

    srand(seed);
    arma_rng::set_seed(seed);

    int nMissGene = parm.missGeneIdx.size();
    parm.nMissGene = nMissGene;
    int nMissMethy = parm.missMethyIdx.size();
    parm.nMissMethy = nMissMethy;
    int nSample = Y.size();
    parm.nSample = nSample;
    int nGene= geneObs.n_cols;
    parm.nGene = nGene;
    int nMethy=methyObs.n_cols;
    parm.nMethy = nMethy;
    parm.noCFlag = false;
    int nC = C.n_cols;
    parm.nC = nC;
    if(parm.nC==1){
        vec tmpSd = arma::stddev(C);
        if(tmpSd(0)==0.0){
            parm.noCFlag=true;
        }
    }

    field<vec> JkAll(nGene);
    for (int k=0; k<nGene; k++){
        JkAll[k] = zeros<vec>(nMethy);
        int tmpCounter=0;
        for (int j=0; j<nMethy; j++){
            if(parm.mapMethy(j)==k){
                JkAll(k)(tmpCounter)=j;
                tmpCounter++;
            }
        }
        JkAll(k).resize(tmpCounter);
    }
    parm.JkAll = JkAll;


    // Initialize parameters.
    mat methyMiss = zeros<mat>(nMissMethy, nMethy);
    mat tmpMean = arma::mean(methyObs);
    for(int j=0; j<nMethy; j++){
        methyMiss.col(j).fill(tmpMean(j));
    }
    mat methy = zeros<mat>(nSample, nMethy);
    for(int i=0; i<(nSample-nMissMethy); i++){
        methy.row(parm.obsMethyIdx(i)) = methyObs.row(i); //assuming sample idx are sorted
    }
    for(int i=0; i<nMissMethy; i++){
        methy.row(parm.missMethyIdx(i)) = methyMiss.row(i); //assuming sample idx are sorted
    }
    parm.methy = methy;

    mat geneMiss = zeros<mat>(nMissGene, nGene);
    mat gene = zeros<mat>(nSample, nGene);
    for(int i=0; i<(nSample-nMissGene); i++){
        gene.row(parm.obsGeneIdx(i)) = geneObs.row(i); //assuming sample idx are sorted
    }
    for(int i=0; i<nMissGene; i++){
        gene.row(parm.missGeneIdx(i)) = geneMiss.row(i); //assuming sample idx are sorted
    }
    parm.gene=gene;

    mat geneM = zeros<mat>(nSample, nGene);
    mat geneMbar = zeros<mat>(nSample, nGene);
    vec omega = randn<vec>(nMethy);
    omega.fill(5);

    for(int k=0; k<nGene; k++){
        int tmpIdx1 = parm.JkAll(k)(0);
        int tmpIdx2 = parm.JkAll(k)(parm.JkAll(k).size()-1);
        geneM.col(k) = methy.cols(tmpIdx1, tmpIdx2) * omega.subvec(tmpIdx1,tmpIdx2); // (N,Jk) x (Jk,1)
        geneMbar.col(k) = gene.col(k) - geneM.col(k);
    }

    parm.geneM = geneM;
    parm.geneMbar = geneMbar;
    parm.omega=omega;

    parm.sigma = as_scalar(randu()*4)+1;
    parm.sigmak = randu<vec>(nGene)*4+1;

    parm.tauM = as_scalar(randu()*4)+1;
    parm.tauMbar = as_scalar(randu()*4)+1;
    parm.tauOmega = as_scalar(randu()*4)+1;

    parm.piM = as_scalar(randu());
    parm.piMbar = as_scalar(randu());
    parm.piOmega = as_scalar(randu());

    double tmpP = as_scalar(randu());
    parm.IM = zeros<vec>(nGene);
    parm.IMbar = zeros<vec>(nGene);
    parm.IOmega = zeros<vec>(nMethy);
    if(tmpP < parm.piM){
        parm.IM.fill(1);
    }
    if(tmpP < parm.piMbar){
        parm.IMbar.fill(1);
    }
    if(tmpP < parm.piOmega){
        parm.IOmega.fill(1);
    }

    parm.gammaM = randn<vec>(nGene)*10;
    parm.gammaMbar = randn<vec>(nGene)*parm.tauMbar;
    parm.gammaC = randn<vec>(nC)*5;
    parm.gammaM.fill(0);
    parm.gammaMbar.fill(0);
    parm.gammaC.fill(10);

    if(nSample==1){
        printf("ERROR: More than one sample is needed.\n");
        exit(1);
    }
    if( nGene<1 || nMethy<1 || nC <1){
        printf("ERROR: No gene/methylation/clinical variable.\n");
        exit(1);
    }

    /// storing iterations ///
    mat omegaAll = zeros<mat>(nItr, nMethy); // cube is more effecient than field when usable
    mat gammaMAll = zeros<mat>(nItr, nGene);
    mat gammaMbarAll = zeros<mat>(nItr, nGene);
    mat sigmakAll = zeros<mat>(nItr, nGene);
    mat IMAll = zeros<mat>(nItr, nGene);
    mat IMbarAll = zeros<mat>(nItr, nGene);
    mat IOmegaAll = zeros<mat>(nItr, nMethy);
    mat gammaCAll = zeros<mat>(nItr, nC);
    vec tauMAll = zeros<vec>(nItr);
    vec tauMbarAll = zeros<vec>(nItr);
    vec tauOmegaAll = zeros<vec>(nItr);
    vec sigmaAll = zeros<vec>(nItr);
    vec piMAll = zeros<vec>(nItr);
    vec piMbarAll = zeros<vec>(nItr);
    vec piOmegaAll = zeros<vec>(nItr);

    parm.omegaAll = omegaAll; // cube is more effecient than field when usable
    parm.gammaMAll = gammaMAll;
    parm.gammaMbarAll = gammaMbarAll;
    parm.sigmakAll = sigmakAll;
    parm.IMAll = IMAll;
    parm.IMbarAll = IMbarAll;
    parm.IOmegaAll = IOmegaAll;
    parm.gammaCAll = gammaCAll;
    parm.tauMAll = tauMAll;
    parm.tauMbarAll = tauMbarAll;
    parm.tauOmegaAll = tauOmegaAll;
    parm.sigmaAll = sigmaAll;
    parm.piMAll = piMAll;
    parm.piMbarAll = piMbarAll;
    parm.piOmegaAll = piOmegaAll;

    ////////// update ///////////
    printf("Running Gibbs sampling...\n");
    parm.Gibbs();
    printf("Gibbs sampling completed. \n");

    Rcpp::List result;

    result["omega"] = parm.omegaAll;
    result["sigma"] = parm.sigmaAll;
    result["sigmak"] = parm.sigmakAll;
    result["gammaM"] = parm.gammaMAll;
    result["gammaMbar"] = parm.gammaMbarAll;
    result["tauM"] = parm.tauMAll;
    result["tauMbar"] = parm.tauMbarAll;
    result["tauOmega"] = parm.tauOmegaAll;
    result["IOmega"] = parm.IOmegaAll;
    result["IM"] = parm.IMAll;
    result["IMbar"] = parm.IMbarAll;
    result["gammaC"] = parm.gammaCAll;
    result["nItrActual"] = parm.nItrActual;
    result["nBurninActual"] = parm.nBurninActual;
result["piM"] =parm.piMAll;
result["piMbar"] =parm.piMbarAll;

    return result;
}

