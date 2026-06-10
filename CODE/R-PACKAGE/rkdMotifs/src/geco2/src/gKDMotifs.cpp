#include "gKDMotifs.h"
#include "gNMF.h"
#include "gPWM.h"

#ifdef GECO_HAS_CUDA
#include "gMethods_cuda.h"
#endif
#include "gMethods_armadillo.h"

#include <math.h>
#include <chrono>
#include <map>
#include <unordered_map>
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>

// nlopt include
#include <nlopt.hpp>
// autodiff include
#include "autodiff/forward/real.hpp"
#include "autodiff/forward/real/eigen.hpp"

using namespace std;
using namespace geco;
using namespace geco::bed;
using namespace geco::methods;
using namespace geco::methods::kmers;
using namespace geco::methods::kdmotifs;
using namespace geco::methods::cluster;
using namespace geco::methods::nmf;

using namespace autodiff;
using namespace nlopt;
using namespace Eigen;


//autodiff & eigen & nlopt stuff
typedef struct {
    ArrayXreal x;
    ArrayXreal y;
} func_data;


autodiff::real j4(const ArrayXreal& x,const ArrayXreal& y, const ArrayXreal& p){
    return (y-p[0]-p[1]*(exp( -((x-p[2])/p[3]).square()) + exp( -((x+p[2])/p[3]).square()))).square().sum();
}

autodiff::real j3(const ArrayXreal& x,const ArrayXreal& y, const ArrayXreal& p){
    return ( y-p[0]-p[1]*(2*exp( -(x/p[2]).square())) ).square().sum();
}

double myvfunc(const std::vector<double> &p, std::vector<double> & grad, void * mfd){
    autodiff::real u;
    func_data & data = *(func_data * ) mfd;
    ArrayXreal ep(p.size());
    for(size_t i=0;i<p.size();i++){
        ep[i]=p[i];
    }
    
    if(grad.size()>0){
        Eigen::VectorXd vec;
        if(p.size()==4){
            vec   = gradient(j4, wrt(ep), at(data.x, data.y, ep), u);
        }else if(p.size()==3){
            vec   = gradient(j3, wrt(ep), at(data.x, data.y, ep), u);
        }else throw gMethodException("Invalid parameters number");
        for(size_t i=0;i<p.size();i++){
            grad[i]=vec[i];
        }
    }else{
        if(p.size()==4){
            u=j4(data.x,data.y,ep);
        }else if(p.size()==3){
            u=j3(data.x,data.y,ep);            
        }else throw gMethodException("Invalid parameters number");
    }
    return (double) u;
}

double AIC(double fmin,double N,double k){
    //see in biblio for reference
    return N*log(fmin/N) +2*(k+1);
}

double myIntegral(double x,double peak, double width){
    return (erf((peak+x)/width)-erf((peak-x)/width))/2;
}


arma::vec getA_internal(const vector<double> & x, const vector<double> & y,double intv,double win,size_t nrep,unsigned int maxSelectedPeak,unsigned int half90interval){
    
    func_data data;
    data.x=ArrayXreal(x.size());
    data.y=ArrayXreal(y.size());
    
    for(size_t i=0;i<x.size();i++){
        data.x[i]=x[i];
        data.y[i]=y[i];
    }
   
    vector<double> lb(4),ub(4);
    
    // baseline delta peak width
    lb[0]=(*min_element(y.begin(),y.end()))*1.01;
    ub[0]=(*max_element(y.begin(),y.end()))*1.01;
    
    lb[1]=0;
    ub[1]=ub[0]-lb[0];
    
    lb[2]=0;
    ub[2]=intv;
    
    lb[3]=win/2;
    ub[3]=intv/2;
    
    nlopt::opt opt(nlopt::LD_MMA, 4);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_min_objective(myvfunc, &data);
    opt.set_ftol_rel(1e-6);
    opt.set_xtol_rel(1e-4);
    
    std::mt19937 rng(0);
    std::uniform_real_distribution<double> gen0(lb[0],ub[0]);
    std::uniform_real_distribution<double> gen1(lb[1],ub[1]);
    std::uniform_real_distribution<double> gen2(lb[2],ub[2]);
    std::uniform_real_distribution<double> gen3(lb[3],ub[3]);
    
    vector<double> par4(4);
    double minf4=HUGE_VAL;
    for(size_t i=0;i<nrep;i++){
        try{
            double pminf;
            vector<double> par0(4,0);
            par0[0] = gen0(rng);
            par0[1] = gen1(rng);
            par0[2] = gen2(rng);
            par0[3] = gen3(rng);
            
            vector<double> tpar=par0;
            auto result = opt.optimize(tpar, pminf);
            if(pminf<minf4){
                minf4=pminf;
                par4=tpar;
            }
        }catch(std::exception &e) {
            
        }
    }

    lb=vector<double>(3);
    ub=vector<double>(3);
    
    // baseline delta peak width
    lb[0]=(*min_element(y.begin(),y.end()))*1.01;
    ub[0]=(*max_element(y.begin(),y.end()))*1.01;
    
    lb[1]=lb[0]-ub[0];
    ub[1]=ub[0]-lb[0];
    
    lb[2]=win/2;
    ub[2]=intv;
    
    nlopt::opt opt3(nlopt::LD_MMA, 3);
    opt3.set_lower_bounds(lb);
    opt3.set_upper_bounds(ub);
    opt3.set_min_objective(myvfunc, &data);
    opt3.set_ftol_rel(1e-6);
    opt3.set_xtol_rel(1e-4);
    
    rng=std::mt19937(0);
    gen0=std::uniform_real_distribution<double>(lb[0],ub[0]);
    gen1=std::uniform_real_distribution<double>(lb[1],ub[1]);
    gen2=std::uniform_real_distribution<double>(lb[2],ub[2]);

    
    vector<double> par3(3);
    double minf3=HUGE_VAL;
    for(size_t i=0;i<nrep;i++){
        try{
            double pminf;
            vector<double> par0(3,0);
            par0[0] = gen0(rng);
            par0[1] = gen1(rng);
            par0[2] = gen2(rng);
            
            vector<double> tpar=par0;
            auto result = opt3.optimize(tpar, pminf);
            if(pminf<minf3){
                minf3=pminf;
                par3=tpar;
            }
        }catch(std::exception &e) {
            
        }
    }
    
    arma::Col<double> ax=arma::conv_to< arma::vec >::from(x);
    arma::Col<double> ay=arma::conv_to< arma::vec >::from(y);
    
    arma::Col<double> p4=(arma::exp( -arma::square( (ax-par4[2])/par4[3]))  +  arma::exp(-arma::square( (ax+par4[2])/par4[3])) ) * par4[1] + par4[0];
    arma::Col<double> p3=2*arma::exp( -arma::square(ax/par3[2]) ) * par3[1] + par3[0];
    
    double srd4 =arma::stddev(ay-p4);
    double srd3 =arma::stddev(ay-p3);

    arma::Col<double> c4=cor(ay,p4);
    arma::Col<double> c3=cor(ay,p3);
    
    arma::Col<double> pp4=abs(p4-par4[0]);
    arma::Col<double> pp3=abs(p3-par3[0]);
    double pm4=max(pp4)/srd4;
    double pm3=max(pp3)/srd3;
    
    double m4 = round(pm4);
    double m3 = round(pm3);

    double my=mean(ay);
    double ny=norm(ay-my);
    double sdr0 = arma::stddev(ay-my);
    
    double aic4=AIC(minf4,y.size(),4);
    double aic3=AIC(minf3,y.size(),3);
    double aic0=AIC(ny,y.size(),1);
    
    
    bool good5,good6;
    good5 = max(p4)>1+srd4;
    if(par3[1]>0){
        good6 = max(p3)>1+srd3;
    }else{
        good6 = min(p3)<1-srd3;
    }

//     cout << "----------------------------------------------------------" << endl;
//     cout << "aic: " << "\t" << aic4 << "\t" << aic3 << endl;        
//     cout << "par4: " << "\t" << par4[0] << "\t" << par4[1] << "\t" << par4[2] << "\t" << par4[3] << endl;
//     cout << "par3: " << "\t" << par3[0] << "\t" << par3[1] << "\t" << 0 << "\t" << par3[2] << endl;
//     cout << "srd"   << "\t" << srd4 << "\t" << srd3 << endl;
//     cout << "max p" << "\t" << max(p4) << "\t" << max(p3) << endl;
//     cout << "min p" << "\t" << min(p4) << "\t" << min(p3) << endl;
//     cout << "good"  << "\t" << good5 << "\t" << good6 << endl;
//     cout << "raw: "<< "\t" << pm4 << "\t" << pm3 << endl;
//     cout << "ratio" << "\t" << m4 << "\t" << m3 << endl;
//     cout << "----------------------------------------------------------" << endl;
    
    arma::vec res(12);
    //initially we set the parameters as if the best model is the mean of the profile
    res[0]=my;      //baseline
    res[1]=0;       //delta
    res[2]=0;       //peak
    res[3]=1e-6;    //width is not set at 0 to allow for model profile calculation
    res[4]=0;       //corr
    res[5]=aic0;    //aic
    res[6]=ny;      //fmin
    res[7]=sdr0;    //residue sd
    res[8]=0;       //where peak reach 95%
    res[9]=0;       //classification
    res[10]=0;       //selected
    res[11]=0;

    
    if( (round(aic4)<round(aic3)) && good5){
        if(m4 >= 3 ){
            res[0]=par4[0];
            res[1]=par4[1];
            res[2]=par4[2];
            res[3]=par4[3];
            res[4]=c4[0];
            res[5]=aic4;
            res[6]=minf4;
            res[7]=srd4;
            res[11]=arma::max(p4)-par4[0];
        }
    }else if(good6){
        if(m3 >= 3){
            res[0]=par3[0];
            res[1]=par3[1];
            res[2]=0;
            res[3]=par3[2];
            res[4]=c3[0];
            res[5]=aic3;
            res[6]=minf3;
            res[7]=srd3;
            if(par3[1]>0){
                res[11]=arma::max(p3)-par3[0];
            }else{
                res[11]=arma::min(p3)-par3[0];
            }
        }
    }
    double pupr = 2* pow(res[2],2) / pow(res[3],2);
    double erre  = (res[1]<=0)?(0):((erf((res[2]+half90interval)/res[3])-erf((res[2]-half90interval)/res[3]))/2);

    if(res[1]!=0){
        double pos=0;
        while(myIntegral(pos,res[2],res[3])<0.95){
            pos++;
        }
        res[8]=pos;
    }else{
        res[8]=INFINITY;
    }
    

    res[10] = erre>0.9 && res[4]>0.9 && res[2]<=maxSelectedPeak;
    
    if(floor(res[2])==0){
        if(res[1]>0){
            res[9]=1;                               //CENTRALLY ENCRICHED
        }else if(res[1]==0){
            if(res[0]>1){
                res[9]=5;                           //CONSTANTLY ENRICHED                
            }else{
                res[9]=7;                           //CONSTANTLY DEPLETED
            }
        }else{
            res[9]=6;                               //CENTRALLY DEPLETED            
        }
    }else{
        if(pupr<1){
            res[9]=2;                               //TWO MERGED PEAK
        }else if(res[1]>0){
            if(2*res[1]*exp(-pow(res[2]/res[3],2))>=res[7]){
                res[9]=3;                           //TWO OVERLAPPING PEAKS
            }else{
                res[9]=4;                           //TWO DISTINCT PEAKS
            }
        }
    }
    return res;
}

#ifdef GECO_HAS_CUDA
template <typename T> 
arma::Mat<T> calculateProfile_cuda(const gKDMotifs<T> & motifs,unsigned int win,unsigned int nwin,const vector<string> & pseq,const vector<string> nseq,bool strict){
    gMThreadException te;
    cuda::autoCudaHandles cudaHandles;
    cuda::cudaDense<T> pmat(nwin,motifs.i_ncols,0);
    cuda::cudaDense<T> nmat(nwin,motifs.i_ncols,0);
    cuda::cudaDense<T> W = cuda::convert<gMDense,cuda::cudaDense,T>(motifs,false,cudaHandles);
    
    vector<gMSparse<T>> pmats(pseq.size());
    vector<gMSparse<T>> nmats(pseq.size());
    
#pragma omp parallel for
    for(size_t i=0;i<pseq.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            pmats[i]=kdmDistrProfile<T>(pseq[i],win,motifs.i_l,motifs.i_k,motifs.i_doubleStrand,motifs.i_estimated,strict);
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    
#pragma omp parallel for
    for(size_t i=0;i<nseq.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            nmats[i]=kdmDistrProfile<T>(nseq[i],win,motifs.i_l,motifs.i_k,motifs.i_doubleStrand,motifs.i_estimated,strict);
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    
    for(size_t i=0;i<pseq.size();i++){
        if(pmats[i].i_ncols>0){
            cuda::multiply_cuda<T>(cuda::convert<gMSparse,cuda::cudaSparse,T>(pmats[i],false,cudaHandles),W,pmat,true,false,1,1,cudaHandles);
        }
    }

    for(size_t i=0;i<nseq.size();i++){
        if(nmats[i].i_ncols>0){
            cuda::multiply_cuda<T>(cuda::convert<gMSparse,cuda::cudaSparse,T>(nmats[i],false,cudaHandles),W,nmat,true,false,1,1,cudaHandles);
        }
    }

    cuda::elementwise_divide(pmat,nmat,pmat);
    return arma::trans(armadillo::convert<gMDense,arma::Mat,T>(cuda::convert<cuda::cudaDense,gMDense,T>(pmat,false,cudaHandles))*((((T) 1)/(T) nseq.size())/(((T) 1)/(T) pseq.size())));
}

template <typename T> 
pair<vector<gMDense<T>>,vector<gMDense<T>>> calculateProfile_cuda_3(const gKDMotifs<T> & motifs,unsigned int win,unsigned int nwin,const vector<string> & pseq,const vector<string> nseq,bool strict){
    gMThreadException te;
    cuda::autoCudaHandles cudaHandles;
    cuda::cudaDense<T> pmat(nwin,motifs.i_ncols,0);
    cuda::cudaDense<T> nmat(nwin,motifs.i_ncols,0);
    cuda::cudaDense<T> W = cuda::convert<gMDense,cuda::cudaDense,T>(motifs,false,cudaHandles);
    
    vector<gMSparse<T>> pmats(pseq.size());
    
#pragma omp parallel for
    for(size_t i=0;i<pseq.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            pmats[i]=kdmDistrProfile<T>(pseq[i],win,motifs.i_l,motifs.i_k,motifs.i_doubleStrand,motifs.i_estimated,strict);
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();

    vector<gMSparse<T>> nmats(nseq.size());    
    
#pragma omp parallel for
    for(size_t i=0;i<nseq.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            nmats[i]=kdmDistrProfile<T>(nseq[i],win,motifs.i_l,motifs.i_k,motifs.i_doubleStrand,motifs.i_estimated,strict);
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    vector<gMDense<T>> pfeat(pseq.size());
    vector<gMDense<T>> nfeat(nseq.size());
    
    for(size_t i=0;i<pseq.size();i++){
        if(pmats[i].i_ncols>0){
            cuda::multiply_cuda<T>(cuda::convert<gMSparse,cuda::cudaSparse,T>(pmats[i],false,cudaHandles),W,pmat,true,false,1,0,cudaHandles);
            pfeat[i]=cuda::convert<cuda::cudaDense,gMDense,T>(pmat,false,cudaHandles);
        }
    }

    for(size_t i=0;i<nseq.size();i++){
        if(nmats[i].i_ncols>0){
            cuda::multiply_cuda<T>(cuda::convert<gMSparse,cuda::cudaSparse,T>(nmats[i],false,cudaHandles),W,nmat,true,false,1,0,cudaHandles);
            nfeat[i]=cuda::convert<cuda::cudaDense,gMDense,T>(nmat,false,cudaHandles);
        }
    }
    auto res=make_pair(pfeat,nfeat);
    return move(res);
}

#endif

template <typename T> 
arma::Mat<T> calculateProfile_arma(const gKDMotifs<T> & motifs,unsigned int win,unsigned int nwin,const vector<string> & pseq,const vector<string> nseq,bool strict){
    gMThreadException te;
    arma::Mat<T> pmat(nwin,motifs.i_ncols,arma::fill::zeros);
    arma::Mat<T> nmat(nwin,motifs.i_ncols,arma::fill::zeros);
    arma::Mat<T> W = armadillo::convert<gMDense,arma::Mat,T>(motifs);
    
    vector<gMSparse<T>> pmats(pseq.size());
    vector<gMSparse<T>> nmats(pseq.size());
    
#pragma omp parallel for
    for(size_t i=0;i<pseq.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            pmats[i]=kdmDistrProfile<T>(pseq[i],win,motifs.i_l,motifs.i_k,motifs.i_doubleStrand,motifs.i_estimated,strict);
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();    
    
#pragma omp parallel for
    for(size_t i=0;i<nseq.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            nmats[i]=kdmDistrProfile<T>(nseq[i],win,motifs.i_l,motifs.i_k,motifs.i_doubleStrand,motifs.i_estimated,strict);
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();        
    
    for(size_t i=0;i<pseq.size();i++){
        pmat += trans(armadillo::convert<gMSparse,arma::SpMat,T>(toHellinger(pmats[i]))) * W;
        checkInterrupt();
    }

    for(size_t i=0;i<nseq.size();i++){
        nmat += trans(armadillo::convert<gMSparse,arma::SpMat,T>(toHellinger(nmats[i]))) * W ;
        checkInterrupt();
    }
    
    return arma::trans((pmat/nmat)*((((T) 1)/(T) nseq.size())/(((T) 1)/(T) pseq.size())));
}


template <typename T> 
pair<gMDense<T>,vector<T>> geco::methods::kdmotifs::getMotifsProfile2(const gKDMotifs<T> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,bool strict){
    geco::bed::gBedData<geco::bed::gBedBaseEntry> pos,neg;
    size_t base=0;
    for(auto r=regions.begin();r!=regions.end();r++){
        bool isF = r->i_isForward;        
        if(labels[base]>0){
//             cout << base << "\t" << ((isF)?("+"):("-")) << endl;
            pos.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,isF));
        }else{
            neg.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,isF));
        }
        base++;
    }
    vector<string> pseq=pos.getSequences(genomeFile);
    vector<string> nseq=neg.getSequences(genomeFile);
//     cout << pseq[0] << endl;
//     cout << nseq[0] << endl;
    
    unsigned int win = 2*halfWin +1;
    unsigned int nwin = 2*(halfInterval-halfWin)+1;

#ifdef GECO_HAS_CUDA
    arma::Mat<T> Y;    
    if(useCuda()){
        Y=calculateProfile_cuda(motifs,win,nwin,pseq,nseq,strict);
    }else{
        Y=calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);
    }
#else
    arma::Mat<T> Y = calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);    
#endif
        
    vector<T> x(Y.n_cols);
    x[0]=-(T)halfInterval+(T)halfWin;
    for(unsigned int i=1;i<x.size();i++){
        x[i]=x[i-1]+1;
    }
    
    return make_pair(armadillo::convert<arma::Mat,gMDense,T>(Y),x);    
}
template pair<gMDense<float>,vector<float>> geco::methods::kdmotifs::getMotifsProfile2(const gKDMotifs<float> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, bool strict);
template pair<gMDense<double>,vector<double>> geco::methods::kdmotifs::getMotifsProfile2(const gKDMotifs<double> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, bool strict);


template<typename T>
pair<pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T>> geco::methods::kdmotifs::getMotifsProfile3(const gKDMotifs<T> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict){
    unsigned int win = 2*halfWin +1;
    unsigned int nwin = 2*(halfInterval-halfWin)+1;
    vector<T> x(nwin);
    x[0]=-(T)halfInterval+(T)halfWin;
    for(unsigned int i=1;i<x.size();i++){
        x[i]=x[i-1]+1;
    }
    
    pair< pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T> > RES;
    
    
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        RES=make_pair(calculateProfile_cuda_3(motifs,win,nwin,pseq,nseq,strict),x);
//         RES.first=calculateProfile_cuda_3(motifs,win,nwin,pseq,nseq,strict);
    }else{
        //Y=calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);
    }
#else
    //arma::Mat<T> Y = calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);    
#endif



    return move(RES);
}
template pair<pair<vector<gMDense<float>>,vector<gMDense<float>>>,vector<float>> geco::methods::kdmotifs::getMotifsProfile3(const gKDMotifs<float> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict);
template pair<pair<vector<gMDense<double>>,vector<gMDense<double>>>,vector<double>> geco::methods::kdmotifs::getMotifsProfile3(const gKDMotifs<double> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict);

template <typename T> 
pair<pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T>> geco::methods::kdmotifs::getMotifsProfile3(const gKDMotifs<T> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,bool strict){
    geco::bed::gBedData<geco::bed::gBedBaseEntry> pos,neg;
    size_t base=0;
    for(auto r=regions.begin();r!=regions.end();r++){
//         bool isF = r->i_isForward;
        char strand = (r->i_isForward)?('+'):('-');
        if(labels[base]>0){
            pos.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,strand));
        }else{
            neg.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,strand));
        }
        base++;
    }
    vector<string> pseq=pos.getSequences(genomeFile);
    vector<string> nseq=neg.getSequences(genomeFile);
    
    return move(getMotifsProfile3<T>(motifs,pseq,nseq,halfInterval,halfWin,strict));
 
//     unsigned int win = 2*halfWin +1;
//     unsigned int nwin = 2*(halfInterval-halfWin)+1;
//     pair<vector<gMDense<T>>,vector<gMDense<T>>> Y;
//     
// #ifdef GECO_HAS_CUDA
//     if(useCuda()){
//         Y=calculateProfile_cuda_3(motifs,win,nwin,pseq,nseq,strict);
//     }else{
//         //Y=calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);
//     }
// #else
//     //arma::Mat<T> Y = calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);    
// #endif
//         
//     vector<T> x(nwin);
//     x[0]=-(T)halfInterval+(T)halfWin;
//     for(unsigned int i=1;i<x.size();i++){
//         x[i]=x[i-1]+1;
//     }
//     
//     return move(make_pair(Y,x));
}
template pair<pair<vector<gMDense<float>>,vector<gMDense<float>>>,vector<float>> geco::methods::kdmotifs::getMotifsProfile3(const gKDMotifs<float> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, bool strict);
template pair<pair<vector<gMDense<double>>,vector<gMDense<double>>>,vector<double>> geco::methods::kdmotifs::getMotifsProfile3(const gKDMotifs<double> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, bool strict);

template<typename T>
pair<pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T>> geco::methods::kdmotifs::getMotifsProfile4(const gKDMotifs<T> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict){
    unsigned int nmot = motifs.i_ncols;
    unsigned int win = 2*halfWin +1;
    unsigned int nwin = 2*(halfInterval-halfWin)+1;
    vector<T> x(nwin);
    x[0]=-(T)halfInterval+(T)halfWin;
    for(unsigned int i=1;i<x.size();i++){
        x[i]=x[i-1]+1;
    }
    
    pair< pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T> > RES;
    
    
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        RES=make_pair(calculateProfile_cuda_3(motifs,win,nwin,pseq,nseq,strict),x);
//         RES.first=calculateProfile_cuda_3(motifs,win,nwin,pseq,nseq,strict);
    }else{
        //Y=calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);
    }
#else
    //arma::Mat<T> Y = calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);    
#endif
    pair< pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T> > RES2=make_pair(make_pair(vector<gMDense<T>>(nmot,gMDense<T>(pseq.size(),nwin)),vector<gMDense<T>>(nmot,gMDense<T>(nseq.size(),nwin))),x);
    vector<gMDense<T>> & pfeat_o=RES.first.first;    
    vector<gMDense<T>> & pfeat_n=RES2.first.first;
    for(unsigned int m=0;m<nmot;++m){
        for(indexType s=0;s<pseq.size();++s){
            for(indexType p=0;p<nwin;++p){
                //mat.i_values.i_ptr[row + col*mat.i_nrows]
                pfeat_n[m].i_values[s + p*pfeat_n[m].i_nrows]=pfeat_o[s].i_values[p + m*pfeat_o[s].i_nrows];
            }
        }
    }
    
    vector<gMDense<T>> & nfeat_o=RES.first.second;    
    vector<gMDense<T>> & nfeat_n=RES2.first.second;
    for(unsigned int m=0;m<nmot;++m){
        for(indexType s=0;s<nseq.size();++s){
            for(indexType p=0;p<nwin;++p){
                //mat.i_values.i_ptr[row + col*mat.i_nrows]
                nfeat_n[m].i_values[s + p*nfeat_n[m].i_nrows]=nfeat_o[s].i_values[p + m*nfeat_o[s].i_nrows];
            }
        }
    }

    return move(RES2);
}
template pair<pair<vector<gMDense<float>>,vector<gMDense<float>>>,vector<float>> geco::methods::kdmotifs::getMotifsProfile4(const gKDMotifs<float> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict);
template pair<pair<vector<gMDense<double>>,vector<gMDense<double>>>,vector<double>> geco::methods::kdmotifs::getMotifsProfile4(const gKDMotifs<double> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict);

template <typename T> 
pair<pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T>> geco::methods::kdmotifs::getMotifsProfile4(const gKDMotifs<T> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,bool strict){
    geco::bed::gBedData<geco::bed::gBedBaseEntry> pos,neg;
    size_t base=0;
    for(auto r=regions.begin();r!=regions.end();r++){
//         bool isF = r->i_isForward;
        char strand = (r->i_isForward)?('+'):('-');
        if(labels[base]>0){
            pos.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,strand));
        }else{
            neg.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,strand));
        }
        base++;
    }
    vector<string> pseq=pos.getSequences(genomeFile);
    vector<string> nseq=neg.getSequences(genomeFile);
    
    return move(getMotifsProfile4<T>(motifs,pseq,nseq,halfInterval,halfWin,strict));
 
//     unsigned int win = 2*halfWin +1;
//     unsigned int nwin = 2*(halfInterval-halfWin)+1;
//     pair<vector<gMDense<T>>,vector<gMDense<T>>> Y;
//     
// #ifdef GECO_HAS_CUDA
//     if(useCuda()){
//         Y=calculateProfile_cuda_3(motifs,win,nwin,pseq,nseq,strict);
//     }else{
//         //Y=calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);
//     }
// #else
//     //arma::Mat<T> Y = calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);    
// #endif
//         
//     vector<T> x(nwin);
//     x[0]=-(T)halfInterval+(T)halfWin;
//     for(unsigned int i=1;i<x.size();i++){
//         x[i]=x[i-1]+1;
//     }
//     
//     return move(make_pair(Y,x));
}
template pair<pair<vector<gMDense<float>>,vector<gMDense<float>>>,vector<float>> geco::methods::kdmotifs::getMotifsProfile4(const gKDMotifs<float> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, bool strict);
template pair<pair<vector<gMDense<double>>,vector<gMDense<double>>>,vector<double>> geco::methods::kdmotifs::getMotifsProfile4(const gKDMotifs<double> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, bool strict);






template <typename T> 
gMDense<T> geco::methods::kdmotifs::getMotifsProfileParameters(const pair<gMDense<T>,vector<T>> & profile,unsigned int halfInterval, unsigned int halfWin,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval){
    unsigned int win = 2*halfWin +1;
    const arma::Mat<T> Y=armadillo::convert<gMDense,arma::Mat,T>(profile.first);
    arma::Mat<T> A(Y.n_rows,12);
    gMThreadException te;
#pragma omp parallel for
    for(arma::uword r=0;r<Y.n_rows;r++){
        if(te.exceptionOccurred()) continue;
        try{
            A.row(r)=arma::conv_to< arma::Col<T> >::from(getA_internal(vector<double>(profile.second.begin(),profile.second.end()),arma::conv_to< vector<double> >::from(Y.row(r)),halfInterval,win,nrep,maxSelectedPeak,half90interval));
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    
    return move(armadillo::convert<arma::Mat,gMDense,T>(A));
    
}
template gMDense<float> geco::methods::kdmotifs::getMotifsProfileParameters(const pair<gMDense<float>,vector<float>> & profile,unsigned int halfInterval, unsigned int halfWin,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval);
template gMDense<double> geco::methods::kdmotifs::getMotifsProfileParameters(const pair<gMDense<double>,vector<double>> & profile,unsigned int halfInterval, unsigned int halfWin,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval);

template <typename T>
pair<gMDense<T>,pair<gMDense<T>,vector<T>>> geco::methods::kdmotifs::getMotifsProfile(const gKDMotifs<T> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval,bool strict){
    geco::bed::gBedData<geco::bed::gBedBaseEntry> pos,neg;
    size_t base=0;
    for(auto r=regions.begin();r!=regions.end();r++){
        bool isF = r->i_isForward;            
        if(labels[base]>0){
//             cout << base << "\t" << ((isF)?("+"):("-")) << endl;            
            pos.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,isF));
        }else{
            neg.push_back(gBedBaseEntry(r->i_chrom,r->i_start+centers[base]-halfInterval,r->i_start+centers[base]+halfInterval+1,r->i_name,r->i_score,isF));
        }
        base++;
    }
    vector<string> pseq=pos.getSequences(genomeFile);
    vector<string> nseq=neg.getSequences(genomeFile);
//     cout << pseq[0] << endl;
//     cout << nseq[0] << endl;
    
    unsigned int win = 2*halfWin +1;
    unsigned int nwin = 2*(halfInterval-halfWin)+1;

#ifdef GECO_HAS_CUDA
    arma::Mat<T> Y;    
    if(useCuda()){
        Y=calculateProfile_cuda(motifs,win,nwin,pseq,nseq,strict);
    }else{
        Y=calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);
    }
#else
    arma::Mat<T> Y = calculateProfile_arma(motifs,win,nwin,pseq,nseq,strict);    
#endif
        
    vector<T> x(Y.n_cols);
    x[0]=-(T)halfInterval+(T)halfWin;
    for(unsigned int i=1;i<x.size();i++){
        x[i]=x[i-1]+1;
    }
    
    arma::Mat<T> A(Y.n_rows,12);
    
    gMThreadException te;
#pragma omp parallel for
    for(arma::uword r=0;r<Y.n_rows;r++){
        if(te.exceptionOccurred()) continue;
        try{
            A.row(r)=arma::conv_to< arma::Col<T> >::from(getA_internal(vector<double>(x.begin(),x.end()),arma::conv_to< vector<double> >::from(Y.row(r)),halfInterval,win,nrep,maxSelectedPeak,half90interval));
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();

    return make_pair(armadillo::convert<arma::Mat,gMDense,T>(A),make_pair(armadillo::convert<arma::Mat,gMDense,T>(Y),x));

/*
    pair<gMDense<T>,vector<T>> profile=getMotifsProfile2(motifs,regions,centers,labels,genomeFile,halfInterval,halfWin,strict);
    gMDense<T> A=getMotifsProfileParameters(profile,halfInterval,halfWin,nrep,maxSelectedPeak,half90interval);
    return make_pair(A,profile);
    */
}
template pair<gMDense<float>,pair<gMDense<float>,vector<float>>> geco::methods::kdmotifs::getMotifsProfile(const gKDMotifs<float> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval,bool strict);
template pair<gMDense<double>,pair<gMDense<double>,vector<double>>> geco::methods::kdmotifs::getMotifsProfile(const gKDMotifs<double> & motifs,const gBedData<gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval,bool strict);

//FINE


template <typename entryTypeSelected,typename entryTypeExcluded>
class regionFinder{
public:    
    const geco::bed::gBedData<entryTypeSelected> & i_selected;
    std::multiset<geco::bed::gBedBaseEntry,geco::bed::comparePositionA> i_excluded;
    const geco::sequence_retrievers::gLocal2bitSequenceRetriever & i_retriever;
    
    struct keyCalc{ 
        
    public:
        unsigned int i_gcbins;
        unsigned int i_rpbins;
        
        keyCalc(unsigned int gcbins,unsigned int rpbins):i_gcbins(gcbins),i_rpbins(rpbins){
        }
        
        inline unsigned int operator() (const std::string & s){
            unsigned int len=s.length();
            unsigned int gc_count=0;
            unsigned int rp_count=0;
            
            for(auto c:s){
                rp_count+=islower(c)!=0;
                char tc=tolower(c);
                if(tc=='c' || tc=='g'){
                    gc_count++;
                }
            }
            unsigned int gc_bin = (gc_count * i_gcbins) / len;
            unsigned int rp_bin = (rp_count * i_rpbins) / len;
            return gc_bin * (i_rpbins+1) + rp_bin;
        }
    } getKey;
    
    std::unordered_map<gSize,std::unordered_map<unsigned int,std::list<entryTypeSelected>>> lendmap;
    
    std::vector<gSize> valid_chroms;
    size_t nchroms;
    
    std::discrete_distribution<gSize> ldistr;
    std::vector<gSize> lcounts;
    std::vector<gSize> lvals;
    bool goOn;

    unsigned int searchSemaphore;
    bool doNotSearch;

    void updateLdistr(){
        gSize TOT=0;
        for(auto dmap=lendmap.begin();dmap!=lendmap.end();dmap++){
            for(auto dco=dmap->second.begin();dco!=dmap->second.end();dco++){
                lcounts[TOT]+=dco->second.size();
            }
            lvals[TOT]=dmap->first;
            //std::cout << dmap->first << " (" << TOT << "): " << lcounts[TOT] << std::endl;
            TOT++;
        }
        
        ldistr=std::discrete_distribution<gSize>(lcounts.begin(),lcounts.end());
    }

    regionFinder(const geco::bed::gBedData<entryTypeSelected> & selected,const geco::bed::gBedData<entryTypeExcluded> & excluded,const geco::sequence_retrievers::gLocal2bitSequenceRetriever & retriever,unsigned int gcbins,unsigned int rpbins):i_selected(selected),i_excluded(geco::bed::comparePositionA(10)),i_retriever(retriever),getKey(gcbins,rpbins){
        for(auto r:excluded){
            i_excluded.insert(r);
        }

        //std::cout << "Building structures...";std::cout.flush();
        
        std::set<std::string> tchroms;
        for(auto r:selected){
            tchroms.insert(r.i_chrom);
            gSize len=r.i_end-r.i_start;
            lendmap[len][getKey(r.getSequence(retriever))].push_back(r);
        }
        for(auto s:tchroms){
            for(gSize i=0;i<retriever.getReferenceCount();i++){
                if(retriever.getReferenceName(i)==s){
                    valid_chroms.push_back(i);
                }
            }
        }
        //std::cout << lendmap.size() << " maps found" << std::endl;
        nchroms=valid_chroms.size();
        
        lcounts=std::vector<gSize>(lendmap.size(),0);
        lvals=std::vector<gSize>(lendmap.size(),0);

        updateLdistr();
    }

    inline geco::bed::gBedBaseEntry searchSuitableRegion(std::mt19937_64 & rng,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & currentset,unsigned int npeaks){
        while(doNotSearch){};
#pragma omp critical
        searchSemaphore++;

        geco::bed::gBedBaseEntry reg;

        if(currentset.size()<npeaks){
            gSize len;
            gSize gchrom;
            gSize rlen;
            bool keynot;
            do{
                len=lvals[ldistr(rng)];
                do{
                    gchrom=valid_chroms[rng() % nchroms];
                    rlen=i_retriever.getReferenceLength(gchrom);
                }while(rlen<=len);
                gPos start=rng() % (rlen-len);
                reg=geco::bed::gBedBaseEntry(i_retriever.getReferenceName(gchrom),start,start+len,".",0,'+');
                keynot=false;
                if(currentset.size()<npeaks){
                    auto ptr=lendmap[len].find(getKey(reg.getSequence(i_retriever)));
                    if(ptr==lendmap[len].end() || ptr->second.size()==0){
                        keynot=true;
                    }else if(i_excluded.find(reg)!=i_excluded.end()){
                        keynot=true;
                    }
                }
                if(goOn && omp_get_thread_num()==0){
                    try{
                        checkInterrupt(0);
                    }catch(gMethodsInterrupt & e){
                        goOn=false;
                        warning(string("regionFinder: user interruption before end"));
                    }
                }
            }while(keynot && goOn);
        }
        
#pragma omp critical
        searchSemaphore--;
        
        if(!goOn){
            return geco::bed::gBedBaseEntry();
        }else{
            return reg;
        }
    }

    inline void update(const geco::bed::gBedBaseEntry & reg,geco::bed::gBedData<geco::bed::gBedBaseEntry> & currentset,geco::bed::gBedData<entryTypeSelected> & currentsel,unsigned int npeaks){
        if(reg.i_end-reg.i_start>0){
#pragma omp critical        
            {
                while(searchSemaphore>0){}
                if(currentset.size()<npeaks){
                    doNotSearch=true;
                    auto key=getKey(reg.getSequence(i_retriever));
                    gSize len=reg.i_end-reg.i_start;
                    auto lmap=lendmap.find(len);
                    auto entry=lmap->second.find(key);
                    if(entry->second.size()>0){
                        if(i_excluded.find(reg)==i_excluded.end()){
                            currentset.push_back(reg);
                            i_excluded.insert(reg);
                            currentsel.push_back(*entry->second.begin());
                            entry->second.pop_front();
                            auto id=find(lvals.begin(),lvals.end(),len);
                            if(id!=lvals.end()){
                                lcounts[distance(lvals.begin(),id)]--;
                                ldistr=std::discrete_distribution<gSize>(lcounts.begin(),lcounts.end());
                            }
                            //std::cout << std::setw(5) << currentset.size() << "\tSEQ:" << std::setw(10) << len << std::setw(10) << key << std::setw(10) << lendmap.size() << std::endl;
                        }//else std::cout << "\tOverlappa con una già presente" << std::endl;
                    }//else std::cout << "\tChiave e lunghezza sono già esaurite" << std::endl;
                }
                doNotSearch=false;
            }
        }
    }

    inline std::pair<geco::bed::gBedData<entryTypeSelected>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> operator() (unsigned int npeaks){
        geco::bed::gBedData<geco::bed::gBedBaseEntry> BKG;
        geco::bed::gBedData<entryTypeSelected> SEL;
        if(npeaks==0){
            npeaks=i_selected.size();
        }else if(npeaks>i_selected.size()){
            throw gMethodException("Requested number of background regions greater than number of selected regions");
        }
            
        doNotSearch=false;
        searchSemaphore=0;
        gMethodsOut msgs(2);
        msgs.Start([=]{std::cout << "\tprogress...";});
        gMThreadException te;
        unsigned p=0;

        int tn=omp_get_thread_num();
        goOn=true;
#pragma omp parallel
        try{
            int tID;
            std::mt19937_64 rng;
#pragma omp critical
            {
                tID=omp_get_thread_num();
//                 cout << "tID: " << tID << " seed:" << getgMehtodsSeed()+tID << endl;
                rng=std::mt19937_64(getgMehtodsSeed()+tID);
            }
            geco::bed::gBedBaseEntry reg;
            do{
                if(goOn){
                    reg=searchSuitableRegion(rng,BKG,npeaks);
                    if(goOn){
                        update(reg,BKG,SEL,npeaks);
//                         if(goOn && tID==tn){
//                             try{
//                                 checkInterrupt(tn);
//                             }catch(gMethodsInterrupt & e){
//                                 goOn=false;
//                                 warning(string("regionFinder: user interruption before end"));
//                             }
//                         }
                        if(goOn){
#pragma omp critical
                            {
                                unsigned np = (BKG.size() * 100) / npeaks ;
                                if(np != p){
                                    p=np;
                                    msgs.End([=]{std::cout << "\t" << p << "%";}); 
                                    msgs.Start([=]{std::cout << "\tprogress...";});
                                }
                            }
                        }
                    }
                }
            }while(BKG.size()<npeaks && goOn);
            
//             do{
//                 reg=searchSuitableRegion(rng,BKG,npeaks);
//                 update(reg,BKG,SEL,npeaks);
//                 if(omp_get_thread_num()==tn){
//                     checkInterrupt(tn);
//                 }
// #pragma omp critical
//                 {
//                     unsigned np = (BKG.size() * 100) / npeaks ;
//                     if(np != p){
//                         p=np;
//                         msgs.End([=]{std::cout << "\t" << p << "%";}); 
//                         msgs.Start([=]{std::cout << "\tprogress...";});
//                     }
//                 }
//             }while(BKG.size()<npeaks && goOn);
            
        }catch(...){
            te.CaptureException();
#pragma omp critical    
            {
                goOn=false;
                std::cout << std::endl;
            }
        }
        te.Rethrow();
        msgs.End([=]{std::cout << "\tcompleted";});
        return std::make_pair(SEL,BKG);
    }
};


            
template <typename entryTypeSelected,typename entryTypeExcluded>
std::pair<geco::bed::gBedData<entryTypeSelected>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> geco::methods::kdmotifs::getBackground(const geco::bed::gBedData<entryTypeSelected> & selected,const geco::bed::gBedData<entryTypeExcluded> & excluded,unsigned int nregions,unsigned int ngcbins,unsigned int nrpbins,const std::string & fileName2bits){
    if(nregions>selected.size()){
        throw gMethodException("Requested number of background regions greater than number of selected regions");
    }    
    gMethodsOut msgs(1);
    msgs.Start([=]{std::cout << "\tInitializing...";});
    geco::sequence_retrievers::gLocal2bitSequenceRetriever retriever(fileName2bits,true);
    regionFinder<entryTypeSelected,entryTypeExcluded> finder(selected,excluded,retriever,ngcbins,nrpbins);
    msgs.End([=]{});
    msgs.Start([=]{std::cout << "\tSearching background regions...";});
    auto res = finder(nregions);
    msgs.End([=]{});
    return move(res);
}
template std::pair<geco::bed::gBedData<geco::bed::gBedBaseEntry>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> geco::methods::kdmotifs::getBackground(const geco::bed::gBedData<geco::bed::gBedBaseEntry> & selected,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & excluded,unsigned int nregions,unsigned int ngcbins,unsigned int nrpbins,const std::string & fileName2bits);
template std::pair<geco::bed::gBedData<geco::bed::gBedBaseEntry>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> geco::methods::kdmotifs::getBackground(const geco::bed::gBedData<geco::bed::gBedBaseEntry> & selected,const geco::bed::gBedData<geco::bed::gNarrowPeak> & excluded,unsigned int nregions,unsigned int ngcbins,unsigned int nrpbins,const std::string & fileName2bits);
template std::pair<geco::bed::gBedData<geco::bed::gNarrowPeak>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> geco::methods::kdmotifs::getBackground(const geco::bed::gBedData<geco::bed::gNarrowPeak> & selected,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & excluded,unsigned int nregions,unsigned int ngcbins,unsigned int nrpbins,const std::string & fileName2bits);
template std::pair<geco::bed::gBedData<geco::bed::gNarrowPeak>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> geco::methods::kdmotifs::getBackground(const geco::bed::gBedData<geco::bed::gNarrowPeak> & selected,const geco::bed::gBedData<geco::bed::gNarrowPeak> & excluded,unsigned int nregions,unsigned int ngcbins,unsigned int nrpbins,const std::string & fileName2bits);







void checkForcompatibility(const gKDMInfo & i1, const gKDMInfo & i2){
    if(
            ( (i1.i_doubleStrand!=i2.i_doubleStrand) || (i1.i_l!=i2.i_l) ) ||
            ( (i1.i_estimated && i2.i_estimated) && (i2.i_k!=i1.i_k) )
    ){
        throw gMethodException("Incompatible kmer distributions");
    }
}

template<typename T> 
gKDistr<T> geco::methods::kdmotifs::toKDistr(const gKDMotifs<T> & kmotifs){
    return std::move( gKDistr<T>( toSparse(kmotifs),kmotifs) );
}
template gKDistr<float> geco::methods::kdmotifs::toKDistr(const gKDMotifs<float> & kmotifs);
template gKDistr<double> geco::methods::kdmotifs::toKDistr(const gKDMotifs<double> & kmotifs);


template<typename T> 
gKDMotifs<T> geco::methods::kdmotifs::toKDMotifs(const gKDistr<T> & kdistr){
    return std::move( gKDMotifs<T>( toDense(kdistr),kdistr) );
}
template gKDMotifs<float> geco::methods::kdmotifs::toKDMotifs(const gKDistr<float> & kdistr);
template gKDMotifs<double> geco::methods::kdmotifs::toKDMotifs(const gKDistr<double> & kdistr);


template<typename T>
gMDense<T> doCentroids(const gMDense<indexType> CLUS,const gMSparse<T> & counts,indexType rank){

    arma::SpMat<T> C=armadillo::convert<gMSparse,arma::SpMat,T>(counts);
    arma::Mat<indexType> cl=armadillo::convert<gMDense,arma::Mat,indexType>(CLUS);
    arma::Mat<T> RES(C.n_rows,rank);
    gMThreadException te;
#pragma omp parallel for
    for(arma::uword i=0;i<rank;i++){
        try{
            if(te.exceptionOccurred()) continue;
            arma::uvec sel=arma::find(cl==i);
            RES.col(i)=arma::sqrt(arma::mean(arma::Mat<T>(square(C.cols(sel))), 1 ));
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return std::move(armadillo::convert<arma::Mat,gMDense,T>(RES));
}

template<typename T>
gMDense<T> doCentroids(const gMDense<indexType> CLUS,const gMDense<T> & counts,indexType rank){
    arma::Mat<T> C=armadillo::convert<gMDense,arma::Mat,T>(counts);
    arma::Mat<indexType> cl=armadillo::convert<gMDense,arma::Mat,indexType>(CLUS);
    arma::Mat<T> RES(C.n_rows,rank);
    gMThreadException te;
#pragma omp parallel for
    for(arma::uword i=0;i<rank;i++){
        if(te.exceptionOccurred()) continue;
        try{
            arma::uvec sel=arma::find(cl==i);
            RES.col(i)=arma::sqrt(arma::mean(arma::Mat<T>(square(C.cols(sel))), 1 ));
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return std::move(armadillo::convert<arma::Mat,gMDense,T>(RES));
}

template <template<typename> class D,typename T>
bool kcountsNMFInit<D,T>::operator ()(const D<T> & V, nmf::gNMF<D,T> & M, indexType rank) const{
    gMethodsOut outers(1);
    outers.Start([=]{cout << "Cluster W0 initialization";});
    gMethodsOut msgs(2);
    if(i_bhattacharyya){
        msgs.Start([=]{cout << "\tCalculating counts distances (Bhattacharyya)...";});
    }else{
        msgs.Start([=]{cout << "\tCalculating counts distances (Hellinger)...";});
    }
    //D<T> tD=multiply(V,V,true,false);
    
    gMDense<T> tD=kdmDistance(V,V,i_bhattacharyya);
    msgs.End([=]{});
    
    
    msgs.Start([=]{cout << "\tClustering (" << getClusteringMethodName(i_method) << ")...";});
    gHCluster dend(tD,i_method);
    msgs.End([=]{});    
    
    msgs.Start([=]{cout << "\tBuilding centroids...";});
    gMDense<indexType> cl=dend.cutree_k(rank);
    M.i_W=doCentroids(cl,V,rank);
    msgs.End([=]{});    
    outers.End([=]{});
    return true;
}
//explicit specializations
template class geco::methods::kdmotifs::kcountsNMFInit<gMDense,float>;
template class geco::methods::kdmotifs::kcountsNMFInit<gMDense,double>;
template class geco::methods::kdmotifs::kcountsNMFInit<gMSparse,float>;
template class geco::methods::kdmotifs::kcountsNMFInit<gMSparse,double>;



template <template<typename> class D,typename T>
bool kcountsNMFInit2<D,T>::operator ()(const dataFeeder<D,T> & feeder, nmf::gNMF2<D,T> & M, indexType rank) const{
    gMethodsOut outers(1);
    outers.Start([=]{cout << "Cluster W0 initialization";});
    gMethodsOut msgs(2);
    D<T> V=feeder(0);
    if(i_bhattacharyya){
        msgs.Start([=]{cout << "\tCalculating counts distances (Bhattacharyya)...";});
    }else{
        msgs.Start([=]{cout << "\tCalculating counts distances (Hellinger)...";});
    }
    //D<T> tD=multiply(V,V,true,false);
    
    gMDense<T> tD=kdmDistance(V,V,i_bhattacharyya);
    msgs.End([=]{});
    
    
    msgs.Start([=]{cout << "\tClustering (" << getClusteringMethodName(i_method) << ")...";});
    gHCluster dend(tD,i_method);
    msgs.End([=]{});    
    
    msgs.Start([=]{cout << "\tBuilding centroids...";});
    gMDense<indexType> cl=dend.cutree_k(rank);
    M.i_W=doCentroids(cl,V,rank);
    msgs.End([=]{});    
    outers.End([=]{});
    return true;
}
//explicit specializations
template class geco::methods::kdmotifs::kcountsNMFInit2<gMDense,float>;
template class geco::methods::kdmotifs::kcountsNMFInit2<gMDense,double>;
template class geco::methods::kdmotifs::kcountsNMFInit2<gMSparse,float>;
template class geco::methods::kdmotifs::kcountsNMFInit2<gMSparse,double>;

template<typename T>
gKDistr<T> geco::methods::kdmotifs::kdmDistr(const std::vector<std::string> & sequences,indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict){
    gMSparse<T> mat;
    if(estimatekmers && ngaps!=0){
        std::cout << "1" << std::endl;
        gLmerEstimator<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
        std::cout << "2" << std::endl;
        mat=counter.countSparse(sequences,strict);
    }else if(ngaps==0){
        gLmerCounter<T> counter(kmerlength,doubleStrand);
        mat=counter.countSparse(sequences,strict);
    }else{
        gGappedKmerCounter<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
        mat=counter.countSparse(sequences,strict);
    }
    std::cout << "3" << std::endl;
    return gKDistr<T>(toHellinger(mat),kmerlength,ngaps,doubleStrand,estimatekmers);    
}
template gKDistr<float> geco::methods::kdmotifs::kdmDistr(const std::vector<std::string> & sequences,indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict);
template gKDistr<double> geco::methods::kdmotifs::kdmDistr(const std::vector<std::string> & sequences,indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict);

/*
template<typename T>
gKDistr<T> geco::methods::kdmotifs::kdmDistrProfile_old(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict){
    size_t slen=sequence.size();
    if(slen<win ) throw gMethodException("kdmDistrProfile: win cannot be longer than sequence");
    size_t nwins=slen-win+1;
    vector<string> sseqs(nwins);
#pragma omp parallel for
    for(size_t pos=0;pos<nwins;pos++){
        sseqs[pos]=sequence.substr(pos,win);
    }
    return kdmDistr<T>(sseqs,kmerlength,ngaps,doubleStrand,estimatekmers,strict);
    
//     gMSparse<T> mat;
//     
//     if(estimatekmers && ngaps!=0){
//         gLmerEstimator<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
//         size_t nwins=sequence.size()-win+1;
//         vector<string> sseqs(nwins);
// #pragma omp parallel for
//         for(size_t pos=0;pos<nwins;pos++){
//             sseqs[pos]=sequence.substr(pos,win);
//         }
//         mat=counter.countSparse(sseqs,strict);
//     }else if(ngaps==0){
//         gLmerCounter<T> counter(kmerlength,doubleStrand);
//         size_t nwins=sequence.size()-win+1;
//         vector<string> sseqs(nwins);
// #pragma omp parallel for
//         for(size_t pos=0;pos<nwins;pos++){
//             sseqs[pos]=sequence.substr(pos,win);
//         }
//         mat=counter.countSparse(sseqs,strict);
//     }else{
//         gGappedKmerCounter<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
//         mat=counter.countWinSparse(sequence,win,strict);
//     }
//     return gKDistr<T>(toHellinger(mat),kmerlength,ngaps,doubleStrand,estimatekmers);    
    
}
template gKDistr<float> geco::methods::kdmotifs::kdmDistrProfile(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict);
template gKDistr<double> geco::methods::kdmotifs::kdmDistrProfile(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict);
*/

template<typename T>
gKDistr<T> geco::methods::kdmotifs::kdmDistrProfile(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict){
    gMSparse<T> mat;
    if(estimatekmers && ngaps!=0){
        throw gMethodException("estimatekmers not yet implemented");
    }else if(ngaps==0){
        gLmerCounter<T> counter(kmerlength,doubleStrand);
        mat=counter.countWinSparse(sequence,win,strict);
    }else{
        gGappedKmerCounter<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
        mat=counter.countWinSparse(sequence,win,strict);
    }
    return gKDistr<T>(toHellinger(mat),kmerlength,ngaps,doubleStrand,estimatekmers);
}
template gKDistr<float> geco::methods::kdmotifs::kdmDistrProfile(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict);
template gKDistr<double> geco::methods::kdmotifs::kdmDistrProfile(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers,bool strict);


template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMotifs(const gKDistr<T> & kdistr,indexType rank, indexType np, initMethod imethod,const gMDense<T> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, T tolerance,indexType maxIterations){
    gNMF<gMSparse,T> factorization;
    if(initOnly){
        switch(imethod){
            case random:
                factorization=NMF((const gMSparse<T> &) kdistr,randomNMFInit<gMSparse,T>(true,false),factorization,rank,tolerance,maxIterations);
                break;
            case nndsvd:
                factorization=NMF((const gMSparse<T> &) kdistr,nndsvdNMFInit<gMSparse,T>(true,false),factorization,rank,tolerance,maxIterations);
                break;
            case cluster:
                factorization=NMF((const gMSparse<T> &) kdistr,kcountsNMFInit<gMSparse,T>(clustering_method,bhattacharyya),factorization,rank,tolerance,maxIterations);
                break;
            case custom:
                factorization=NMF((const gMSparse<T> &) kdistr,customNMFInit<gMSparse,T>(W0,gMDense<T>()),factorization,rank,tolerance,maxIterations);
                break;
        }
    }else{
        if(np==0){
            if(rank==0){
                switch(imethod){
                    case random:
                        factorization=NMF((const gMSparse<T> &) kdistr,randomNMFInit<gMSparse,T>(true,false),PNMF<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                    case nndsvd:
                        factorization=NMF((const gMSparse<T> &) kdistr,nndsvdNMFInit<gMSparse,T>(true,false),PNMF<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                    case cluster:
                        factorization=NMF((const gMSparse<T> &) kdistr,kcountsNMFInit<gMSparse,T>(clustering_method,bhattacharyya),PNMF<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                    case custom:
                        factorization=NMF((const gMSparse<T> &) kdistr,customNMFInit<gMSparse,T>(W0,gMDense<T>()),PNMF<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                }
            }else{
                switch(imethod){
                    case random:
                        factorization=NMF((const gMSparse<T> &) kdistr,randomNMFInit<gMSparse,T>(true,false),PNMF<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                    case nndsvd:
                        factorization=NMF((const gMSparse<T> &) kdistr,nndsvdNMFInit<gMSparse,T>(true,false),PNMF<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                    case cluster:
                        factorization=NMF((const gMSparse<T> &) kdistr,kcountsNMFInit<gMSparse,T>(clustering_method,bhattacharyya),PNMF<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                    case custom:
                        factorization=NMF((const gMSparse<T> &) kdistr,customNMFInit<gMSparse,T>(W0,gMDense<T>()),PNMF<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                }
            }
        }else{
            switch(imethod){
                case random:
                    factorization=NMF((const gMSparse<T> &) kdistr,randomNMFInit<gMSparse,T>(true,false),OPNMF<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
                case nndsvd:
                    factorization=NMF((const gMSparse<T> &) kdistr,nndsvdNMFInit<gMSparse,T>(true,false),OPNMF<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
                case cluster:
                    factorization=NMF((const gMSparse<T> &) kdistr,kcountsNMFInit<gMSparse,T>(clustering_method,bhattacharyya),OPNMF<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
                case custom:
                    factorization=NMF((const gMSparse<T> &) kdistr,customNMFInit<gMSparse,T>(W0,gMDense<T>()),OPNMF<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
            }
        }
    }
    return std::move(gKDMotifs<T>(factorization.i_W,kdistr));
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMotifs(const gKDistr<float> & kdistr,indexType rank, indexType np, initMethod imethod,const gMDense<float> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, float tolerance,indexType maxIterations);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMotifs(const gKDistr<double> & kdistr,indexType rank, indexType np, initMethod imethod,const gMDense<double> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, double tolerance,indexType maxIterations);

template<typename T>
class myFeeder{
public:
    vector<string> i_sequences;
    gKDMInfo i_info;
    myFeeder(const vector<string> & sequences,const gKDMInfo & info):i_sequences(sequences),i_info(info){
        
    }
    
    gMSparse<T> & operator()(indexType & p){
        if(p==0){
            return kdmDistr<T>(i_sequences,i_info.i_l,i_info.i_k,i_info.i_doubleStrand,i_info.i_estimated,false);
        }else{
            return kdmDistr<T>(i_sequences,i_info.i_l,i_info.i_k,i_info.i_doubleStrand,i_info.i_estimated,false);
        }
    }
};    


template<template<typename> class D, typename T>
D<T> myback(const boost::any & data,indexType p){
    pair<vector<string>,gKDMInfo> cdata=boost::any_cast<pair<vector<string>,gKDMInfo>>(data);
    const vector<string> & sequences = cdata.first;
    const gKDMInfo & info = cdata.second;
    if(p==0){
        return kdmDistr<T>(sequences,info.i_l,info.i_k,info.i_doubleStrand,info.i_estimated,false);
    }else{
        vector<string> samp(p);
        std::sample(sequences.begin(), sequences.end(), samp.begin(),p,std::mt19937{std::random_device{}()});
        return kdmDistr<T>(samp,info.i_l,info.i_k,info.i_doubleStrand,info.i_estimated,false);
    }
}


template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMotifs2(const vector<string> & sequences,const gKDMInfo & info,indexType rank, indexType np, initMethod imethod,const gMDense<T> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, T tolerance,indexType maxIterations){
    dataFeeder<gMSparse,T> feeder(make_pair(sequences,info),myback<gMSparse,T>);
    gNMF2<gMSparse,T> factorization;
    
    if(initOnly){
        gMSparse<T> kdistr=kdmDistr<T>(sequences,info.i_l,info.i_k,info.i_doubleStrand,info.i_estimated,false);
        switch(imethod){
            case random:
                factorization=NMF2(feeder,randomNMFInit2<gMSparse,T>(true,false),factorization,rank,tolerance,maxIterations);
                break;
            case nndsvd:
                factorization=NMF2(feeder,nndsvdNMFInit2<gMSparse,T>(true,false),factorization,rank,tolerance,maxIterations);
                break;
            case cluster:
                factorization=NMF2(feeder,kcountsNMFInit2<gMSparse,T>(clustering_method,bhattacharyya),factorization,rank,tolerance,maxIterations);
                break;
            case custom:
                factorization=NMF2(feeder,customNMFInit2<gMSparse,T>(W0,gMDense<T>()),factorization,rank,tolerance,maxIterations);
                break;
        }
    }else{
        if(np==0){
            gMSparse<T> kdistr=kdmDistr<T>(sequences,info.i_l,info.i_k,info.i_doubleStrand,info.i_estimated,false);
            if(rank==0){
                switch(imethod){
                    case random:
                        factorization=NMF2(feeder,randomNMFInit2<gMSparse,T>(true,false),PNMF2<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                    case nndsvd:
                        factorization=NMF2(feeder,nndsvdNMFInit2<gMSparse,T>(true,false),PNMF2<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                    case cluster:
                        factorization=NMF2(feeder,kcountsNMFInit2<gMSparse,T>(clustering_method,bhattacharyya),PNMF2<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                    case custom:
                        factorization=NMF2(feeder,customNMFInit2<gMSparse,T>(W0,gMDense<T>()),PNMF2<gMSparse,T>(true,true),kdistr.i_ncols *0.5 ,tolerance,maxIterations);
                        break;
                }
            }else{
                switch(imethod){
                    case random:
                        factorization=NMF2(feeder,randomNMFInit2<gMSparse,T>(true,false),PNMF2<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                    case nndsvd:
                        factorization=NMF2(feeder,nndsvdNMFInit2<gMSparse,T>(true,false),PNMF2<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                    case cluster:
                        factorization=NMF2(feeder,kcountsNMFInit2<gMSparse,T>(clustering_method,bhattacharyya),PNMF2<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                    case custom:
                        factorization=NMF2(feeder,customNMFInit2<gMSparse,T>(W0,gMDense<T>()),PNMF2<gMSparse,T>(false,true),rank ,tolerance,maxIterations);
                        break;
                }
            }
        }else{
            switch(imethod){
                case random:
                    factorization=NMF2(feeder,randomNMFInit2<gMSparse,T>(true,false),OPNMF2<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
                case nndsvd:
                    factorization=NMF2(feeder,nndsvdNMFInit2<gMSparse,T>(true,false),OPNMF2<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
                case cluster:
                    factorization=NMF2(feeder,kcountsNMFInit2<gMSparse,T>(clustering_method,bhattacharyya),OPNMF2<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
                case custom:
                    factorization=NMF2(feeder,customNMFInit2<gMSparse,T>(W0,gMDense<T>()),OPNMF2<gMSparse,T>(np,true),rank,tolerance,maxIterations);
                    break;
            }
        }
    }
    return move(gKDMotifs<T>(factorization.i_W,info));
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMotifs2(const vector<string> & sequences,const gKDMInfo & info,indexType rank, indexType np, initMethod imethod,const gMDense<float> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, float tolerance,indexType maxIterations);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMotifs2(const vector<string> & sequences,const gKDMInfo & info,indexType rank, indexType np, initMethod imethod,const gMDense<double> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, double tolerance,indexType maxIterations);


template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMotifs(const gKDMotifs<T> & motifs,indexType rank, indexType np, initMethod imethod,const gMDense<T> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, T tolerance,indexType maxIterations){
    return kdmMotifs(toKDistr(motifs),rank,np,imethod,W0,clustering_method,bhattacharyya,initOnly,tolerance,maxIterations);
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMotifs(const gKDMotifs<float> & motifs,indexType rank, indexType np, initMethod imethod,const gMDense<float> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, float tolerance,indexType maxIterations);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMotifs(const gKDMotifs<double> & motifs,indexType rank, indexType np, initMethod imethod,const gMDense<double> & W0,cluster::gHClusterMethod clustering_method, bool bhattacharyya, bool initOnly, double tolerance,indexType maxIterations);

template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmEstimateMotifs(const gKDMotifs<T> & motifs){
    gLmerEstimator<T> counter(motifs.i_l,motifs.i_l-motifs.i_k,motifs.i_doubleStrand);
    gMDense<T> mat=counter.countDense(motifs);
    return std::move(gKDMotifs<T>(toHellinger(mat),motifs.i_l,motifs.i_k,motifs.i_doubleStrand,true));
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmEstimateMotifs(const gKDMotifs<float> & motifs);
template gKDMotifs<double> geco::methods::kdmotifs::kdmEstimateMotifs(const gKDMotifs<double> & motifs);


template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMotifsFromPWM(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand,indexType padSize){
    if(ngaps==0){
        gLmerCounter<T> counter(kmerlength,doubleStrand);
        gMDense<T> mat=counter.countPWM(PWMs,padSize);
        return std::move(gKDMotifs<T>(toHellinger(mat),kmerlength,0,doubleStrand,false));
    }else{
        gGappedKmerCounter<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
        gMDense<T> mat=counter.countPWM(PWMs,padSize);
        return std::move(gKDMotifs<T>(toHellinger(mat),kmerlength,ngaps,doubleStrand,false));
    }
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMotifsFromPWM(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand,indexType padSize);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMotifsFromPWM(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand,indexType padSize);

template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMotifsFromPWM2(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand){
    if(ngaps==0){
        gLmerCounter<T> counter(kmerlength,doubleStrand);
        gMDense<T> mat=counter.countPWM2(PWMs);
        return std::move(gKDMotifs<T>(toHellinger(mat),kmerlength,0,doubleStrand,false));
    }else{
        gGappedKmerCounter<T> counter(kmerlength,kmerlength-ngaps,doubleStrand);
        gMDense<T> mat=counter.countPWM2(PWMs);
        return std::move(gKDMotifs<T>(toHellinger(mat),kmerlength,ngaps,doubleStrand,false));
    }
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMotifsFromPWM2(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMotifsFromPWM2(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand);


template<typename T>
vector<T> geco::methods::kdmotifs::kdmGetPWMThresholds(const geco::kmers::gPWMSet & PWMs,T alpha,T initialGranularity,bool forcedGranularity,T maxGranularity){
    vector<T> thresholds(PWMs.size());
    gMethodsOut inner(2);
    gMThreadException te;
#pragma omp parallel for
    for(gSize i=0;i<PWMs.size();i++){
        try{
            if(te.exceptionOccurred()) continue;
            inner.Start([=]{std::cout << i+1 << "\t" << PWMs.getName(i) << " ...";});
            thresholds[i]=PWMs[i].computeThreshold(alpha,initialGranularity,forcedGranularity,maxGranularity);
            inner.End([=]{ std::cout << thresholds[i];});
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return std::move(thresholds);
}
template vector<float> geco::methods::kdmotifs::kdmGetPWMThresholds(const geco::kmers::gPWMSet & PWMs,float alpha,float initialGranularity,bool forcedGranularity,float maxGranularity);
template vector<double> geco::methods::kdmotifs::kdmGetPWMThresholds(const geco::kmers::gPWMSet & PWMs,double alpha,double initialGranularity,bool forcedGranularity,double maxGranularity);



template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMotifsReduce(const gKDMotifs<T>& kmotifs,indexType rank, T tolerance,indexType maxIterations){
//     gNMF<gMSparse,T> factorization;
//     gMSparse<T> M=toSparse(kmotifs);
//     if(rank==0){
//         factorization=NMF(M,kcountsNMFInit<gMSparse,T>(METHOD_METR_WARD_D2,true),PNMF<gMSparse,T>(true,true),kmotifs.i_ncols * 0.8,tolerance,maxIterations);
//     }else{
//         factorization=NMF(M,kcountsNMFInit<gMSparse,T>(METHOD_METR_WARD_D2,true),PNMF<gMSparse,T>(false,true),rank,tolerance,maxIterations);
//     }
//     return move(gKDMotifs<T>(factorization.i_W,kmotifs));
    bool i_bhattacharyya=false;
    gHClusterMethod i_method=METHOD_METR_WARD_D2;
    gMethodsOut msgs(2);
    if(i_bhattacharyya){
        msgs.Start([=]{cout << "\tCalculating counts distances (Bhattacharyya)...";});
    }else{
        msgs.Start([=]{cout << "\tCalculating counts distances (Hellinger)...";});
    }
    
    gMDense<T> tD=kdmDistance(kmotifs,kmotifs,i_bhattacharyya);
    msgs.End([=]{});
    
    
    msgs.Start([=]{cout << "\tClustering (" << getClusteringMethodName(i_method) << ")...";});
    gHCluster dend(tD,i_method);
    msgs.End([=]{});    
    
    msgs.Start([=]{cout << "\tBuilding centroids...";});
    gMDense<indexType> cl=dend.cutree_k(rank);
    msgs.End([=]{});    
    return std::move(gKDMotifs<T>(doCentroids(cl,kmotifs,rank),kmotifs));
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMotifsReduce(const gKDMotifs<float> & kmotifs,indexType rank, float tolerance,indexType maxIterations);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMotifsReduce(const gKDMotifs<double>& kmotifs,indexType rank, double tolerance,indexType maxIterations);

    
template<typename T>
gKDMFeatures<T> geco::methods::kdmotifs::kdmFeatures(const gKDMotifs<T> & kmotifs,const gKDistr<T> & kdistr){
    //here we should check for compatibility between kmotifs and kdistr.
//     if(kmotifs.i_estimated){
//         if( (kdistr.i_l==kmotifs.i_l & kdistr.i_estimated)  | (
//     }else
    return gKDMFeatures<T>(multiply((const gMSparse<T> &) kdistr,(const gMDense<T> &)kmotifs,true,false),kmotifs);
}
template gKDMFeatures<float> geco::methods::kdmotifs::kdmFeatures(const gKDMotifs<float> & kmotifs,const gKDistr<float> & kdistr);
template gKDMFeatures<double> geco::methods::kdmotifs::kdmFeatures(const gKDMotifs<double> & kmotifs,const gKDistr<double> & kdistr);


template<typename T>
gKDMFeatures<T> geco::methods::kdmotifs::kdmFeatures(const gKDMotifs<T> & kmotifs,const vector<std::string> & sequences, bool strict){
    gKDistr<T> counts=kdmDistr<T>(sequences,kmotifs.i_l,kmotifs.i_k,kmotifs.i_doubleStrand,kmotifs.i_estimated,strict);
    return gKDMFeatures<T>(multiply((const gMSparse<T> &)counts,(const gMDense<T> &)kmotifs,true,false),kmotifs);
}
template gKDMFeatures<float> geco::methods::kdmotifs::kdmFeatures(const gKDMotifs<float> & kmotifs,const vector<string> & sequences, bool strict);
template gKDMFeatures<double> geco::methods::kdmotifs::kdmFeatures(const gKDMotifs<double> & kmotifs,const vector<string> & sequences, bool strict);

template<typename T>
gKDMFeatures<T> geco::methods::kdmotifs::kdmFeaturesProfile(const gKDMotifs<T> & kmotifs,const string & sequence, indexType win, bool strict){
    return kdmFeatures(kmotifs,kdmDistrProfile<T>(sequence,win,kmotifs.i_l,kmotifs.i_k,kmotifs.i_doubleStrand,kmotifs.i_estimated,strict));
}
template gKDMFeatures<float> geco::methods::kdmotifs::kdmFeaturesProfile(const gKDMotifs<float> & kmotifs,const string & sequence, indexType win, bool strict);
template gKDMFeatures<double> geco::methods::kdmotifs::kdmFeaturesProfile(const gKDMotifs<double> & kmotifs,const string & sequence, indexType win, bool strict);



template<typename T>
gMDense<T> geco::methods::kdmotifs::predict(const gKDModel<T> & model,const gKDMFeatures<T> & features){
    arma::Mat<T> tmp1=armadillo::convert<gMDense,arma::Mat,T>(multiply((const gMDense<T> &) features,(const gMDense<T> &)model.i_coefficients,false,false));
    tmp1.each_row()+=arma::Row<T>(model.i_intercept);
    tmp1.transform([](T val) {return exp(val)/(1+exp(val));} );
    
    return std::move(armadillo::convert<arma::Mat,gMDense,T>(tmp1));
}
template gMDense<float> geco::methods::kdmotifs::predict(const gKDModel<float> & model,const gKDMFeatures<float> & features);
template gMDense<double> geco::methods::kdmotifs::predict(const gKDModel<double> & model,const gKDMFeatures<double> & features);

template<typename T>
gMDense<T> geco::methods::kdmotifs::predict(const gKDModel<T> & model,const gKDistr<T> & kdistr){
    return predict(model,kdmFeatures(model,kdistr));
    
}
template gMDense<float> geco::methods::kdmotifs::predict(const gKDModel<float> & model,const gKDistr<float> & kdistr);
template gMDense<double> geco::methods::kdmotifs::predict(const gKDModel<double> & model,const gKDistr<double> & kdistr);


template<typename T>
gMDense<T> geco::methods::kdmotifs::predict(const gKDModel<T> & model,const vector<string> & sequences,bool strict){
    return predict(model,kdmDistr<T>(sequences,model.i_l,model.i_k,model.i_doubleStrand,model.i_estimated,strict));
}
template gMDense<float> geco::methods::kdmotifs::predict(const gKDModel<float> & model,const vector<string> & sequences,bool strict);
template gMDense<double> geco::methods::kdmotifs::predict(const gKDModel<double> & model,const vector<string> & sequences,bool strict);



template<typename T,typename entryType>
gMDense<T> geco::methods::kdmotifs::predict(const gKDModel<T> & model,const geco::bed::gBedData<entryType> & bedData,const string & genomeFile,bool strict){
    return predict(model,bedData.getSequences( genomeFile ),strict);
}
template gMDense<float> geco::methods::kdmotifs::predict(const gKDModel<float> & model,const bed::gBedData<bed::gBedBaseEntry> & bedData,const string & fileName,bool strict);
template gMDense<float> geco::methods::kdmotifs::predict(const gKDModel<float> & model,const bed::gBedData<bed::gNarrowPeak> & bedData,const string & fileName,bool strict);
template gMDense<double> geco::methods::kdmotifs::predict(const gKDModel<double> & model,const bed::gBedData<bed::gBedBaseEntry> & bedData,const string & fileName,bool strict);
template gMDense<double> geco::methods::kdmotifs::predict(const gKDModel<double> & model,const bed::gBedData<bed::gNarrowPeak> & bedData,const string & fileName,bool strict);


            







template< template<typename> class D, typename T>
gMDense<T> geco::methods::kdmotifs::kdmSimilarity(const D<T> & D1,const D<T> & D2){
//     if(std::is_same<D<T>,gKDistr<T>>::value){
//         return toDense(multiply((const gMSparse<T> &)D1,(const gMSparse<T> &)D2,true,false));
//     }else if(std::is_same<D<T>,gKDMotifs<T>>::value){
//         return multiply((const gMDense<T> &)D1,(const gMDense<T> &)D2,true,false);
//     }
    if(std::is_convertible<D<T>,gMSparse<T>>::value){
        return toDense(multiply((const gMSparse<T> &)D1,(const gMSparse<T> &)D2,true,false));
    }else if(std::is_convertible<D<T>,gMDense<T>>::value){
        return multiply((const gMDense<T> &)D1,(const gMDense<T> &)D2,true,false);
    }
    
}
template gMDense<float> geco::methods::kdmotifs::kdmSimilarity(const gKDistr<float> & D1,const gKDistr<float> & D2);
template gMDense<double> geco::methods::kdmotifs::kdmSimilarity(const gKDistr<double> & D1,const gKDistr<double> & D2);
template gMDense<float> geco::methods::kdmotifs::kdmSimilarity(const gKDMotifs<float> & D1,const gKDMotifs<float> & D2);
template gMDense<double> geco::methods::kdmotifs::kdmSimilarity(const gKDMotifs<double> & D1,const gKDMotifs<double> & D2);


template< template<typename> class D, typename T>
gMDense<T> geco::methods::kdmotifs::kdmDistance(const D<T> & D1,const D<T> & D2,bool bhattacharyya){
    gMDense<T> res;
    arma::Mat<T> S=armadillo::convert<gMDense,arma::Mat,T>(kdmSimilarity(D1,D2));
    if(bhattacharyya){ //bhattacharyya distance -log(BC)
        arma::uvec sel1=arma::find(S==0);
        arma::uvec sel2=arma::find(S>0);
        T mm=arma::min(S(sel2));
        S.elem(sel1)+=mm;
        res=armadillo::convert<arma::Mat,gMDense,T>( -arma::log(S));
        
    }else{ //hellinger distance sqrt(1-BC)
        arma::Mat<T> B= 1-S;
        arma::uvec sel=arma::find(S>1);
        B.elem(sel)*=0;
        res=armadillo::convert<arma::Mat,gMDense,T>( arma::sqrt(B));
    }
    return std::move(res);
}
template gMDense<float> geco::methods::kdmotifs::kdmDistance(const gKDistr<float> & D1,const gKDistr<float> & D2,bool bhattacharyya);
template gMDense<double> geco::methods::kdmotifs::kdmDistance(const gKDistr<double> & D1,const gKDistr<double> & D2,bool bhattacharyya);
template gMDense<float> geco::methods::kdmotifs::kdmDistance(const gKDMotifs<float> & D1,const gKDMotifs<float> & D2,bool bhattacharyya);
template gMDense<double> geco::methods::kdmotifs::kdmDistance(const gKDMotifs<double> & D1,const gKDMotifs<double> & D2,bool bhattacharyya);



#ifdef GECO_HAS_CUDA
namespace geco{
    namespace methods{
        namespace cuda{
            template <typename T>
            cudaDense<T> & isGreater(const cudaDense<T> & A, cudaDense<T> & B, cudaDense<T> & C,indexType nrep);
        }
    }
};
#endif


template<typename T>
std::pair<gMDense<T>,gMDense<T>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<T> & models1,const gKDModel<T> & models2,indexType nrep){
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        cuda::autoCudaHandles cudaHandles;
        cuda::cudaDense<T> S12(models1.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> S11(models1.i_ncols,models1.i_ncols);
        cuda::cudaDense<T> S22(models2.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> C1 = cuda::convert<gMDense,cuda::cudaDense,T>(models1.i_coefficients,false,cudaHandles);
        cuda::cudaDense<T> C2 = cuda::convert<gMDense,cuda::cudaDense,T>(models2.i_coefficients,false,cudaHandles);
        {
            cuda::cudaDense<T> M1 = cuda::convert<gMDense,cuda::cudaDense,T>(models1,false,cudaHandles);
            cuda::cudaDense<T> M2 = cuda::convert<gMDense,cuda::cudaDense,T>(models2,false,cudaHandles);
            cuda::multiply_cuda<T>(M1,M2,S12,true,false,1,0,cudaHandles);
            cuda::multiply_cuda<T>(M1,M1,S11,true,false,1,0,cudaHandles);
            cuda::multiply_cuda<T>(M2,M2,S22,true,false,1,0,cudaHandles);
        }
        cuda::cudaDense<T> PARTIAL12(models1.i_coefficients.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> N(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        //cuda::multiply_cuda<T>(C1,S12,PARTIAL12,true,false,1,0,cudaHandles);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C1,S12,PARTIAL12,true,false,1,0,cudaHandles),C2,N,false,false,1,0,cudaHandles);
        cuda::cudaDense<T> PARTIAL11(models1.i_coefficients.i_ncols,models1.i_ncols);
        cuda::cudaDense<T> diag1(models1.i_coefficients.i_ncols,1);
        cuda::cudaDense<T> D1(models1.i_coefficients.i_ncols,models1.i_coefficients.i_ncols);
        //cuda::multiply_cuda<T>(C1,S11,PARTIAL11,true,false,1,0,cudaHandles);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C1,S11,PARTIAL11,true,false,1,0,cudaHandles),C1,D1,false,false,1,0,cudaHandles);
        cuda::sqrt(cuda::get_diag<T>(D1,diag1),diag1);
        cuda::cudaDense<T> PARTIAL22(models2.i_coefficients.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> diag2(models2.i_coefficients.i_ncols,1);
        cuda::cudaDense<T> D2(models2.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        //cuda::multiply_cuda<T>(C2,S22,PARTIAL22,true,false,1,0,cudaHandles);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C2,S22,PARTIAL22,true,false,1,0,cudaHandles),C2,D2,false,false,1,0,cudaHandles);
        cuda::sqrt(cuda::get_diag<T>(D2,diag2),diag2);
        cuda::cudaDense<T> D(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        cuda::multiply_cuda<T>(diag1,diag2,D,false,true,1,0,cudaHandles);
        cuda::cudaDense<T> CC(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        cuda::elementwise_divide(N,D,CC);
        cuda::cudaDense<T> PP(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols,0);
        if(nrep>0){
            cuda::cudaDense<T> tC1(C1.i_nrows,C1.i_ncols);
            gMethodsOut inner(2);    
            inner.Start([=]{std::cout << "\titerations [" << (1) << "-" << (10) << "]...";});
            for(indexType i = 0; i<nrep;i++){
                tC1 = select_rows(C1,randperm(C1.i_nrows,C1.i_nrows));
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC1,S12,PARTIAL12,true,false,1,0,cudaHandles),C2,N,false,false,1,0,cudaHandles);
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC1,S11,PARTIAL11,true,false,1,0,cudaHandles),tC1,D1,false,false,1,0,cudaHandles);
                cuda::sqrt(cuda::get_diag<T>(D1,diag1),diag1);
                
                cuda::multiply_cuda<T>(diag1,diag2,D,false,true,1,0,cudaHandles);
                cuda::elementwise_divide(N,D,N);
                cuda::isGreater(CC,N,PP,nrep);
                
                checkInterrupt();                
                if((i+1) % 10 ==0){
                    inner.End([=]{});
                    inner.Start([=]{std::cout << "\titerations [" << (i+2) << "-" << (i+11) << "]...";});
                }
            }
            inner.End([=]{});
        }
        arma::Mat<T> CC0=armadillo::convert<gMDense,arma::Mat,T>(cuda::convert<cuda::cudaDense,gMDense,T>(PP,false,cudaHandles));
        CC0 /= nrep;
        return move(make_pair(move(cuda::convert<cuda::cudaDense,gMDense,T>(CC,false,cudaHandles)),move(armadillo::convert<arma::Mat,gMDense,T>(CC0))));
        //return move(make_pair(move(cuda::convert<cuda::cudaDense,gMDense,T>(CC,false,cudaHandles)),move(cuda::convert<cuda::cudaDense,gMDense,T>(PP,false,cudaHandles))));
    }
#endif
    arma::Mat<T> S12,S11,S22;
    {
        arma::Mat<T> M1=armadillo::convert<gMDense,arma::Mat,T>(models1);
        arma::Mat<T> M2=armadillo::convert<gMDense,arma::Mat,T>(models2);
        
        S12 =  arma::trans(M1) * M2;
        S11 =  arma::trans(M1) * M1;
        S22 =  arma::trans(M2) * M2;
    }
    arma::Mat<T> C1 = armadillo::convert<gMDense,arma::Mat,T>(models1.i_coefficients);
    arma::Mat<T> C2 = armadillo::convert<gMDense,arma::Mat,T>(models2.i_coefficients);
    
//     arma::Mat<T> S12_C2 = S12 * C2;
//     arma::Mat<T> S11_C1 = S11 * C1;
//     arma::Mat<T> S22_C2 = S22 * C2;
    
    arma::Mat<T>  N = arma::trans(C1) * S12 * C2;
    arma::Mat<T>  tD1 = arma::trans(C1) * S11 * C1;
    arma::Mat<T>  tD2 = arma::trans(C2) * S22 * C2;
    
    
    arma::Col<T> D1 = arma::sqrt(tD1.diag());
    arma::Row<T> D2 = trans(arma::sqrt(tD2.diag()));
   
    arma::Mat<T> CC = N / (D1 * D2);
    arma::Mat<T> CC0(CC.n_rows,CC.n_cols,arma::fill::zeros);
    
    if(nrep>0){
        gMethodsOut inner(2);    
        inner.Start([=]{std::cout << "\titerations [" << (1) << "-" << (10) << "]...";});
        for(indexType i = 0; i<nrep;i++){
            arma::Mat<T> tC1 = C1.rows(arma::randperm(C1.n_rows));
            tD1 = arma::trans(tC1) * S11 * tC1;
            D1 = arma::sqrt(tD1.diag());
            
            arma::Mat<T> NCC = (arma::trans(tC1) * S12 * C2)/(D1 * D2);
            
            arma::uvec sel = arma::find(CC>0);
            CC0.elem(sel) = CC0.elem(sel) + (NCC.elem(sel) >= CC.elem(sel));
            sel = arma::find(CC<=0);
            CC0.elem(sel) = CC0.elem(sel) + (NCC.elem(sel) <= CC.elem(sel));
            
            checkInterrupt();                
            if((i+1) % 10 ==0){
                inner.End([=]{});
                inner.Start([=]{std::cout << "\titerations [" << (i+2) << "-" << (i+11) << "]...";});
            }
        }
        inner.End([=]{});
  }
  
  return move(make_pair(move(armadillo::convert<arma::Mat,gMDense,T>(CC)),move(armadillo::convert<arma::Mat,gMDense,T>(CC0/nrep))));
}
template std::pair<gMDense<float>,gMDense<float>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<float> & Models1,const gKDModel<float> & Models2,indexType nrep);
template std::pair<gMDense<double>,gMDense<double>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<double> & Models1,const gKDModel<double> & Models2,indexType nrep);



template<typename T>
std::pair<gMDense<T>,gMDense<T>> geco::methods::kdmotifs::kdmModelsCorrelation2(const gKDModel<T> & models1,const gKDModel<T> & models2,indexType nrep){
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        cuda::autoCudaHandles cudaHandles;
        cuda::cudaDense<T> S12(models1.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> S11(models1.i_ncols,models1.i_ncols);
        cuda::cudaDense<T> S22(models2.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> C1 = cuda::convert<gMDense,cuda::cudaDense,T>(models1.i_coefficients,false,cudaHandles);
        cuda::cudaDense<T> C2 = cuda::convert<gMDense,cuda::cudaDense,T>(models2.i_coefficients,false,cudaHandles);
        {
            cuda::cudaDense<T> M1 = cuda::convert<gMDense,cuda::cudaDense,T>(models1,false,cudaHandles);
            cuda::cudaDense<T> M2 = cuda::convert<gMDense,cuda::cudaDense,T>(models2,false,cudaHandles);
            cuda::multiply_cuda<T>(M1,M2,S12,true,false,1,0,cudaHandles);
            cuda::multiply_cuda<T>(M1,M1,S11,true,false,1,0,cudaHandles);
            cuda::multiply_cuda<T>(M2,M2,S22,true,false,1,0,cudaHandles);
        }
        cuda::cudaDense<T> PARTIAL12(models1.i_coefficients.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> N(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        //cuda::multiply_cuda<T>(C1,S12,PARTIAL12,true,false,1,0,cudaHandles);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C1,S12,PARTIAL12,true,false,1,0,cudaHandles),C2,N,false,false,1,0,cudaHandles);
        cuda::cudaDense<T> PARTIAL11(models1.i_coefficients.i_ncols,models1.i_ncols);
        cuda::cudaDense<T> diag1(models1.i_coefficients.i_ncols,1);
        cuda::cudaDense<T> D1(models1.i_coefficients.i_ncols,models1.i_coefficients.i_ncols);
        //cuda::multiply_cuda<T>(C1,S11,PARTIAL11,true,false,1,0,cudaHandles);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C1,S11,PARTIAL11,true,false,1,0,cudaHandles),C1,D1,false,false,1,0,cudaHandles);
        cuda::sqrt(cuda::get_diag<T>(D1,diag1),diag1);
        cuda::cudaDense<T> PARTIAL22(models2.i_coefficients.i_ncols,models2.i_ncols);
        cuda::cudaDense<T> diag2(models2.i_coefficients.i_ncols,1);
        cuda::cudaDense<T> D2(models2.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        //cuda::multiply_cuda<T>(C2,S22,PARTIAL22,true,false,1,0,cudaHandles);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C2,S22,PARTIAL22,true,false,1,0,cudaHandles),C2,D2,false,false,1,0,cudaHandles);
        cuda::sqrt(cuda::get_diag<T>(D2,diag2),diag2);
        cuda::cudaDense<T> D(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        cuda::multiply_cuda<T>(diag1,diag2,D,false,true,1,0,cudaHandles);
        cuda::cudaDense<T> CC(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols);
        cuda::elementwise_divide(N,D,CC);
        cuda::cudaDense<T> PP(models1.i_coefficients.i_ncols,models2.i_coefficients.i_ncols,0);
        if(nrep>0){
            cout << "doing correlation2 cuda" << endl;
            cuda::cudaDense<T> tC1(C1.i_nrows,C1.i_ncols);
            cuda::cudaDense<T> tC2(C2.i_nrows,C2.i_ncols);
            gMethodsOut inner(2);    
            inner.Start([=]{std::cout << "\titerations [" << (1) << "-" << (10) << "]...";});
            for(indexType i = 0; i<nrep;i++){
                tC1 = select_rows(C1,randperm(C1.i_nrows,C1.i_nrows));
//                 cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC1,S12,PARTIAL12,true,false,1,0,cudaHandles),C2,N,false,false,1,0,cudaHandles);
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC1,S11,PARTIAL11,true,false,1,0,cudaHandles),tC1,D1,false,false,1,0,cudaHandles);
                cuda::sqrt(cuda::get_diag<T>(D1,diag1),diag1);
                
                tC2 = select_rows(C2,randperm(C2.i_nrows,C2.i_nrows));
//                 cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC1,S12,PARTIAL12,true,false,1,0,cudaHandles),C2,N,false,false,1,0,cudaHandles);
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC2,S22,PARTIAL22,true,false,1,0,cudaHandles),tC2,D2,false,false,1,0,cudaHandles);
                cuda::sqrt(cuda::get_diag<T>(D2,diag2),diag2);
                
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC1,S12,PARTIAL12,true,false,1,0,cudaHandles),tC2,N,false,false,1,0,cudaHandles);
                
                cuda::multiply_cuda<T>(diag1,diag2,D,false,true,1,0,cudaHandles);
                cuda::elementwise_divide(N,D,N);
                cuda::isGreater(CC,N,PP,nrep);
                
                checkInterrupt();                
                if((i+1) % 10 ==0){
                    inner.End([=]{});
                    inner.Start([=]{std::cout << "\titerations [" << (i+2) << "-" << (i+11) << "]...";});
                }
            }
            inner.End([=]{});
        }
        arma::Mat<T> CC0=armadillo::convert<gMDense,arma::Mat,T>(cuda::convert<cuda::cudaDense,gMDense,T>(PP,false,cudaHandles));
        CC0 /= nrep;
        return move(make_pair(move(cuda::convert<cuda::cudaDense,gMDense,T>(CC,false,cudaHandles)),move(armadillo::convert<arma::Mat,gMDense,T>(CC0))));
        //return move(make_pair(move(cuda::convert<cuda::cudaDense,gMDense,T>(CC,false,cudaHandles)),move(cuda::convert<cuda::cudaDense,gMDense,T>(PP,false,cudaHandles))));
    }
#endif
    arma::Mat<T> S12,S11,S22;
    {
        arma::Mat<T> M1=armadillo::convert<gMDense,arma::Mat,T>(models1);
        arma::Mat<T> M2=armadillo::convert<gMDense,arma::Mat,T>(models2);
        
        S12 =  arma::trans(M1) * M2;
        S11 =  arma::trans(M1) * M1;
        S22 =  arma::trans(M2) * M2;
    }
    arma::Mat<T> C1 = armadillo::convert<gMDense,arma::Mat,T>(models1.i_coefficients);
    arma::Mat<T> C2 = armadillo::convert<gMDense,arma::Mat,T>(models2.i_coefficients);
    
    
    arma::Mat<T>  N = arma::trans(C1) * S12 * C2;
    arma::Mat<T>  tD1 = arma::trans(C1) * S11 * C1;
    arma::Mat<T>  tD2 = arma::trans(C2) * S22 * C2;
    
    
    arma::Col<T> D1 = arma::sqrt(tD1.diag());
    arma::Row<T> D2 = trans(arma::sqrt(tD2.diag()));
   
    arma::Mat<T> CC = N / (D1 * D2);
    arma::Mat<T> CC0(CC.n_rows,CC.n_cols,arma::fill::zeros);
    
    if(nrep>0){
        cout << "doing correlation2 arma" << endl;
        gMethodsOut inner(2);    
        inner.Start([=]{std::cout << "\titerations [" << (1) << "-" << (10) << "]...";});
        for(indexType i = 0; i<nrep;i++){
            arma::uvec perm1=arma::randperm(C1.n_rows);
            arma::Mat<T> tC1 = C1.rows(perm1);
            tD1 = arma::trans(tC1) * S11 * tC1;
            D1 = arma::sqrt(tD1.diag());
            
            arma::uvec perm2=arma::randperm(C2.n_rows);
            arma::Mat<T> tC2 = C2.rows(perm2);            
            tD2 = arma::trans(tC2) * S22 * tC2;
            D2 = arma::sqrt(tD2.diag());
            
            arma::Mat<T> NCC = (arma::trans(tC1) * S12 * tC2)/(D1 * D2);
            
            arma::uvec sel = arma::find(CC>0);
            CC0.elem(sel) = CC0.elem(sel) + (NCC.elem(sel) >= CC.elem(sel));
            sel = arma::find(CC<=0);
            CC0.elem(sel) = CC0.elem(sel) + (NCC.elem(sel) <= CC.elem(sel));
            
            checkInterrupt();                
            if((i+1) % 10 ==0){
                inner.End([=]{});
                inner.Start([=]{std::cout << "\titerations [" << (i+2) << "-" << (i+11) << "]...";});
            }
        }
        inner.End([=]{});
  }
  
  return move(make_pair(move(armadillo::convert<arma::Mat,gMDense,T>(CC)),move(armadillo::convert<arma::Mat,gMDense,T>(CC0/nrep))));
}
template std::pair<gMDense<float>,gMDense<float>> geco::methods::kdmotifs::kdmModelsCorrelation2(const gKDModel<float> & Models1,const gKDModel<float> & Models2,indexType nrep);
template std::pair<gMDense<double>,gMDense<double>> geco::methods::kdmotifs::kdmModelsCorrelation2(const gKDModel<double> & Models1,const gKDModel<double> & Models2,indexType nrep);





template<typename T>            
std::pair<gMDense<T>,gMDense<T>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<T> & models,indexType nrep){
    //return kdmModelsCorrelation(models,models,nrep);
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        cuda::autoCudaHandles cudaHandles;
        cuda::cudaDense<T> S(models.i_ncols,models.i_ncols);
        cuda::cudaDense<T> C = cuda::convert<gMDense,cuda::cudaDense,T>(models.i_coefficients,false,cudaHandles);

        {
            cuda::cudaDense<T> M = cuda::convert<gMDense,cuda::cudaDense,T>(models,false,cudaHandles);
            cuda::multiply_cuda<T>(M,M,S,true,false,1,0,cudaHandles);
        }
        
        cuda::cudaDense<T> PARTIAL(models.i_coefficients.i_ncols,models.i_ncols);
        cuda::cudaDense<T> N(models.i_coefficients.i_ncols,models.i_coefficients.i_ncols);
        cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(C,S,PARTIAL,true,false,1,0,cudaHandles),C,N,false,false,1,0,cudaHandles);
        
        cuda::cudaDense<T> diag(models.i_coefficients.i_ncols,1);
        cuda::sqrt(cuda::get_diag<T>(N,diag),diag);
        cuda::cudaDense<T> diago=diag;

        cuda::cudaDense<T> D(models.i_coefficients.i_ncols,models.i_coefficients.i_ncols);
        cuda::multiply_cuda<T>(diag,diag,D,false,true,1,0,cudaHandles);
        
        cuda::cudaDense<T> CC(models.i_coefficients.i_ncols,models.i_coefficients.i_ncols);
        cuda::elementwise_divide(N,D,CC);
        
        cuda::cudaDense<T> PP(models.i_coefficients.i_ncols,models.i_coefficients.i_ncols,0);
        if(nrep>0){
            cuda::cudaDense<T> tC(C.i_nrows,C.i_ncols);
            gMethodsOut inner(2);    
            inner.Start([=]{std::cout << "\titerations [" << (1) << "-" << (10) << "]...";});
            for(indexType i = 0; i<nrep/2;i++){
                tC = select_rows(C,randperm(C.i_nrows,C.i_nrows));
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC,S,PARTIAL,true,false,1,0,cudaHandles),C,N,false,false,1,0,cudaHandles);
                cuda::multiply_cuda<T>(cuda::multiply_cuda<T>(tC,S,PARTIAL,true,false,1,0,cudaHandles),tC,D,false,false,1,0,cudaHandles);
                cuda::sqrt(cuda::get_diag<T>(D,diag),diag);
                
                cuda::multiply_cuda<T>(diag,diago,D,false,true,1,0,cudaHandles);
                cuda::elementwise_divide(N,D,N);
                cuda::isGreater(CC,N,PP,nrep);
                
                checkInterrupt();                
                if((i+1) % 10 ==0){
                    inner.End([=]{});
                    inner.Start([=]{std::cout << "\titerations [" << (i+2) << "-" << (i+11) << "]...";});
                }
            }
            inner.End([=]{});
        }
        arma::Mat<T> CC0=armadillo::convert<gMDense,arma::Mat,T>(cuda::convert<cuda::cudaDense,gMDense,T>(PP,false,cudaHandles));
        CC0 = (CC0 + trans(CC0)) / nrep;
        return move(make_pair(move(cuda::convert<cuda::cudaDense,gMDense,T>(CC,false,cudaHandles)),move(armadillo::convert<arma::Mat,gMDense,T>(CC0))));
    }
#endif
    arma::Mat<T> S;
    {
        arma::Mat<T> M=armadillo::convert<gMDense,arma::Mat,T>(models);

        
        S =  arma::trans(M) * M;
    }
    arma::Mat<T> C = armadillo::convert<gMDense,arma::Mat,T>(models.i_coefficients);
    
    arma::Mat<T>  N = arma::trans(C) * S * C;
    arma::Mat<T>  tD = arma::trans(C) * S * C;
    arma::Col<T> D1 = arma::sqrt(N.diag());
    arma::Row<T> D2 = trans(arma::sqrt(N.diag()));
   
    arma::Mat<T> CC = N / (D1 * D2);
    arma::Mat<T> CC0(CC.n_rows,CC.n_cols,arma::fill::zeros);
    
    if(nrep>0){
        gMethodsOut inner(2);    
        inner.Start([=]{std::cout << "\titerations [" << (1) << "-" << (10) << "]...";});
        for(indexType i = 0; i<nrep/2; i++){
            arma::Mat<T> tC = C.rows(arma::randperm(C.n_rows));
            N = arma::trans(tC) * S * C;
            tD = arma::trans(tC) * S * tC;
            D1 = arma::sqrt(tD.diag());
            
            arma::Mat<T> NCC = N / (D1 * D2);
            
            arma::uvec sel = arma::find(CC>0);
            CC0.elem(sel) = CC0.elem(sel) + (NCC.elem(sel) >= CC.elem(sel));
            sel = arma::find(CC<=0);
            CC0.elem(sel) = CC0.elem(sel) + (NCC.elem(sel) <= CC.elem(sel));
            
            checkInterrupt();                
            if((i+1) % 10 ==0){
                inner.End([=]{});
                inner.Start([=]{std::cout << "\titerations [" << (i+2) << "-" << (i+11) << "]...";});
            }
        }
        inner.End([=]{});
  }
  CC0 = (CC0 + trans(CC0))/nrep;
  return move(make_pair(move(armadillo::convert<arma::Mat,gMDense,T>(CC)),move(armadillo::convert<arma::Mat,gMDense,T>(CC0))));    
    
}
template std::pair<gMDense<float>,gMDense<float>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<float> & models,indexType nrep);
template std::pair<gMDense<double>,gMDense<double>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<double> & models,indexType nrep);



template<typename T>            
std::pair<gMDense<T>,gMDense<T>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<T> & models,const gKDMotifs<T> & motifs,indexType nrep){
    gKDModel<T> models2(motifs,armadillo::convert<arma::Mat,gMDense,T>(arma::Mat<T>(motifs.i_ncols,motifs.i_ncols,arma::fill::eye)),vector<T>(motifs.i_ncols,0));
    return kdmModelsCorrelation(models,models2,nrep);
}
template std::pair<gMDense<float>,gMDense<float>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<float> & models,const gKDMotifs<float> & motifs,indexType nrep);
template std::pair<gMDense<double>,gMDense<double>> geco::methods::kdmotifs::kdmModelsCorrelation(const gKDModel<double> & models,const gKDMotifs<double> & motifs,indexType nrep);



template<typename T>
gKDistr<T> geco::methods::kdmotifs::kdmSelectKDistr(const gKDistr<T> & d,const std::vector<indexType> & selection){
    return gKDistr<T>(armadillo::convert<arma::SpMat,gMSparse,T>(armadillo::convert<gMSparse,arma::SpMat,T>(d).cols(arma::uvec(selection))),d);

}
template gKDistr<float> geco::methods::kdmotifs::kdmSelectKDistr(const gKDistr<float> & d,const std::vector<indexType> & selection);
template gKDistr<double> geco::methods::kdmotifs::kdmSelectKDistr(const gKDistr<double> & d,const std::vector<indexType> & selection);

template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmSelectKDMotifs(const gKDMotifs<T> & m,const std::vector<indexType> & selection){
    return gKDMotifs<T>(armadillo::convert<arma::Mat,gMDense,T>(armadillo::convert<gMDense,arma::Mat,T>(m).cols(arma::uvec(selection))),m);    
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmSelectKDMotifs(const gKDMotifs<float> & m,const std::vector<indexType> & selection);
template gKDMotifs<double> geco::methods::kdmotifs::kdmSelectKDMotifs(const gKDMotifs<double> & m,const std::vector<indexType> & selection);



template<typename T>
gKDistr<T> geco::methods::kdmotifs::kdmMergeKDistr(const gKDistr<T> & d1,const gKDistr<T> & d2){
    checkForcompatibility(d1,d2);
    return gKDistr<T>(armadillo::convert<arma::SpMat,gMSparse,T>(arma::join_rows(armadillo::convert<gMSparse,arma::SpMat,T>(d1),armadillo::convert<gMSparse,arma::SpMat,T>(d2))),d1);
}
template gKDistr<float> geco::methods::kdmotifs::kdmMergeKDistr(const gKDistr<float> & d1,const gKDistr<float> & d2);
template gKDistr<double> geco::methods::kdmotifs::kdmMergeKDistr(const gKDistr<double> & d1,const gKDistr<double> & d2);


template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMergeKDMotifs(const gKDMotifs<T> & d1,const gKDMotifs<T> & d2){
    checkForcompatibility(d1,d2);
    return gKDMotifs<T>(armadillo::convert<arma::Mat,gMDense,T>(arma::join_rows(armadillo::convert<gMDense,arma::Mat,T>(d1),armadillo::convert<gMDense,arma::Mat,T>(d2))),d1);
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMergeKDMotifs(const gKDMotifs<float> & d1,const gKDMotifs<float> & d2);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMergeKDMotifs(const gKDMotifs<double> & d1,const gKDMotifs<double> & d2);


template<typename T>
gKDMotifs<T> geco::methods::kdmotifs::kdmMergeKDMotifs(const gKDMotifs<T> & d1,const gMDense<indexType> & labels,indexType rank){
    return move(gKDMotifs<T>(doCentroids(labels,d1,rank),d1));
}
template gKDMotifs<float> geco::methods::kdmotifs::kdmMergeKDMotifs(const gKDMotifs<float> & d1,const gMDense<indexType> & labels,indexType rank);
template gKDMotifs<double> geco::methods::kdmotifs::kdmMergeKDMotifs(const gKDMotifs<double> & d1,const gMDense<indexType> & labels,indexType rank);



template<typename T>
pair<gMDense<T>,pair<gMDense<indexType>,gMDense<T>>> geco::methods::kdmotifs::kdmSilhouette(const gMDense<T> & dist,const gMDense<indexType> & clusterLabels){
    return geco::methods::cluster::getSilhouette(dist,clusterLabels);
}
template pair<gMDense<float>,pair<gMDense<indexType>,gMDense<float>>> geco::methods::kdmotifs::kdmSilhouette(const gMDense<float> & dist,const gMDense<indexType> & clusterLabels);
template pair<gMDense<double>,pair<gMDense<indexType>,gMDense<double>>> geco::methods::kdmotifs::kdmSilhouette(const gMDense<double> & dist,const gMDense<indexType> & clusterLabels);

template<typename T>
pair<gMDense<T>,pair<gMDense<T>,gMDense<indexType>>> geco::methods::kdmotifs::kdmGetW0(const gMSparse<T> & counts,cluster::gHClusterMethod clustering_method,bool goodclusters, indexType minRank, indexType maxRank,bool bhattacharyya){
    arma::SpMat<T> C=armadillo::convert<gMSparse,arma::SpMat,T>(counts);
    gMDense<T> D=kdmDistance(counts,counts,bhattacharyya);
    gHCluster dend(D,clustering_method);
    arma::Mat<indexType> ngood(D.i_nrows,1,arma::fill::zeros);
    arma::Mat<T> silh(D.i_nrows,1,arma::fill::zeros);
    
#pragma omp parallel for
    for(indexType k=1;k<D.i_nrows;k++){
        std::pair<gMDense<T>,std::pair<gMDense<indexType>,gMDense<T>>> res=geco::methods::cluster::getSilhouette(D,dend.cutree_k(k+1));
        arma::Col<T> sc=armadillo::convert<gMDense,arma::Mat,T>(res.second.second).col(0);
        if(goodclusters){
            arma::uvec pp=find(sc>0);
            ngood(k) = pp.n_elem;
        }else{
            ngood(k) = k+1;
        }
        silh(k) = arma::mean(sc);        
    }
    
    arma::uvec s=find(ngood>=minRank);
    indexType h=s(arma::index_max(silh(s)))+1;
    indexType nc=min(h,maxRank);

    gMDense<indexType> mcl=dend.cutree_k(nc);
    arma::Mat<indexType> cl=armadillo::convert<gMDense,arma::Mat,indexType>(mcl);
    arma::Mat<T> ss=armadillo::convert<gMDense,arma::Mat,T>(geco::methods::cluster::getSilhouette(D,mcl).first);
    arma::Col<T> mclus(nc,arma::fill::zeros);
#pragma omp parallel for    
    for(arma::uword i=0;i<nc;i++){
        mclus(i)=arma::mean(ss(arma::find(cl==i)));
    }
    arma::uvec goods;
    if(goodclusters){
        goods=find(mclus>0);
    }else{
        goods=arma::linspace<arma::uvec>(0, nc-1,nc);
    }
    indexType arank=goods.n_elem;
//     cout << "last=" << goods(arank-1) << endl;    
//     cout << "arank=" << arank << "\tprev arank=" << ngood[h-1] << "\th=" << h << "\tnc=" << nc << endl;
    arma::Mat<T> W0(C.n_rows,arank);
    gMThreadException te;
#pragma omp parallel for
    for(arma::uword i=0;i<arank;i++){
        try{
            if(te.exceptionOccurred()) continue;
            arma::uvec sel=arma::find(cl==goods(i));
            W0.col(i)=arma::sqrt(arma::mean(arma::Mat<T>(square(C.cols(sel))), 1 ));
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
//     cout << "nc=" << nc << "\tarank=" << arank << endl;
    
    gMDense<T> silhconv=armadillo::convert<arma::Mat,gMDense,T>(silh);
    gMDense<indexType> ngoodconv=armadillo::convert<arma::Mat,gMDense,indexType>(ngood);
    
    return make_pair(move(armadillo::convert<arma::Mat,gMDense,T>(W0)),make_pair(silhconv,ngoodconv));
}
template pair<gMDense<float>,pair<gMDense<float>,gMDense<indexType>>>  geco::methods::kdmotifs::kdmGetW0(const gMSparse<float> & counts,cluster::gHClusterMethod clustering_method,bool goodclusters, indexType minRank, indexType maxRank,bool bhattacharyya);
template pair<gMDense<double>,pair<gMDense<double>,gMDense<indexType>>> geco::methods::kdmotifs::kdmGetW0(const gMSparse<double> & counts,cluster::gHClusterMethod clustering_method,bool goodclusters, indexType minRank, indexType maxRank,bool bhattacharyya);

/// BEGIN UTILITIES 

#include <numeric>

template<typename T>
T binomial_test(unsigned int x, unsigned int n, T p, tailType tail){
    boost::math::binomial_distribution<T> dist(n, p);    
    T pv=1;
    if(x>0){

        switch (tail) {
            case tailType::LOWER:
                pv = cdf(dist, x); // P(X ≤ x)
                break;

            case tailType::UPPER:
                pv = cdf(complement(dist,x - 1)); // P(X ≥ x)
                break;

            case tailType::TWO_TAILED: 
                if(p == 0){
                    pv=(x==0);
                }else if(p == 1){
                    pv=(x==n);
                }else{
                    T relErr = 1.0 + 1.0e-07;
                    T d = pdf(dist, x);
                    unsigned int m = n * p;
                    auto cond = [dist,d,relErr](unsigned int acc,unsigned int x){ 
                        return acc + (pdf(dist,x) <= d * relErr);
                    };
                    if(x == m){
                        pv=1.0;
                    }else if(x < m){
                        unsigned int start=ceil(m);
                        vector<unsigned int> vec(n-start+1);
                        iota(vec.begin(),vec.end(),start);
                        unsigned int y = std::accumulate(vec.begin(),vec.end(), 0, cond);
//                         cout << "l\t" << n << "\t" << y << "\t" << n-y << endl;
                        pv = cdf(dist,x) + ((y<n)?cdf(complement(dist,n - y)):0);
                    } else {
                        unsigned int end=floor(m);
                        vector<unsigned int> vec(end+1);
                        iota(vec.begin(),vec.end(),0);
                        
                        unsigned int y = std::accumulate(vec.begin(),vec.end(), 0, cond);
//                         cout << "g\t" << n << "\t" << y << "\t" << y-1 << endl;
                        pv = ((y>0)?cdf(dist,y-1):0) + cdf(complement(dist,x - 1));
                    }
                }    
                break;
                
            default:
                throw std::invalid_argument("Unknown tail type");
        }
    }
    return pv;
}

template<typename T>
std::tuple<T, T, T, T> fisherExactTest(int a, int b, int c, int d) {
    int N = a + b + c + d;
    int r = a+b;
    int n = a+c;
    T lowerP=1, upperP=1,twoTailedP=1,oddr=1; 
//     cout << "\tN=" << N << " n=" << n << " r=" << r << " n+r-N=" << n+r-N << endl;
    try{
        boost::math::hypergeometric_distribution<T> hg(r, n, N);

    //     T lowerP = boost::math::cdf(hg, a);
        if(a>0){
            if(a-1 > n+r-N){
    //             upperP = 1.0 - boost::math::cdf(hg, a - 1);
                upperP = boost::math::cdf(complement(hg, a - 1));
            }
        }
    //     T twoTailedP = min<T>(1.0, lowerP + upperP);
        
    //     cout << "\tCalcolo odd ratio" << endl;
        oddr = (b == 0 || c == 0) ? numeric_limits<T>::infinity() : (static_cast<T>(a) * d) / (b * c);
    }catch(exception &e){
        throw gMethodException(e.what());
    }catch(boost::exception & e){
        throw gMethodException("Something wrong occurred in binom_test");
    }
    return {oddr, lowerP, upperP, twoTailedP};
}

template<typename T>
struct PValue {
public:
    T p;
    indexType index;
};

template<typename T>
bool comparePValues(const PValue<T> & a, const PValue<T> & b) {
    return a.p < b.p;
}

template<typename T>
std::vector<T> benjaminiHochbergCorrection(const std::vector<T>& pvalues) {
    indexType n = pvalues.size();
    std::vector<PValue<T>> sortedPValues(n);
    std::vector<T> adjustedPValues(n);
//     cout << "\tAssocia i p-value agli indici originali" << endl;
    for (indexType i = 0; i < n; i++) {
        sortedPValues[i] = {pvalues[i], i};
    }
//     cout << "\tOrdina i p-value in ordine crescente" << endl;
    sort(sortedPValues.begin(), sortedPValues.end(), comparePValues<T>);
//     cout << "\tApplica la correzione di Benjamini-Hochberg" << endl;
    std::vector<T> temp(n);
    for (indexType i = 0; i < n; i++) {
        temp[i] = (sortedPValues[i].p * n) / (i + 1);
    }
//     cout << "\tGarantisce la monotonicità decrescente" << endl;
    for (int i = n - 2; i >= 0; i--) {
        temp[i] = min(temp[i], temp[i + 1]);
    }
//     cout << "\tRiordina i p-value corretti nelle posizioni originali" << endl;
    for (indexType i = 0; i < n; i++) {
        adjustedPValues[sortedPValues[i].index] = min<T>(temp[i], 1.0);
    }
    return move(adjustedPValues);
}

template <typename T>
tuple<T, T, vector<pair<T, T>>> performance(const vector<T>& predictions, const vector<int>& labels) {
    vector<pair<T, int>> scores;
    for (size_t i = 0; i < predictions.size(); ++i) {
        scores.emplace_back(predictions[i], labels[i]);
    }
    
    // Sort by predicted scores in descending order
    sort(scores.begin(), scores.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });
    
    vector<pair<T, T>> accuracy_per_threshold;
    T tp = 0, fp = 0, prev_tp = 0, prev_fp = 0;
    T auc = 0, pr_auc = 0;
    int total_pos = count(labels.begin(), labels.end(), 1);
    int total_neg = labels.size() - total_pos;
    int correct = total_pos;
    
    for (const auto& [score, label] : scores) {
        accuracy_per_threshold.emplace_back(score, static_cast<T>(correct) / labels.size());
        
        if (label == 1) {
            tp += 1;
            --correct;
        } else {
            fp += 1;
            ++correct;
            auc += (fp - prev_fp) * ((tp + prev_tp) / 2.0);
            pr_auc += (tp / (tp + fp)) * (tp - prev_tp) / total_pos;
            prev_tp = tp;
            prev_fp = fp;
        }
    }
    
    auc /= (total_pos * total_neg);
    return {auc, pr_auc, accuracy_per_threshold};
}
/// END UTILITIES 

template<typename T>
pair<gMDense<T>,gMDense<indexType>> geco::methods::kdmotifs::kdmRankInfo(const gMDense<T> & Distance,cluster::gHClusterMethod clustering_method){
    gHCluster dend(Distance,clustering_method);
    arma::Mat<indexType> ngood(Distance.i_nrows,1,arma::fill::zeros);
    arma::Mat<T> silh(Distance.i_nrows,1,arma::fill::zeros);
//     int comp=0;
#pragma omp parallel for
    for(indexType k=1;k<Distance.i_nrows;k++){
        std::pair<gMDense<T>,std::pair<gMDense<indexType>,gMDense<T>>> res=geco::methods::cluster::getSilhouette(Distance,dend.cutree_k(k+1));
        arma::Col<T> sc=armadillo::convert<gMDense,arma::Mat,T>(res.second.second).col(0);
        arma::uvec pp=find(sc>0);
        ngood(k) = pp.n_elem;
        silh(k) = arma::mean(sc);
// #pragma omp critical
//         {
//             cout << ++comp << endl;
//         }
    }
    return make_pair(armadillo::convert<arma::Mat,gMDense,T>(silh),armadillo::convert<arma::Mat,gMDense,indexType>(ngood));
}
template pair<gMDense<float>,gMDense<indexType>> geco::methods::kdmotifs::kdmRankInfo(const gMDense<float> & Distance,cluster::gHClusterMethod clustering_method);
template pair<gMDense<double>,gMDense<indexType>> geco::methods::kdmotifs::kdmRankInfo(const gMDense<double> & Distance,cluster::gHClusterMethod clustering_method);


template<typename T> 
motifProfileInfo<T> kdmProfileInfo_internal(const pair<pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T>> & motifProfile, T tolerance){
    motifProfileInfo<T> ret;
    indexType nmot = motifProfile.first.first.size();
    indexType npositions = motifProfile.first.first[0].i_ncols;
    indexType npseqs = motifProfile.first.first[0].i_nrows;
    indexType nnseqs = motifProfile.first.second[0].i_nrows;
    
    ret.pfeat=motifProfile.first.first;
    ret.nfeat=motifProfile.first.second;
    ret.x=motifProfile.second;
    
    vector<int> labels1(npseqs,1);
    vector<int> ln1(nnseqs,0);
    labels1.insert( labels1.end(), ln1.begin(), ln1.end() );
    
    vector<int> labels2(npseqs,0);
    vector<int> ln2(nnseqs,1);
    labels2.insert( labels2.end(), ln2.begin(), ln2.end() );
    
    ret.thresholds=vector<T>(nmot);
    ret.aurocs=vector<T>(nmot);
    
    arma::Mat<int> phits(nmot,npseqs,arma::fill::zeros);
    arma::Mat<int> nhits(nmot,nnseqs,arma::fill::zeros);
    arma::Mat<int> pcounts(nmot,npositions-2,arma::fill::zeros);
    arma::Mat<int> ncounts(nmot,npositions-2,arma::fill::zeros);
    
#pragma omp parallel for
    for(indexType m=0;m<nmot;++m){
        arma::Mat<T> apfeat=armadillo::convert<gMDense,arma::Mat,T>(ret.pfeat[m]);
        arma::Mat<T> anfeat=armadillo::convert<gMDense,arma::Mat,T>(ret.nfeat[m]);
        
        apfeat=arma::round(apfeat / tolerance) * tolerance;
        anfeat=arma::round(anfeat / tolerance) * tolerance;
       
        
        arma::Col<T> ap_max = arma::max(apfeat, 1);
        arma::Col<T> an_max = arma::max(anfeat, 1);
        
        arma::uvec ap_ind = arma::index_max(apfeat, 1);
        arma::uvec an_ind = arma::index_max(anfeat, 1);
        
        
        //a and b contain the maximum feature for each sequence
        vector<T> a=arma::conv_to< vector<T> >::from(ap_max);
        vector<T> b=arma::conv_to< vector<T> >::from(an_max);
        a.insert( a.end(), b.begin(), b.end() );
        auto [auroc1, aupr1, accuracy1] = performance<T>(a, labels1);
        auto [auroc2, aupr2, accuracy2] = performance<T>(a, labels2);
        //auto res2 = performance<T>(a, labels2);
        if(auroc1>auroc2){
            ret.thresholds[m]=(*max_element(accuracy1.begin(), accuracy1.end(), [](const auto& a, const auto& b){return a.second > b.second;})).first;
            ret.aurocs[m] = auroc1;
        }else{
            ret.thresholds[m]=(*max_element(accuracy2.begin(), accuracy2.end(), [](const auto& a, const auto& b){return a.second > b.second;})).first;
            ret.aurocs[m] = - auroc2;
        }
        
        arma::uvec am={m};
        arma::uvec good_ap=arma::find(ap_max>ret.thresholds[m]);
        arma::uvec good_an=arma::find(an_max>ret.thresholds[m]);
        arma::uvec good_ap_ind=ap_ind.elem(good_ap).as_row();
        arma::uvec good_an_ind=an_ind.elem(good_an).as_row();
        
        for(arma::uword i=0;i<good_ap_ind.n_elem;++i){
            T val=ap_max(good_ap(i));
            
            arma::uvec sele=find((val-apfeat.row(good_ap[i]))<tolerance);
            good_ap_ind[i]=sele(0);
            
            arma::uword mpos=good_ap_ind[i];
            arma::uword start=mpos;
            arma::uword end=mpos;
            bool goodstart=true,goodend=true;
            while(goodstart || goodend){
                if(goodstart){
                    if(start>0){
                        if((val-apfeat(good_ap(i),start-1))<=tolerance){
                            start--;
                        }else{
                            goodstart=false;
                        }
                    }else{
                        goodstart=false;
                    }
                }
                if(goodend){
                    if(end+1<npositions){
                        if((val-apfeat(good_ap(i),end+1))<=tolerance){
                            end++;
                        }else{
                            goodend=false;
                        }
                    }else{
                        goodend=false;
                    }
                    
                }
            }
//             cout << "m=" << m << " POS seq=" << seq+1 << " npos=" << npositions <<  " val=" << val << " mpos=" << mpos+1 <<" [" << start+1 << "," << end+1 << "]" << endl;
            good_ap_ind(i)=arma::mean(arma::linspace<arma::uvec>(start,end,end-start+1));
        }
        
        
        for(arma::uword i=0;i<good_an_ind.n_elem;++i){
            T val=an_max(good_an(i));
            arma::uvec sele=find((val-anfeat.row(good_an[i]))<tolerance);
            good_an_ind[i]=sele(0);

            arma::uword mpos=good_an_ind[i];
            arma::uword start=mpos;
            arma::uword end=mpos;
            bool goodstart=true,goodend=true;
            while(goodstart || goodend){
                if(goodstart){
                    if(start>0){
                        if((val-anfeat(good_an(i),start-1))<=tolerance){
                            start--;
                        }else{
                            goodstart=false;
                        }
                    }else{
                        goodstart=false;
                    }
                }
                if(goodend){
                    if(end+1<npositions){
                        if((val-anfeat(good_an(i),end+1))<=tolerance){
                            end++;
                        }else{
                            goodend=false;
                        }
                    }else{
                        goodend=false;
                    }
                }
            }
//             cout << "m=" << m << " NEG seq=" << seq+1 << " npos=" << npositions <<  " val=" << val << " mpos=" << mpos+1 <<" [" << start+1 << "," << end+1 << "]" << endl;
            good_an_ind(i)=arma::mean(arma::linspace<arma::uvec>(start,end,end-start+1));
        }
        
        phits.submat(am,good_ap)=arma::conv_to<arma::Row<int>>::from(good_ap_ind) + 1;
        nhits.submat(am,good_an)=arma::conv_to<arma::Row<int>>::from(good_an_ind) + 1;
        
        for(indexType s=0;s<npseqs;++s){
            indexType phit=phits(m,s);
            if(phit>1 & phit<npositions){
                pcounts(m,phit-2)++;
            }
        }
        
        for(indexType s=0;s<nnseqs;++s){
            indexType nhit=nhits(m,s);
            if(nhit>1 & nhit<npositions){
                ncounts(m,nhit-2)++;
            }
        }
    }
    
    ret.phits=armadillo::convert<arma::Mat,gMDense,int>(phits);
    ret.nhits=armadillo::convert<arma::Mat,gMDense,int>(nhits);
    
    ret.pcounts=armadillo::convert<arma::Mat,gMDense,int>(pcounts);
    ret.ncounts=armadillo::convert<arma::Mat,gMDense,int>(ncounts);
    
    return move(ret);
}

template<typename T>
motifProfileInfo<T> geco::methods::kdmotifs::kdmProfileInfo(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,T tolerance,bool strict){
    return kdmProfileInfo_internal<T>(getMotifsProfile4(motifs,regions,centers,labels,genomeFile,halfInterval,halfWin,strict),tolerance);
}
template motifProfileInfo<float> geco::methods::kdmotifs::kdmProfileInfo(const gKDMotifs<float> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,float tolerance, bool strict);
template motifProfileInfo<double> geco::methods::kdmotifs::kdmProfileInfo(const gKDMotifs<double> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, double tolerance, bool strict);

template<typename T>
motifProfileInfo<T> geco::methods::kdmotifs::kdmProfileInfo(const gKDMotifs<T> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,T tolerance,bool strict){
    //auto motifProfile=getMotifsProfile4(motifs,pseq,nseq,halfInterval,halfWin,strict);
    return kdmProfileInfo_internal<T>(getMotifsProfile4(motifs,pseq,nseq,halfInterval,halfWin,strict),tolerance);
}
template motifProfileInfo<float> geco::methods::kdmotifs::kdmProfileInfo(const gKDMotifs<float> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,float tolerance, bool strict);
template motifProfileInfo<double> geco::methods::kdmotifs::kdmProfileInfo(const gKDMotifs<double> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,double tolerance, bool strict);


/* QUESTA ERA UNA PROVA PER L'APPROSSIMAZIONE
template<typename T> 
motifProfileInfo<T> kdmProfileInfo_internal2(const pair<pair<vector<gMDense<T>>,vector<gMDense<T>>>,vector<T>> & motifProfile, T tolerance){
    motifProfileInfo<T> ret;
    indexType nmot = motifProfile.first.first.size();
    indexType npositions = motifProfile.first.first[0].i_ncols;
    indexType npseqs = motifProfile.first.first[0].i_nrows;
    indexType nnseqs = motifProfile.first.second[0].i_nrows;
    
    ret.pfeat=motifProfile.first.first;
    ret.nfeat=motifProfile.first.second;
    ret.x=motifProfile.second;
    
    vector<int> labels1(npseqs,1);
    vector<int> ln1(nnseqs,0);
    labels1.insert( labels1.end(), ln1.begin(), ln1.end() );
    
    vector<int> labels2(npseqs,0);
    vector<int> ln2(nnseqs,1);
    labels2.insert( labels2.end(), ln2.begin(), ln2.end() );
    
    ret.thresholds=vector<T>(nmot);
    ret.aurocs=vector<T>(nmot);
    
    arma::Mat<int> phits(nmot,npseqs,arma::fill::zeros);
    arma::Mat<int> nhits(nmot,nnseqs,arma::fill::zeros);
    arma::Mat<int> pcounts(nmot,npositions-2,arma::fill::zeros);
    arma::Mat<int> ncounts(nmot,npositions-2,arma::fill::zeros);
    
#pragma omp parallel for
    for(indexType m=0;m<nmot;++m){
        arma::Mat<T> apfeat=armadillo::convert<gMDense,arma::Mat,T>(ret.pfeat[m]);
        arma::Mat<T> anfeat=armadillo::convert<gMDense,arma::Mat,T>(ret.nfeat[m]);

        
        apfeat=arma::round(apfeat / tolerance) * tolerance;
        anfeat=arma::round(anfeat / tolerance) * tolerance;

        
        
        arma::Col<T> ap_max = arma::max(apfeat, 1);
        arma::uvec ap_ind = arma::index_max(apfeat, 1);        

        
        for(arma::uword i=0;i<ap_max.n_elem;++i){
            T val=ap_max(i);
            
            arma::uvec sele=find((val-apfeat.row(i))<tolerance);
            ap_ind[i]=sele(0);
            
            arma::uword mpos=ap_ind[i];
            arma::uword start=mpos;
            arma::uword end=mpos;
            bool goodstart=true,goodend=true;
            while(goodstart || goodend){
                if(goodstart){
                    if(start>0){
                        if((val-apfeat(i,start-1))<=tolerance){
                            start--;
                        }else{
                            goodstart=false;
                        }
                    }else{
                        goodstart=false;
                    }
                }
                if(goodend){
                    if(end+1<npositions){
                        if((val-apfeat(i,end+1))<=tolerance){
                            end++;
                        }else{
                            goodend=false;
                        }
                    }else{
                        goodend=false;
                    }
                    
                }
            }
//             cout << "m=" << m << " POS seq=" << seq+1 << " npos=" << npositions <<  " val=" << val << " mpos=" << mpos+1 <<" [" << start+1 << "," << end+1 << "]" << endl;
            ap_ind(i)=arma::mean(arma::linspace<arma::uvec>(start,end,end-start+1));
            ap_max(i)=arma::mean(apfeat(i,arma::span(start,end)));
        }

        
        arma::Col<T> an_max = arma::max(anfeat, 1);
        arma::uvec an_ind = arma::index_max(anfeat, 1);
        
        for(arma::uword i=0;i<an_ind.n_elem;++i){
            T val=an_max(i);
            arma::uvec sele=find((val-anfeat.row(i))<tolerance);
            an_ind[i]=sele(0);

            arma::uword mpos=an_ind[i];
            arma::uword start=mpos;
            arma::uword end=mpos;
            bool goodstart=true,goodend=true;
            while(goodstart || goodend){
                if(goodstart){
                    if(start>0){
                        if((val-anfeat(i,start-1))<=tolerance){
                            start--;
                        }else{
                            goodstart=false;
                        }
                    }else{
                        goodstart=false;
                    }
                }
                if(goodend){
                    if(end+1<npositions){
                        if((val-anfeat(i,end+1))<=tolerance){
                            end++;
                        }else{
                            goodend=false;
                        }
                    }else{
                        goodend=false;
                    }
                }
            }
//             cout << "m=" << m << " NEG seq=" << seq+1 << " npos=" << npositions <<  " val=" << val << " mpos=" << mpos+1 <<" [" << start+1 << "," << end+1 << "]" << endl;
            an_ind(i)=arma::mean(arma::linspace<arma::uvec>(start,end,end-start+1));
            an_max(i)=arma::mean(anfeat(i,arma::span(start,end)));
        }
        //an_max=apfeat[an_id]
        
        
        
        //a and b contain the maximum feature for each sequence
        vector<T> a=arma::conv_to< vector<T> >::from(ap_max);
        vector<T> b=arma::conv_to< vector<T> >::from(an_max);
        a.insert( a.end(), b.begin(), b.end() );
        auto [auroc1, aupr1, accuracy1] = performance<T>(a, labels1);
        auto [auroc2, aupr2, accuracy2] = performance<T>(a, labels2);
        auto res2 = performance<T>(a, labels2);
        if(auroc1>auroc2){
            ret.thresholds[m]=(*max_element(accuracy1.begin(), accuracy1.end(), [](const auto& a, const auto& b){return a.second > b.second;})).first;
            ret.aurocs[m] = auroc1;
        }else{
            ret.thresholds[m]=(*max_element(accuracy2.begin(), accuracy2.end(), [](const auto& a, const auto& b){return a.second > b.second;})).first;
            ret.aurocs[m] = - auroc2;
        }
        
        arma::uvec am={m};
        arma::uvec good_ap=arma::find(ap_max-ret.thresholds[m]>tolerance);
        arma::uvec good_an=arma::find(an_max-ret.thresholds[m]>tolerance);
        arma::uvec good_ap_ind=ap_ind.elem(good_ap).as_row();
        arma::uvec good_an_ind=an_ind.elem(good_an).as_row();
        
        
        phits.submat(am,good_ap)=arma::conv_to<arma::Row<int>>::from(good_ap_ind) + 1;
        nhits.submat(am,good_an)=arma::conv_to<arma::Row<int>>::from(good_an_ind) + 1;
        
        for(indexType s=0;s<npseqs;++s){
            indexType phit=phits(m,s);
            if(phit>1 & phit<npositions){
                pcounts(m,phit-2)++;
            }
        }
        
        for(indexType s=0;s<nnseqs;++s){
            indexType nhit=nhits(m,s);
            if(nhit>1 & nhit<npositions){
                ncounts(m,nhit-2)++;
            }
        }
    }
    
    ret.phits=armadillo::convert<arma::Mat,gMDense,int>(phits);
    ret.nhits=armadillo::convert<arma::Mat,gMDense,int>(nhits);
    
    ret.pcounts=armadillo::convert<arma::Mat,gMDense,int>(pcounts);
    ret.ncounts=armadillo::convert<arma::Mat,gMDense,int>(ncounts);
    
    return move(ret);
}

template<typename T>
motifProfileInfo<T> geco::methods::kdmotifs::kdmProfileInfo2(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,T tolerance,bool strict){
    return kdmProfileInfo_internal2<T>(getMotifsProfile4(motifs,regions,centers,labels,genomeFile,halfInterval,halfWin,strict),tolerance);
}
template motifProfileInfo<float> geco::methods::kdmotifs::kdmProfileInfo2(const gKDMotifs<float> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin,float tolerance, bool strict);
template motifProfileInfo<double> geco::methods::kdmotifs::kdmProfileInfo2(const gKDMotifs<double> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const vector<unsigned int> & centers,const vector<unsigned int> & labels,const string & genomeFile,unsigned int halfInterval,unsigned int halfWin, double tolerance, bool strict);

template<typename T>
motifProfileInfo<T> geco::methods::kdmotifs::kdmProfileInfo2(const gKDMotifs<T> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,T tolerance,bool strict){
    //auto motifProfile=getMotifsProfile4(motifs,pseq,nseq,halfInterval,halfWin,strict);
    return kdmProfileInfo_internal2<T>(getMotifsProfile4(motifs,pseq,nseq,halfInterval,halfWin,strict),tolerance);
}
template motifProfileInfo<float> geco::methods::kdmotifs::kdmProfileInfo2(const gKDMotifs<float> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,float tolerance, bool strict);
template motifProfileInfo<double> geco::methods::kdmotifs::kdmProfileInfo2(const gKDMotifs<double> & motifs,const vector<string> & pseq, const vector<string> & nseq, unsigned int halfInterval,unsigned int halfWin,double tolerance, bool strict);
*/


/*
template<typename T>
vector< vector<vector<T> >> geco::methods::kdmotifs::motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,bool symmetric){
    indexType nmot=pcounts.i_nrows;
    indexType npos;
    arma::Mat<int> ap,an;
    if(symmetric){
        indexType onpos=pcounts.i_ncols;
        if(!(onpos & 1)){
            throw gMethodException("Symmetric enrichment requires an odd number of positions");
        }else{
            arma::Mat<int> oap=armadillo::convert<gMDense,arma::Mat,int>(pcounts);
            arma::Mat<int> oan=armadillo::convert<gMDense,arma::Mat,int>(ncounts);
            
            npos=(onpos-1)/2+1;
            indexType right_start = (onpos-1)/2 + 1;
            indexType right_end = right_start+(npos-1)-1;
            indexType left_start = 0;
            indexType left_end = left_start+(npos-1) -1;
            
            ap=arma::Mat<int>(nmot,npos,arma::fill::zeros);
            ap.col(0)=oap.col(right_start-1);
            ap.cols(1,npos-1) = oap.cols(right_start,right_end) + arma::reverse(oap.cols(left_start,left_end),1);
       
            an=arma::Mat<int>(nmot,npos,arma::fill::zeros);
            an.col(0)=oan.col(right_start-1);
            an.cols(1,npos-1) = oan.cols(right_start,right_end) + arma::reverse(oan.cols(left_start,left_end),1);
        }
    }else{
        npos=pcounts.i_ncols;
        ap=armadillo::convert<gMDense,arma::Mat,int>(pcounts);
        an=armadillo::convert<gMDense,arma::Mat,int>(ncounts);
    }
    indexType totrows=(npos*npos-npos)/2 + npos;
    vector< vector<vector<T>> > res(nmot,vector<vector<T>>(14,vector<T>(totrows)));
 
    gMThreadException te;    
#pragma omp parallel for
    for(indexType m=0;m<nmot;m++){
        try{
            vector<vector<T>> & resmot=res[m];
            auto motap = ap.row(m);
            auto motan = an.row(m);
            int np = arma::sum(motap);
            int nn = arma::sum(motan);
            
            indexType cc=0;
            for(indexType win=1; win<=npos; win++){
                T probp = ((T)win)/((T)npos);
                boost::math::binomial_distribution<T> bd1(np, probp);

                for(indexType pos=0;pos<(npos-win+1);pos++){
                    int max_pos = pos + motap.cols(pos,pos+win-1).index_max();
                    auto a=arma::linspace(pos,pos+win-1,win);
                    auto b=motap.cols(pos,pos+win-1);
                    auto c = a % b;
                    
                    int np_in = arma::sum(b);
                    T mean_pos = ((T) arma::sum(c)/(T) np_in);
//                     cout << "------------------" << endl;
//                     cout << c.n_elem << endl;
//                     cout << "------------------" << endl;
                  
//                     cout << "motif=" << m << "\twin=" << win << "\tpos=" << pos << "\tsum(c)=" << arma::sum(c) << "\tsum(b)=" << np_in << "\tratio=" << arma::sum(c)/np_in << endl;
                    int nn_in = arma::sum(motan.cols(pos,pos+win-1));                    
                    int np_out = np - np_in;
                    int nn_out = nn - nn_in;
                    T prob_in = ((T) nn_in)/((T) nnseqs);
                    boost::math::binomial_distribution<T> bd2(npseqs, prob_in);
                    
                    T cen_pv=(np_in>0)?cdf(complement(bd1,np_in-1)):1;
                    T cen = (((T) np_in)/((T)np))/probp;

                    T enr_pv=(np_in>0)?cdf(complement(bd2,np_in-1)):1;
                    T enr = (((T) np_in)/((T)npseqs))/prob_in;
                    
                    //auto [oddratio, pLower, pUpper, pTwoTailed] = fisherExactTest<T>(np_in, np_out, nn_in, nn_out);
                
                    resmot[0][cc]=pos+1;
                    resmot[1][cc]=pos+win;
                    resmot[2][cc]=max_pos+1;
                    resmot[3][cc]=mean_pos+1;
                    
                    
                    resmot[4][cc]=np_in;
                    resmot[5][cc]=np_out;
                    resmot[6][cc]=nn_in;
                    resmot[7][cc]=nn_out;
                    
                    resmot[8][cc]=cen;
                    resmot[9][cc]=cen_pv;
                   
                    resmot[11][cc]=enr;
                    resmot[12][cc]=enr_pv;

                    cc++;
                }
            }
            resmot[10] = benjaminiHochbergCorrection(resmot[9]);
            resmot[13] = benjaminiHochbergCorrection(resmot[12]);
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return move(res);
}
template vector< vector<vector<float>> > geco::methods::kdmotifs::motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,bool symmetric);
template vector< vector<vector<double>> > geco::methods::kdmotifs::motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,bool symmetric);
*/

template<typename T>
vector< vector<vector<T> >> geco::methods::kdmotifs::motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,tailType tail,bool symmetric){
    indexType nmot=pcounts.i_nrows;
    indexType npos;
    arma::Mat<int> ap,an;
    if(symmetric){
        indexType onpos=pcounts.i_ncols;
        if(!(onpos & 1)){
            throw gMethodException("Symmetric enrichment requires an odd number of positions");
        }else{
            arma::Mat<int> oap=armadillo::convert<gMDense,arma::Mat,int>(pcounts);
            arma::Mat<int> oan=armadillo::convert<gMDense,arma::Mat,int>(ncounts);
            
            npos=(onpos-1)/2+1;
            indexType right_start = (onpos-1)/2 + 1;
            indexType right_end = right_start+(npos-1)-1;
            indexType left_start = 0;
            indexType left_end = left_start+(npos-1) -1;
            
            ap=arma::Mat<int>(nmot,npos,arma::fill::zeros);
            ap.col(0)=oap.col(right_start-1);
            ap.cols(1,npos-1) = oap.cols(right_start,right_end) + arma::reverse(oap.cols(left_start,left_end),1);
       
            an=arma::Mat<int>(nmot,npos,arma::fill::zeros);
            an.col(0)=oan.col(right_start-1);
            an.cols(1,npos-1) = oan.cols(right_start,right_end) + arma::reverse(oan.cols(left_start,left_end),1);
        }
    }else{
        npos=pcounts.i_ncols;
        ap=armadillo::convert<gMDense,arma::Mat,int>(pcounts);
        an=armadillo::convert<gMDense,arma::Mat,int>(ncounts);
    }
    indexType totrows=(npos*npos-npos)/2 + npos;
    vector< vector<vector<T>> > res(nmot,vector<vector<T>>(14,vector<T>(totrows)));
 
    gMThreadException te;    
#pragma omp parallel for
    for(indexType m=0;m<nmot;m++){
        try{
            vector<vector<T>> & resmot=res[m];
            auto motap = ap.row(m);
            auto motan = an.row(m);
            int np = arma::sum(motap);
            int nn = arma::sum(motan);
            
            indexType cc=0;
            for(indexType win=1; win<=npos; win++){
                T probp = ((T)win)/((T)npos);
//                 boost::math::binomial_distribution<T> bd1(np, probp);

                for(indexType pos=0;pos<(npos-win+1);pos++){
                    int max_pos = pos + motap.cols(pos,pos+win-1).index_max();
                    auto a=arma::linspace(pos,pos+win-1,win);
                    auto b=motap.cols(pos,pos+win-1);
                    auto c = a % b;
                    
                    int np_in = arma::sum(b);
                    T mean_pos = ((T) arma::sum(c)/(T) np_in);
//                     cout << "------------------" << endl;
//                     cout << c.n_elem << endl;
//                     cout << "------------------" << endl;
                  
//                     cout << "motif=" << m << "\twin=" << win << "\tpos=" << pos << "\tsum(c)=" << arma::sum(c) << "\tsum(b)=" << np_in << "\tratio=" << arma::sum(c)/np_in << endl;
                    int nn_in = arma::sum(motan.cols(pos,pos+win-1));                    
                    int np_out = np - np_in;
                    int nn_out = nn - nn_in;
                    T prob_in = ((T) nn_in)/((T) nnseqs);
                  

                    T cen_pv=binomial_test<T>(np_in,np,probp,tail);
                    T cen = (((T) np_in)/((T)np))/probp;


                    T enr_pv=binomial_test<T>(np_in,npseqs,prob_in,tail);
                    T enr = (((T) np_in)/((T)npseqs))/prob_in;
                    
                    //auto [oddratio, pLower, pUpper, pTwoTailed] = fisherExactTest<T>(np_in, np_out, nn_in, nn_out);
                
                    resmot[0][cc]=pos+1;
                    resmot[1][cc]=pos+win;
                    resmot[2][cc]=max_pos+1;
                    resmot[3][cc]=mean_pos+1;
                    
                    
                    resmot[4][cc]=np_in;
                    resmot[5][cc]=np_out;
                    resmot[6][cc]=nn_in;
                    resmot[7][cc]=nn_out;
                    
                    resmot[8][cc]=cen;
                    resmot[9][cc]=cen_pv;
                   
                    resmot[11][cc]=enr;
                    resmot[12][cc]=enr_pv;

                    cc++;
                }
            }
            resmot[10] = benjaminiHochbergCorrection(resmot[9]);
            resmot[13] = benjaminiHochbergCorrection(resmot[12]);
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return move(res);
}
template vector< vector<vector<float>> > geco::methods::kdmotifs::motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,tailType tail,bool symmetric);
template vector< vector<vector<double>> > geco::methods::kdmotifs::motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,tailType tail,bool symmetric);


