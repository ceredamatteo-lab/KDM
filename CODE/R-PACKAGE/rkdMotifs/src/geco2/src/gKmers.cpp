#include "gKmers.h"
#include <fstream>
#include <set>
#include <list>
#include <algorithm>
#include <queue>

#ifdef GECO_HAS_CUDA
#include "gKmers_cuda.h"
#include "gMethods_cuda.h"
#endif
#include "gMethods_armadillo.h"

#include <nlopt.hpp>

using namespace std;
using namespace geco;
using namespace geco::methods;
using namespace geco::methods::kmers;

// // autodiff include
// #include <gMethods/autodiff/reverse/var.hpp>
// #include <gMethods/autodiff/reverse/var/eigen.hpp>
// using namespace autodiff;
// 
// var f(const ArrayXvar& x,const ArrayXMatrixXvar & C,const ArrayXvar & W){
//     
//     return W * Eigen::sqrt(x[C]); // sqrt(sum([xi * xi for i = 1:5]))
// }



static unsigned char cmap[256]={4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,0,4,1,4,4,4,2,0,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

int64_t bin(int64_t x){
    if(x==0) return 1;
    int64_t r = x;
    for(int64_t i = r-1;i>0;i--){
		r = r*i;
    }
    return r;
}

int64_t choose(int64_t n, int64_t k){
  int64_t result,fn,fk,fnk,p;
  if(n>0){
    fn = bin(n);
    fk = bin(k);
    fnk = bin(n-k);
    result = (double)fn/((double)(fk*fnk));
  }else{
    fn=bin(k-n-1);
    fk=bin(k);
    fnk=bin(-n-1);
    p=pow(-1,k);
    result= ((double)p)*((double)fn)/ ( ((double)fk) * ((double)fnk));
  }
  return result;
}

template<typename T>
struct isSmall{
    T i_nelm;
    isSmall(T nelm):i_nelm(nelm){
    }
    bool operator() (const T & val){
        return (val * i_nelm) < (T) 1;
        //return val < std::numeric_limits<T>::epsilon();
    }
};


vector<unsigned char > seqCode(const string & sequence){
    gMThreadException te;
    vector<unsigned char> tres(sequence.length());
    size_t sl=sequence.length();
    
#pragma omp parallel for
    for(indexType i=0;i<sl;i++){
        tres[i] = cmap[(size_t)sequence[i]];
    }
    
    return tres;
}

map<size_t,list<unsigned char>>  seqCode(const std::string  & seq,indexType sID,bool strict){
    
    vector<unsigned char> tres=seqCode(seq);
    // size_t sl=seq.length();    
    size_t sl=tres.size();    
    map<size_t,list<unsigned char>> res;
    size_t pos=0;

    while(pos<sl){
        while(pos<sl && tres[pos]==4){
            pos++;
        }
        auto plist=res.insert(make_pair(pos,list<unsigned char>()));
        while(pos < sl){
            if(tres[pos]!=4){
                (*plist.first).second.push_back(tres[pos]);
                pos++;
            }
        }
    }
    if(res.size()>1){
        string msg=string("gLmerCounter: sequence ")+seq+" ("+to_string(sID+1)+") " + string(" contains invalid characters");
        if(strict){
            throw gMethodException(msg);
        }else{
            warning(msg);
        }
    }
    return res;
}

double ICpwm(const arma::Mat<double> & pwm){
 arma::Col<double> ics(pwm.n_cols);
 for(indexType i=0;i<pwm.n_cols;i++){
   arma::Col<double> w=pwm.col(i);
   arma::Col<indexType> sel=find(w>0);
   ics(i) = 2+arma::sum( w(sel) % arma::log2(w(sel)) );
 }
 ics/=(double)2;
 for(indexType i=0;i<8;i++){
     ics(i)=-ics(i)/8;
 }
 for(indexType i=8;i<33;i++){
     ics(i)=ics(i)/25;
 }
 for(indexType i=33;i<41;i++){
     ics(i)=-ics(i)/8;
 }
 
 return arma::sum(ics);
} 



//BEGIN ROBA PER CALCOLARE PWM da W senza GAP

double pwmfunc(unsigned n, const double *x, double *grad, void * fundata);
void pwmConstraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);

class pwmdata{
public:
    indexType i_nlmers;
    indexType i_nDSlmers;
    indexType i_l;
    indexType i_ncols;
    indexType i_ncols_padded;
    indexType i_nsteps;
    arma::Col<double> i_w;
    indexType ncall=0;
    
    const vector<indexType> & i_lmerDSLmer;
    
    arma::Cube<indexType> IND;
    arma::Col<indexType> R;
    double * paddedpwm;
    double * pwm;
    nlopt::opt  & i_opt;
    
    
    template<typename T> 
    pwmdata():i_nlmers(0),i_l(0),i_ncols(0),i_ncols_padded(0),i_nsteps(0){
        paddedpwm=NULL;
        pwm=NULL;
    }
    
    pwmdata(indexType l,indexType nlmers,indexType nDSlmers,indexType ncols,const arma::Col<double> & w,const vector<double> & pwm0,double tol, indexType maxIter,const vector<indexType> & lmerDSLmer,nlopt::opt & myopt):i_nlmers(nlmers),i_nDSlmers(nDSlmers),i_l(l),i_ncols(ncols),i_ncols_padded(i_ncols + 2*(i_l-1)),i_nsteps(i_ncols_padded-i_l+1),i_w(w),i_lmerDSLmer(lmerDSLmer),IND(i_l,i_nsteps,i_nlmers),i_opt(myopt){
        R=arma::Col<indexType>({0,1,2});
        paddedpwm=new double[4*i_ncols_padded];
        pwm = paddedpwm+4*(i_l-1);
        for(indexType i=0;i<4*i_ncols_padded;i++){
            paddedpwm[i]=0.25;
        }
        copy(pwm0.begin(),pwm0.end(),pwm);
        arma::Col<indexType> DISP(i_l);
        for(indexType i=0;i<i_l;i++){
            DISP[i] = i*4;
        }
        vector<string> alpha({"A","C","G","T"});
#pragma omp parallel for
        for(indexType lmer=0;lmer<i_nlmers;lmer++){
            arma::Col<indexType> C(i_l);
            indexType b=lmer;
            for(indexType pos=0;pos<i_l;pos++){
                C(i_l-pos-1)=b & 0x00000003;
                b >>= 2;
            }
            for(indexType step=0;step<i_nsteps;step++){
                IND.slice(lmer).col(step)=C+DISP+step*4;
            }
        }
    }
    
    ~pwmdata(){
        if(paddedpwm!=NULL){
            delete [] paddedpwm;
        }
    }
};

double pwmfunc(unsigned n, const double *x, double *grad, void * fundata){
    pwmdata & data=*(pwmdata *) fundata;
    arma::Mat<double> M0(x,3,data.i_ncols);
    arma::Mat<double> M(4,data.i_ncols);        
    M.rows(data.R)=M0;
    M.row(3)=arma::clamp(-arma::sum(M0)+1,0.0,1.0);
    copy(M.memptr(),M.memptr()+4*data.i_ncols,data.pwm);
    arma::Col<double> mat(data.paddedpwm,data.i_ncols_padded*4,false,false);
    
    arma::Col<double> tmp(data.i_nlmers);
#pragma omp parallel for
    for(indexType lmer=0;lmer<data.i_nlmers;lmer++){
        double ptot=0;
        arma::Mat<indexType> MI=data.IND.slice(lmer);
        for(indexType step=0;step<data.i_nsteps;step++){
            ptot+=arma::prod(mat(MI.col(step)));
        }
        tmp(lmer)=ptot/(double)data.i_nsteps;
    }
    arma::Col<double> westimate((data.i_lmerDSLmer.size()>0)?(data.i_nDSlmers):(data.i_nlmers),arma::fill::zeros);
    if(data.i_lmerDSLmer.size()>0){
        for(indexType lmer=0;lmer<data.i_nlmers;lmer++){
            westimate(data.i_lmerDSLmer[lmer])+=tmp[lmer];
        }
    }else{
        westimate=tmp;
    }
    westimate/=arma::sum(westimate);
    double res=arma::dot(arma::sqrt(westimate),data.i_w);
    
    data.ncall++;
    if(data.ncall % 10==0){
        cout << "iteration: " << data.ncall << "\t" << res << endl;
//             inner.End([=]{ std::cout << "\tchange: " << change;});
//             inner.Start([=]{std::cout << "\titerations [" << (data.ncall+1) << "-" << (data.ncall+10) << "]...";});
        try{
            checkInterrupt();
        }catch(gMethodsInterrupt &e){
            warning("User intrruption before convergence");
            data.i_opt.force_stop();
        }
    }
    return res;
}

void pwmConstraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data){
    arma::Mat<double> M0(x,3,m);
    arma::Col<double> R=sum(M0)-1;
//     R.print();
    copy(R.begin(),R.end(),result);
}
//END ROBA PER CALCOLARE PWM da W senza GAP



template<typename T>
gLmerCounter<T>::gLmerCounter(indexType l,bool doubleStrand):i_doubleStrand(doubleStrand),i_l(l),i_nlmers(pow(4,l)){
//     auto prime_time = std::chrono::high_resolution_clock::now();    
//     std::chrono::duration<double> dur;
    
//     cout << "N. lmers: \t" << i_nlmers << endl;
    if(i_doubleStrand){
        //MAPPA lmers -> dslmers
//         cout << "building lmers -> dslmers...";cout.flush();
//         prime_time = std::chrono::high_resolution_clock::now();
        {
            arma::Col<indexType> list1(i_nlmers,arma::fill::zeros);
            arma::Col<indexType> list2(i_nlmers,arma::fill::zeros);
            indexType v,c;
            for(indexType i=0;i<i_nlmers;i++){
                v=i;
                c=0;
                for(indexType h=0;h<(i_l-1);h++){
                    c|= ( (~v) & 3);
                    v>>=2;
                    c<<=2;
                }
                c|= ( (~v) & 3);
                indexType b=min(i,c);
                list1[i]=b;
                list2[b]=1;
            }
            arma::uvec pos=arma::find(list2);
            arma::Col<indexType> i_values=list1(pos);
            list2.elem(pos)=arma::regspace<arma::Col<indexType>>(0,pos.n_elem);
            arma::Col<indexType> i_map=list2(list1);
            i_lmerDSLmer = vector<indexType>(i_map.begin(),i_map.end());
            i_DSlmers = vector<indexType>(i_values.begin(),i_values.end());
            i_nDSlmers=i_values.size();
        }
//         dur=std::chrono::high_resolution_clock::now()-prime_time;
//         cout << "\t(" << dur.count() << ")" << endl;
//         cout << "N. double strand lmers: \t" << i_nDSlmers << endl;
    }else{
        i_nDSlmers=0;
    }
}

template<typename T>
gMDense<T> gLmerCounter<T>::countDense(const vector<string> & sequences,bool strict) const{
    auto prime_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;    
    
    gMThreadException te;    
    indexType i_mask=pow(4,i_l)-1;
    indexType nrows=i_doubleStrand?i_nDSlmers:i_nlmers;    
    
    prime_time = std::chrono::high_resolution_clock::now();
    vector<T> lmers(sequences.size() * nrows);
    dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "time1: " << dur.count()*1000 << "ms" << endl;
    prime_time = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
    for(indexType s=0;s<sequences.size();s++){
        if(te.exceptionOccurred()) continue;
        try{
            indexType base=s*nrows;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                auto bb=u->second.begin();
                if(u->second.size()>=i_l){
                    nnzero=true;
                    indexType lmer=*bb++;
                    for(indexType i=1;i<i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    if(i_doubleStrand){
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            lmers[base + i_lmerDSLmer[lmer & i_mask]]++;
//                             cout << ".";
                        }
                    }else{
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            lmers[base + (lmer & i_mask)]++;
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gLmerCounter: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
//     cout << endl;
    dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "time2: " << dur.count()*1000 << "ms" << endl;
    te.Rethrow();
    return move(gMDense<T>(nrows,sequences.size(),lmers.data()));
}

template<typename T>
gMSparse<T> gLmerCounter<T>::countSparse(const vector<string> & sequences,bool strict) const{
    gMThreadException te;    
    indexType i_mask=pow(4,i_l)-1;
    indexType nrows=i_doubleStrand?i_nDSlmers:i_nlmers;    

    auto prime_time = std::chrono::high_resolution_clock::now();
    vector<map<indexType,T>> lmers(sequences.size());    
#pragma omp parallel for
    for(indexType s=0;s<sequences.size();s++){
        if(te.exceptionOccurred()) continue;
        try{
            map<indexType,T> tmers;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                auto bb=u->second.begin();
                if(u->second.size()>=i_l){
                    nnzero=true;
                    indexType lmer=*bb;
                    bb++;
                    for(indexType i=1;i<i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    if(i_doubleStrand){
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            //lmers[s][i_lmerDSLmer[lmer & i_mask]]++;
                            tmers[i_lmerDSLmer[lmer & i_mask]]++;
                        }
                    }else{
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            //lmers[s][lmer & i_mask]++;
                            tmers[lmer & i_mask]++;
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gLmerCounter: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            lmers[s]=tmers;
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "time: " << dur.count()*1000 << "ms" << endl;
    te.Rethrow();
    return move(gMSparse<T>(lmers,nrows));
}

template<typename T>
gMDense<T> gLmerCounter<T>::countWinDense(const string & sequence,indexType win,bool strict) const{
    auto prime_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;    
    
    gMDense<T> RES;
    indexType i_mask=pow(4,i_l)-1;
    indexType nrows=i_doubleStrand?i_nDSlmers:i_nlmers;    
    indexType ncols=sequence.size()-win+1;
    
    prime_time = std::chrono::high_resolution_clock::now();
    vector<T> lmers(ncols * nrows);
    dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "time1: " << dur.count()*1000 << "ms" << endl;
    prime_time = std::chrono::high_resolution_clock::now();    
    
    vector<unsigned char> tres=seqCode(sequence);
    
    if(find(tres.begin(),tres.end(),4)==tres.end()){
        queue<indexType> qw;
        indexType mlmer;
        indexType base=0;
        indexType base2=0;
        

        auto bb=tres.begin();
        indexType lmer=*bb++;
        for(indexType i=1;i<i_l-1;i++){
            lmer <<= 2;
            lmer|=*bb++;
        }

        if(i_doubleStrand){
            for(indexType pos=i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=i_lmerDSLmer[lmer & i_mask];
                qw.push(mlmer);
                lmers[base2+mlmer]++;
            }
            while(base<ncols-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=i_lmerDSLmer[lmer & i_mask];
                qw.push(mlmer);
                lmers[base2+mlmer]++;
//                 lmers[base2+i_lmerDSLmer[lmer & i_mask]]++;

//                 indexType a=base2;
//                 indexType b=(base+1)*nrows;
// #pragma omp parallel for
//                 for(indexType i=0;i<nrows;i++){
//                     lmers[b+i]=lmers[a+i];
//                 }
                copy_n(lmers.begin()+base2,nrows,lmers.begin()+(base+1)*nrows);
                base2 = (++base) * nrows;
                lmers[base2+qw.front()]--;
                qw.pop();
//                 cout << ".";
            }
            lmer <<= 2;
            lmer|=*bb++;
            lmers[base2+i_lmerDSLmer[lmer & i_mask]]++;
        }else{
            for(indexType pos=i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=lmer & i_mask;
                qw.push(mlmer);
                lmers[base2+mlmer]++;
            }
            while(base<ncols-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=lmer & i_mask;
                qw.push(mlmer);
                lmers[base2+mlmer]++;
                copy_n(lmers.begin()+base2,nrows,lmers.begin()+(base+1)*nrows);
                base2=(++base) * nrows;
                lmers[base2+qw.front()]--;
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            lmers[base2+(lmer & i_mask)]++;
        }
//         cout << endl;
        dur=std::chrono::high_resolution_clock::now()-prime_time;
//         cout << "time2: " << dur.count()*1000 << "ms" << endl;

        RES=gMDense<T>(nrows,ncols,lmers.data());
    }else{ //Invalid charaters we do it thes slow way
        cout << "The slow way" << endl;
        vector<string> seqs(ncols);
        for(indexType p=0;p<ncols;p++){
            seqs[p] = sequence.substr(p,win);
        }
        RES=countDense(seqs,strict);
    }
    return move(RES);
}


template<typename T>
gMSparse<T> gLmerCounter<T>::countWinSparse(const string & sequence,indexType win,bool strict) const{
    
    gMSparse<T> RES;    
    indexType i_mask=pow(4,i_l)-1;
    indexType nrows=i_doubleStrand?i_nDSlmers:i_nlmers;    
    indexType ncols=sequence.size()-win+1;
    
    auto prime_time = std::chrono::high_resolution_clock::now();
    vector<map<indexType,T>> lmers(ncols);    
    vector<unsigned char> tres=seqCode(sequence);
    
    if(find(tres.begin(),tres.end(),4)==tres.end()){
        map<indexType,T> ores;
        queue<indexType> qw;
        indexType mlmer;
        indexType base=0;
    
        auto bb=tres.begin();
        indexType lmer=*bb++;
        for(indexType i=1;i<i_l-1;i++){
            lmer <<= 2;
            lmer|=*bb++;
        }

        if(i_doubleStrand){
            for(indexType pos=i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer = i_lmerDSLmer[lmer & i_mask];
                qw.push(mlmer);
                ores[mlmer]++;
            }
            while(base<ncols-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer = i_lmerDSLmer[lmer & i_mask];
                qw.push(mlmer);
                ores[mlmer]++;
                lmers[base]=ores;
                base++;
                if( (--ores[qw.front()]) == 0){
                    ores.erase(qw.front());
                }
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            ores[i_lmerDSLmer[lmer & i_mask]]++;
            lmers[base]=ores;
        }else{
            for(indexType pos=i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                qw.push(lmer & i_mask);
                ores[qw.back()]++;
            }
            while(base<ncols-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=lmer & i_mask;
                qw.push(mlmer);
                ores[mlmer]++;
                lmers[base]=ores;
                base++;
                if((--ores[qw.front()]) == 0){
                    ores.erase(qw.front());
                }
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            ores[lmer & i_mask]++;
            lmers[base]=ores;
        }
        std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
//         cout << "time: " << dur.count()*1000 << "ms" << endl;
        
        RES=gMSparse<T>(lmers,nrows);
    }else{ //Invalid charaters we do it thes slow way
        cout << "The slow way" << endl;
        vector<string> seqs(ncols);
        for(indexType p=0;p<ncols;p++){
            seqs[p] = sequence.substr(p,win);
        }
        RES=countSparse(seqs,strict);
    }
    return move(RES);
}

template<typename T>
gMDense<T> gLmerCounter<T>::countPWM(const geco::kmers::gPWMSet & PWMs,indexType padSize) const{
    if(padSize>=i_l){
        throw gException("countPWM: padSize must be lower than kmer length");
    }
    //gSize padSize=i_l-1;
    indexType nrows=i_doubleStrand?i_nDSlmers:i_nlmers;    
    vector<T> lmers(PWMs.size() * nrows,0);
    

	gMThreadException te;
#pragma omp parallel for
	for(indexType a=0;a<PWMs.size();a++){
        indexType base= a * nrows;
        if(te.exceptionOccurred()) continue;
		try{
			auto m=PWMs.begin();
			advance(m,a);
			
			gArray<T> background=m->second.getBackground();
			gMatrix<T> pwm=m->second.getPWM();

			gArray<T> pads;
			for(indexType i=0;i<padSize;i++){
				pads.concatenate(background);
			}
		
			gArray<T> full=pads;
			for(indexType i=0;i<pwm.getColsNum();i++){
				full.concatenate(pwm.getCol(i));
			}
			full.concatenate(pads);
			gMatrix<T> padded(2*padSize+pwm.getColsNum(),4,full);

			
			const T * apads=padded.getConstData();
			indexType nstep=padded.getRowsNum()-i_l+1;
			indexType cols[i_l];

			for(indexType lmer=0;lmer<i_nlmers;lmer++){
				indexType b=lmer;
				for(indexType i=0;i<i_l;i++){
					cols[i_l-i-1]=b & 0x00000003;
					b >>= 2;
				}
				T ptot=0;
				for(indexType p=0;p<nstep;p++){
					T prob=1;
					for(indexType i=0;i<i_l;i++){
						prob*=apads[4*(p+i)+cols[i]];
					}
					ptot+=prob;
				}
				ptot/=nstep;
				if(i_doubleStrand){
                    lmers[base + i_lmerDSLmer[lmer]]+=ptot;
                }else{
                    lmers[base + lmer]+=ptot;
                }
			}
		}catch(...){
			te.CaptureException();
		}
	}
	te.Rethrow();
    //replace_if(lmers.begin(), lmers.end(),isSmall<T>(nrows),0.0);
 	return move(gMDense<T>(nrows,PWMs.size(),lmers.data()));
    
}

template<typename T>
gMDense<T> gLmerCounter<T>::getPWM(const gMDense<T> & W,indexType pwmlength,const vector<double> & pwm0,double tol,indexType maxiter) const{
    nlopt::opt myopt(nlopt::LN_COBYLA,pwmlength*3);
    arma::Col<double> w=arma::conv_to<arma::Mat<double>>::from(armadillo::convert<gMDense,arma::Mat,T>(W)).col(0);
    pwmdata data(i_l,i_nlmers,i_nDSlmers,pwmlength,w,pwm0,tol,maxiter,i_lmerDSLmer,myopt);
    
    myopt.set_max_objective(pwmfunc, &data);
    myopt.set_stopval(1.0-tol);
    myopt.set_ftol_rel(tol);
    myopt.set_maxeval(maxiter);
    myopt.set_lower_bounds(0.0);
    myopt.set_upper_bounds(1.0);
    myopt.add_inequality_mconstraint(pwmConstraint, NULL, vector<double>(pwmlength,1e-8));
    
   
    double minf;
    arma::Mat<double> M0(pwm0.data(),4,data.i_ncols);
    arma::Mat<double> M1=M0.rows(data.R);
    vector<double> v(M1.begin(),M1.end());
    myopt.optimize(v, minf);
    arma::Mat<double> M2(v.data(),3,data.i_ncols);
    arma::Mat<double> M(4,data.i_ncols);
    M.rows(data.R)=M2;
    M.row(3)=1-arma::sum(M2);
    arma::Col<indexType> bad=find(M.row(3)<0);
    if(bad.n_elem>0){
        for(indexType i=0;i<bad.n_elem;i++){
            arma::Col<double> c=M.col(bad(i));
            M.col(bad(i))=(c-c(3))/sum(c-c(3));
        }
        warning("Negative entries in pwm due to approximation in constraints: corrected");        
    }
    vector<T> ret(4*pwmlength);
    copy(M.memptr(),M.memptr()+4*pwmlength,ret.begin());
    return move(armadillo::convert<arma::Mat,gMDense,T>( arma::conv_to<arma::Mat<T>>::from(M)));
}

template<typename T>
gMDense<T> gLmerCounter<T>::countPWM2(const geco::kmers::gPWMSet & PWMs) const{
	gSize padSize=i_l-1;
    indexType nrows=i_doubleStrand?i_nDSlmers:i_nlmers;    
    vector<T> lmers(PWMs.size() * nrows,0);
    

	gMThreadException te;
#pragma omp parallel for
	for(indexType a=0;a<PWMs.size();a++){
        indexType base= a * nrows;
        if(te.exceptionOccurred()) continue;
		try{
			auto m=PWMs.begin();
			advance(m,a);
			
			gArray<T> background=m->second.getBackground();
			gMatrix<T> pwm=m->second.getPWM();

			gArray<T> pads;
			for(indexType i=0;i<padSize;i++){
				pads.concatenate(background);
			}
		
			gArray<T> full=pads;
			for(indexType i=0;i<pwm.getColsNum();i++){
				full.concatenate(pwm.getCol(i));
			}
			full.concatenate(pads);
			gMatrix<T> padded(2*padSize+pwm.getColsNum(),4,full);

			
			const T * apads=padded.getConstData();
			indexType nstep=padded.getRowsNum()-i_l+1;
			indexType cols[i_l];

			for(indexType lmer=0;lmer<i_nlmers;lmer++){
				indexType b=lmer;
				for(indexType i=0;i<i_l;i++){
					cols[i_l-i-1]=b & 0x00000003;
					b >>= 2;
				}
				
                gArray<T> probs(1,nstep,false);
				for(indexType p=0;p<nstep;p++){
                    T prob=1;
					for(indexType i=0;i<i_l;i++){
						prob*=apads[4*(p+i)+cols[i]];
					}
					probs.setValue(p,prob,false);
				}
				T ptot=probs.getMax()[0];
				if(i_doubleStrand){
                    lmers[base + i_lmerDSLmer[lmer]]+=ptot;
                }else{
                    lmers[base + lmer]+=ptot;
                }
			}
		}catch(...){
			te.CaptureException();
		}
	}
	te.Rethrow();
    replace_if(lmers.begin(), lmers.end(),isSmall<T>(nrows),0.0);
 	return move(gMDense<T>(nrows,PWMs.size(),lmers.data()));
    
}


template<typename T>
pair< vector<std::string>,vector<string> > gLmerCounter<T>::getStrings() const{
    pair<vector<string>,vector<string>> ret;
    if(i_doubleStrand){
        vector<string> falpha={"A","C","G","T"};
        vector<string> ralpha={"T","G","C","A"};
        ret.first=vector<string>(i_nDSlmers);
        ret.second=vector<string>(i_nDSlmers);
#pragma omp parallel for        
        for(indexType dslmer=0;dslmer<i_nDSlmers;dslmer++){
            indexType lmer=i_DSlmers[dslmer];
            for(indexType i=0;i<i_l;i++){
                ret.first[dslmer] = falpha[lmer & 3] + ret.first[dslmer];
                ret.second[dslmer] = ret.second[dslmer] + ralpha[lmer & 3];
                lmer >>= 2;
            }
//             cout << ret.first[dslmer] << "\t" << ret.second[dslmer] << endl;
        }
    }else{
        vector<string> falpha={"A","C","G","T"};
        ret.first=vector<string>(i_nlmers);
        ret.second=vector<string>(i_nlmers);
#pragma omp parallel for                
        for(indexType lmer=0; lmer< i_nlmers; lmer++){
            indexType lmcode=lmer;
            for(indexType i=0;i<i_l;i++){
                ret.first[lmer] = falpha[lmcode & 3] + ret.first[lmer];
                lmcode >>= 2;
            }
//             cout << ret.first[lmer] << endl;
        }
    }
    return ret;
}

template class geco::methods::kmers::gLmerCounter<float>;
template class geco::methods::kmers::gLmerCounter<double>;
template class geco::methods::kmers::gLmerCounter<indexType>;
template class geco::methods::kmers::gLmerCounter<unsigned short>;




//BEGIN ROBA PER CALCOLARE PWM da W CON GAP

double pwmfuncG(unsigned n, const double *x, double *grad, void * fundata);
void pwmConstraintG(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data);

class pwmdataG{
public:
    indexType i_l;    
    indexType i_nlmers;
    indexType i_nDSgkmers;
    indexType i_ngkmers;
    indexType i_nmasks;
    
    indexType i_ncols;
    indexType i_ncols_padded;
    indexType i_nsteps;
    arma::Col<double> i_w;
    indexType ncall=0;
    
    const vector<indexType> & i_gkmerDSGkmer;
    const vector<indexType> & i_mlmerGkmer;
    
    arma::Cube<indexType> IND;
    arma::Col<indexType> R;
    double * paddedpwm;
    double * pwm;
    nlopt::opt  & i_opt;
    
    bool i_doubleStrand;    
    vector<arma::Col<indexType>> IND2;
//     gMethodsOut i_outer=gMethodsOut(1);    
//     gMethodsOut i_inner=gMethodsOut(2);
    
    
    template<typename T> 
    pwmdataG():i_l(0),i_nlmers(0),i_ncols(0),i_ncols_padded(0),i_nsteps(0){
        paddedpwm=NULL;
        pwm=NULL;
    }
    
    pwmdataG(indexType l,indexType nlmers,indexType nDSgkmers,indexType ngkmers,indexType nmasks,indexType ncols,const arma::Col<double> & w,const vector<double> & pwm0,const vector<indexType> & gkmerDSGkmer,const vector<indexType> & mlmerGkmer,nlopt::opt & myopt):i_l(l),i_nlmers(nlmers),i_nDSgkmers(nDSgkmers),i_ngkmers(ngkmers),i_nmasks(nmasks),i_ncols(ncols),i_ncols_padded(i_ncols + 2*(i_l-1)),i_nsteps(i_ncols_padded-i_l+1),i_w(w),i_gkmerDSGkmer(gkmerDSGkmer),i_mlmerGkmer(mlmerGkmer),IND(i_l,i_nsteps,i_nlmers),i_opt(myopt){
        R=arma::Col<indexType>({0,1,2});
        paddedpwm=new double[4*i_ncols_padded];
        pwm = paddedpwm+4*(i_l-1);
        for(indexType i=0;i<4*i_ncols_padded;i++){
            paddedpwm[i]=0.25;
        }
        copy(pwm0.begin(),pwm0.end(),pwm);
        arma::Col<indexType> DISP(i_l);
        for(indexType i=0;i<i_l;i++){
            DISP[i] = i*4;
        }
        vector<string> alpha({"A","C","G","T"});
        gMThreadException te;
#pragma omp parallel for
        for(indexType lmer=0;lmer<i_nlmers;lmer++){
            try{
                arma::Col<indexType> C(i_l);
                indexType b=lmer;
                for(indexType pos=0;pos<i_l;pos++){
                    C(i_l-pos-1)=b & 0x00000003;
                    b >>= 2;
                }
                for(indexType step=0;step<i_nsteps;step++){
                    IND.slice(lmer).col(step)=C+DISP+step*4;
                }
            }catch(...){
                te.CaptureException();
            }
        }
        te.Rethrow();
        
        vector<list<indexType>> INDTMP;
        if(i_nDSgkmers>0){
            i_doubleStrand=true;
            INDTMP=vector<list<indexType>>(i_nDSgkmers);
            for(indexType lmer=0;lmer<i_nlmers;lmer++){
                indexType mlmer = lmer * i_nmasks;
                for(indexType mask=0;mask<i_nmasks;mask++){
                    INDTMP[i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]].push_back(lmer);
                }                   
            }
        }else{
            i_doubleStrand=false;
            INDTMP=vector<list<indexType>>(i_ngkmers);
            for(indexType lmer=0;lmer<i_nlmers;lmer++){
                indexType mlmer = lmer * i_nmasks;
                for(indexType mask=0;mask<i_nmasks;mask++){
                    INDTMP[i_mlmerGkmer[mlmer++]].push_back(lmer);
                }                   
            }
        }
        
        IND2=vector<arma::Col<indexType>>(INDTMP.size());
        for(indexType i=0;i<INDTMP.size();i++){
            IND2[i]=arma::Col<indexType>(INDTMP[i].size());
            copy(INDTMP[i].begin(),INDTMP[i].end(),IND2[i].begin());
        }
        
    }
    
    ~pwmdataG(){
        if(paddedpwm!=NULL){
            delete [] paddedpwm;
        }
    }
};

double pwmfuncG(unsigned n, const double *x, double *grad, void * fundata){
    pwmdataG & data=*(pwmdataG *) fundata;
    data.ncall++;        
    cout << "iteration: " << data.ncall << "...\t";cout.flush();
    arma::Mat<double> M0(x,3,data.i_ncols);
    arma::Mat<double> M(4,data.i_ncols);        
    M.rows(data.R)=M0;
    M.row(3)=arma::clamp(-arma::sum(M0)+1,0.0,1.0);
    copy(M.memptr(),M.memptr()+4*data.i_ncols,data.pwm);
    arma::Col<double> mat(data.paddedpwm,data.i_ncols_padded*4,false,false);
    arma::Col<double> tmp(data.i_nlmers);
    
    gMThreadException te;
#pragma omp parallel for
    for(indexType lmer=0;lmer<data.i_nlmers;lmer++){
        try{
            double ptot=0;
            arma::Mat<indexType> MI=data.IND.slice(lmer);
            for(indexType step=0;step<data.i_nsteps;step++){
                ptot+=arma::prod(mat(MI.col(step)));
            }
            tmp(lmer)=ptot/(double)data.i_nsteps;
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    
    indexType nelm=data.i_doubleStrand?data.i_nDSgkmers:data.i_ngkmers;
    arma::Col<double> westimate(nelm,arma::fill::zeros);
    
#pragma omp parallel for
    for(indexType i=0;i<nelm;i++){
        westimate[i]=arma::sum(tmp(data.IND2[i]));
    }
    westimate/=arma::sum(westimate);
    
//     bool i_doubleStrand = data.i_gkmerDSGkmer.size()>0;
//     //indexType nelm=i_doubleStrand?data.i_nDSgkmers:data.i_ngkmers;
//     arma::Col<double> westimate(nelm,arma::fill::zeros);
//     
//     if(i_doubleStrand){
//         for(indexType lmer=0;lmer<data.i_nlmers;lmer++){
//             indexType mlmer = lmer * data.i_nmasks;
//             for(indexType mask=0;mask<data.i_nmasks;mask++){
//                 westimate[data.i_gkmerDSGkmer[data.i_mlmerGkmer[mlmer++]]]+=tmp[lmer];
//             }
//         }
//     }else{
//         for(indexType lmer=0;lmer<data.i_nlmers;lmer++){
//             indexType mlmer = lmer * data.i_nmasks;
//             for(indexType mask=0;mask<data.i_nmasks;mask++){
//                 westimate[data.i_mlmerGkmer[mlmer++]]+=tmp[lmer];
//             }                   
//         }
//     }
//     westimate/=arma::sum(westimate);

//     double res2=arma::dot(arma::sqrt(westimate2),data.i_w);
    double res=arma::dot(arma::sqrt(westimate),data.i_w);
    //double ic=ICpwm(M);

    if(data.ncall % 10==0){
//         cout << "iteration: " << data.ncall << "\t" << res << endl;
//         data.i_inner.End([=]{ std::cout << "\t" << res << "\t" << res2 << "\t" << res-res2;});
//         data.i_inner.Start([=]{ std::cout << "\titeration: " << data.ncall+1;});
        try{
            checkInterrupt();
        }catch(gMethodsInterrupt &e){
            warning("User intrruption before convergence");
            data.i_opt.force_stop();
        }
    }
    cout << res << endl;    
    return res;
}

void pwmConstraintG(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data){
    arma::Mat<double> M0(x,3,m);
    arma::Col<double> R=sum(M0)-1;
//     R.print();
    copy(R.begin(),R.end(),result);
}
//END ROBA PER CALCOLARE PWM da W CON GAP



template<typename T>
gGappedKmerCounter<T>::gGappedKmerCounter(geco::methods::indexType l,geco::methods::indexType k,bool doubleStrand):gLmerCounter<T>(l,doubleStrand),i_k(k),i_nmasks(choose(l,k)),i_nkmers(pow(4,k)),i_ngkmers(i_nmasks * i_nkmers),i_nmlmers(i_nmasks * gLmerCounter<T>::i_nlmers),i_masks(i_nmasks),i_mlmerGkmer(i_nmlmers){
    auto prime_time = std::chrono::high_resolution_clock::now();    
    std::chrono::duration<double> dur;
    
    if(k>l | k<1){
        throw gMethodException("error: k must be between 1 and l");
    }
//     cout << "N. Masks: \t" << i_nmasks << endl;
//     cout << "N. kmers: \t" << i_nkmers << endl;
//     cout << "N. gkmers:\t" << i_ngkmers << endl;
//     cout << "N. mlmers: \t" << i_nmlmers << endl;
    
    //MASKS MAP
//     cout << "building masks map";cout.flush();
    prime_time = std::chrono::high_resolution_clock::now();
    vector<indexType> i_mask_nmask(gLmerCounter<T>::i_nlmers);
    {
        indexType lsn=1;
        vector<indexType> selector(l);
        std::fill(selector.begin(), selector.begin() + k, (indexType) 1);
        indexType u=0;
        do{
            i_masks[u]=0;
            for (indexType  i = 0; i < l; i++) {
                if (selector[i]){
                    i_masks[u]|=(lsn << i);
                }
            }
            i_mask_nmask[i_masks[u]]=u;
            u++;
        }while(prev_permutation(selector.begin(), selector.end()));
//     for(indexType i=0;i<i_nmasks;i++){
//         cout << "A\t" << i << "\t" << i_masks[i] << endl;
//     }
//     
//     for(indexType i=0;i<i_nlmers;i++){
//         cout << i << "\t" << i_mask_nmask[i] << endl;
//     }
//     cout << endl;
    }
//     dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "\t(" << dur.count() << ")" << endl;
    
    //MAPPA mlmers -> gkmers
//     cout << "building mlmers -> gkmers map...";cout.flush();
    prime_time = std::chrono::high_resolution_clock::now();
    {
        indexType lsn=3;
#pragma omp parallel for
        for(indexType i=0;i<i_nmlmers;i++){
            indexType a = i/i_nmasks;
            indexType mask=i % i_nmasks;
            indexType b=i_masks[mask];
            indexType rk=0;
            for(indexType p=0;p<(l-1);p++){
                if(b & 1){
                    rk |= (a & 3);
                }
                b >>= 1;
                a >>= 2;
                if(b & 1){
                    rk <<= 2;
                }
            }
            if(b & 1){
                rk |= (a & 3);
            }
            
            indexType fk=0;
            for(indexType p=0;p<k-1;p++){
                fk|=(rk & lsn);
                fk <<= 2;
                rk >>= 2;
            }
            fk|=(rk & lsn);
            indexType gkmer=fk * i_nmasks + mask;
            i_mlmerGkmer[i] = gkmer;
        }
    }
//     dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "\t(" << dur.count() << ")" << endl;
    

    if(gLmerCounter<T>::i_doubleStrand){
        //MAPPA gkmers -> dsgkmers
//         cout << "building gkmers -> dsgkmers...";cout.flush();
        prime_time = std::chrono::high_resolution_clock::now();
        {
            arma::Col<indexType> list1(i_ngkmers,arma::fill::zeros);
            arma::Col<indexType> list2(i_ngkmers,arma::fill::zeros);

            for(indexType gkmer=0;gkmer<i_ngkmers;gkmer++){
                indexType kmer=gkmer/i_nmasks;
                indexType rckmer=0;
                for(indexType h=0;h<(i_k-1);h++){
                    rckmer|= ( (~kmer) & 3);
                        kmer>>=2;
                        rckmer<<=2;
                }
                rckmer|= ( (~kmer) & 3);
                
                
                indexType nmask=gkmer % i_nmasks;
                indexType mask=i_masks[nmask];
                indexType rmask=0;
                for(indexType h=0;h<gLmerCounter<T>::i_l-1;h++){
                    rmask |= mask & 1;
                    mask >>= 1;
                    rmask <<= 1;
                }
                rmask |= mask & 1;
                indexType rnmask=i_mask_nmask[rmask];
                indexType rgkmer = rckmer * i_nmasks + rnmask;
                
                indexType b=min(gkmer,rgkmer);
                list1[gkmer]=b;
                list2[b]=1;
                
                //cout << gkmer << "\t" << gkmer/i_nmasks << "\t" << nmask << "\t" << i_masks[nmask] << "\t" << rgkmer << "\t" << rckmer << "\t" << rnmask << "\t" << rmask << endl;
    //             cout << printgkmer(gkmer,i_l,i_masks) << "\t" << printgkmer(rgkmer,i_l,i_masks) << endl;
                
            }
            
            arma::uvec pos=arma::find(list2);
            arma::Col<indexType> i_values=list1(pos);
            list2.elem(pos)=arma::regspace<arma::Col<indexType>>(0,pos.n_elem);
            arma::Col<indexType> i_map=list2(list1);
            i_gkmerDSGkmer = vector<indexType>(i_map.begin(),i_map.end());
            i_DSgkmers=vector<indexType>(i_values.begin(),i_values.end());
            i_nDSgkmers=i_values.size();
        }    
//         dur=std::chrono::high_resolution_clock::now()-prime_time;
//         cout << "\t(" << dur.count() << ")" << endl;
    }else{
        i_nDSgkmers=0;
    }
//     cout << "N. double strand gapped kmers: \t" << i_nDSgkmers << endl;
}

template<typename T>
gMDense<T> gGappedKmerCounter<T>::countDense(const std::vector<std::string> & sequences,bool strict) const{
    indexType i_mask=pow(4,gLmerCounter<T>::i_l)-1;
    indexType nrows=gLmerCounter<T>::i_doubleStrand?i_nDSgkmers:i_ngkmers;
    
    auto prime_time = std::chrono::high_resolution_clock::now();
    vector<T> gkmers(sequences.size() * nrows);
    gMThreadException te;
#pragma omp parallel for
    for(indexType s=0;s<sequences.size();s++){
        if(te.exceptionOccurred()) continue;
        try{
            indexType mlmer;
            indexType base=s*nrows;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                if(u->second.size()>=gLmerCounter<T>::i_l){
                    nnzero=true;
                    auto bb=u->second.begin();
                    indexType lmer=*bb++;
                    for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    if(gLmerCounter<T>::i_doubleStrand){
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            mlmer=(lmer & i_mask) * i_nmasks;
                            for(indexType mask=0;mask<i_nmasks;mask++){
                                gkmers[base+i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                            }
                        }
                    }else{
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            mlmer=(lmer & i_mask) * i_nmasks;
                            for(indexType mask=0;mask<i_nmasks;mask++){
                                gkmers[base+i_mlmerGkmer[mlmer++]]++;
                            }
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gGappedKmerCounter: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "time: " << dur.count()*1000 << "ms" << endl;
    te.Rethrow();
    return move(gMDense<T>(nrows,sequences.size(),gkmers.data()));
}

template<typename T>
gMSparse<T> gGappedKmerCounter<T>::countSparse(const std::vector<std::string> & sequences,bool strict) const{
    gMThreadException te;
    indexType i_mask=pow(4,gLmerCounter<T>::i_l)-1;
    indexType nrows=gLmerCounter<T>::i_doubleStrand?i_nDSgkmers:i_ngkmers;
    auto prime_time = std::chrono::high_resolution_clock::now();    
    vector<map<indexType,T>> gkmers(sequences.size());
#pragma omp parallel for
    for(indexType s=0;s<sequences.size();s++){
        if(te.exceptionOccurred()) continue;
        try{
            map<indexType,T> tres;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                auto bb=u->second.begin();
                if(u->second.size()>=gLmerCounter<T>::i_l){
                    nnzero=true;
                    indexType lmer=*bb;
                    bb++;
                    for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    if(gLmerCounter<T>::i_doubleStrand){
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            
                            indexType mlmer=(lmer & i_mask) * i_nmasks;
                            for(indexType mask=0;mask<i_nmasks;mask++){
                                tres[i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                            }
                        }
                    }else{
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            
                            indexType mlmer=(lmer & i_mask) * i_nmasks;
                            for(indexType mask=0;mask<i_nmasks;mask++){
                                tres[i_mlmerGkmer[mlmer++]]++;
                            }
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gGappedKmerCounter: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            gkmers[s]=tres;
            //checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "time: " << dur.count()*1000 << "ms" << endl;
    te.Rethrow();
    return move(gMSparse<T>(gkmers,nrows));
}

template<typename T>
gMDense<T> gGappedKmerCounter<T>::countWinDense(const string & sequence,indexType win,bool strict) const{
    gMDense<T> RES;    
    indexType i_mask=pow(4,gLmerCounter<T>::i_l)-1;
    indexType nrows=gLmerCounter<T>::i_doubleStrand?i_nDSgkmers:i_ngkmers;
    indexType ncols=sequence.size()-win+1;
    
    auto prime_time = std::chrono::high_resolution_clock::now();    
    vector<T> gkmers(ncols * nrows);
    vector<unsigned char> tres=seqCode(sequence);
    
    if(find(tres.begin(),tres.end(),4)==tres.end()){
        queue<indexType> qw;
        indexType mlmer;
        indexType base=0;
        indexType base2=0;

        auto bb=tres.begin();
        indexType lmer=*bb++;
        for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
            lmer <<= 2;
            lmer|=*bb++;
        }

        if(gLmerCounter<T>::i_doubleStrand){
            for(indexType pos=gLmerCounter<T>::i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    gkmers[base2+i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                }
            }
            while(base<ncols-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    gkmers[base2+i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                }                
                copy_n(gkmers.begin()+base2,nrows,gkmers.begin()+(base+1)*nrows);
                base2 = (++base)*nrows;
                mlmer=qw.front();
                for(indexType mask=0;mask<i_nmasks;mask++){
                    gkmers[base2+i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]--;
                }                
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            mlmer=(lmer & i_mask) * i_nmasks;
            for(indexType mask=0;mask<i_nmasks;mask++){
                gkmers[base2+i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
            }
        }else{
            for(indexType pos=gLmerCounter<T>::i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    gkmers[base2+i_mlmerGkmer[mlmer++]]++;
                }
            }
            while(base<ncols-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    gkmers[base2+i_mlmerGkmer[mlmer++]]++;
                }
                copy_n(gkmers.begin()+base*nrows,nrows,gkmers.begin()+(base+1)*nrows);
                base2=(++base) * nrows;
                mlmer=qw.front();
                for(indexType mask=0;mask<i_nmasks;mask++){
                    gkmers[base2+i_mlmerGkmer[mlmer++]]--;
                }                
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            mlmer=(lmer & i_mask) * i_nmasks;
            for(indexType mask=0;mask<i_nmasks;mask++){
                gkmers[base2+i_mlmerGkmer[mlmer++]]++;
            }
        }
        std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
//         cout << "time: " << dur.count()*1000 << "ms" << endl;
        RES=gMDense<T>(nrows,ncols,gkmers.data());
    }else{ //Invalid charaters we do it thes slow way
        cout << "The slow way" << endl;
        vector<string> seqs(ncols);
        for(indexType p=0;p<ncols;p++){
            seqs[p] = sequence.substr(p,win);
        }
        RES=countDense(seqs,strict);
    }

    return move(RES);
}

template<typename T>
gMSparse<T> gGappedKmerCounter<T>::countWinSparse(const string & sequence,indexType win,bool strict) const{
    gMSparse<T> RES;    
    indexType i_mask=pow(4,gLmerCounter<T>::i_l)-1;
    indexType nrows=gLmerCounter<T>::i_doubleStrand?i_nDSgkmers:i_ngkmers;
    indexType nwins=sequence.size()-win+1;
    
    auto prime_time = std::chrono::high_resolution_clock::now();
    vector<map<indexType,T>> gkmers(nwins);
    vector<unsigned char> tres=seqCode(sequence);
    
    if(find(tres.begin(),tres.end(),4)==tres.end()){
        map<indexType,T> ores;
        queue<indexType> qw;
        indexType mlmer;
        indexType base=0;

        auto bb=tres.begin();
        indexType lmer=*bb++;
        for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
            lmer <<= 2;
            lmer|=*bb++;
        }
        if(gLmerCounter<T>::i_doubleStrand){
            for(indexType pos=gLmerCounter<T>::i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    ores[i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                }
            }
            while(base<nwins-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    ores[i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                }
                gkmers[base]=ores;
                base++;
                mlmer=qw.front();
                for(indexType mask=0;mask<i_nmasks;mask++){
                    indexType i=i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]];
                    if((--ores[i]) == 0){
                        ores.erase(i);
                    }
                }                
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            mlmer=(lmer & i_mask) * i_nmasks;
            for(indexType mask=0;mask<i_nmasks;mask++){
                ores[i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
            }
            gkmers[base]=ores;
        }else{
            for(indexType pos=gLmerCounter<T>::i_l-1;pos<win-1;pos++){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    ores[i_mlmerGkmer[mlmer++]]++;
                }
            }
            while(base<nwins-1){
                lmer <<= 2;
                lmer|=*bb++;
                mlmer=(lmer & i_mask) * i_nmasks;
                qw.push(mlmer);
                for(indexType mask=0;mask<i_nmasks;mask++){
                    ores[i_mlmerGkmer[mlmer++]]++;
                }
                gkmers[base]=ores;
                base++;
                mlmer=qw.front();
                for(indexType mask=0;mask<i_nmasks;mask++){
                    indexType i=i_mlmerGkmer[mlmer++];
                    if((--ores[i]) == 0){
                        ores.erase(i);
                    }
                }                
                qw.pop();
            }
            lmer <<= 2;
            lmer|=*bb++;
            mlmer=(lmer & i_mask) * i_nmasks;
            for(indexType mask=0;mask<i_nmasks;mask++){
                ores[i_mlmerGkmer[mlmer++]]++;
            }
            gkmers[base]=ores;
        }
        std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
//         cout << "time: " << dur.count()*1000 << "ms" << endl;
        RES=gMSparse<T>(gkmers,nrows);
    }else{ //Invalid charaters we do it the slow way
        cout << "The slow way" << endl;
        vector<string> sequences(nwins);
        for(indexType p=0;p<nwins;p++){
            sequences[p] = sequence.substr(p,win);
        }
        vector<map<indexType,T>> gkmers(nwins);
        for(indexType s=0;s<sequences.size();s++){

            map<indexType,T> tres;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                auto bb=u->second.begin();
                if(u->second.size()>=gLmerCounter<T>::i_l){
                    nnzero=true;
                    indexType lmer=*bb;
                    bb++;
                    for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    if(gLmerCounter<T>::i_doubleStrand){
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            
                            indexType mlmer=(lmer & i_mask) * i_nmasks;
                            for(indexType mask=0;mask<i_nmasks;mask++){
                                tres[i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]++;
                            }
                        }
                    }else{
                        while(bb!=u->second.end()){
                            lmer <<= 2;
                            lmer|=*bb++;
                            
                            indexType mlmer=(lmer & i_mask) * i_nmasks;
                            for(indexType mask=0;mask<i_nmasks;mask++){
                                tres[i_mlmerGkmer[mlmer++]]++;
                            }
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gGappedKmerCounter: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            gkmers[s]=tres;
        }
        RES=gMSparse<T>(gkmers,nrows); 
        cout << "done" << endl;
    }
    return move(RES);
}

template<typename T>
gMDense<T> gGappedKmerCounter<T>::countPWM(const geco::kmers::gPWMSet & PWMs,indexType padSize) const{
    if(padSize>=gLmerCounter<T>::i_l){
        throw gException("countPWM: padSize must be lower than kmer length");
    }
    
//     gSize padSize=gLmerCounter<T>::i_l-1;
    indexType i_mask=pow(4,gLmerCounter<T>::i_l)-1;
    indexType nrows=gLmerCounter<T>::i_doubleStrand?i_nDSgkmers:i_ngkmers;
    vector<T> gkmers(PWMs.size() * nrows);
    
	gMThreadException te;
#pragma omp parallel for
	for(indexType a=0;a<PWMs.size();a++){
        indexType base = a * nrows;        
        if(te.exceptionOccurred()) continue;
		try{
			auto m=PWMs.begin();
			advance(m,a);
			
			gArray<T> background=m->second.getBackground();
			gMatrix<T> pwm=m->second.getPWM();
			
			gArray<T> pads;
			for(indexType i=0;i<padSize;i++){
				pads.concatenate(background);
			}
		
			gArray<T> full=pads;
			for(indexType i=0;i<pwm.getColsNum();i++){
				full.concatenate(pwm.getCol(i));
			}
			full.concatenate(pads);
			gMatrix<T> padded(2*padSize+pwm.getColsNum(),4,full);
			
			const T * apads=padded.getConstData();
			indexType nstep=padded.getRowsNum()-gLmerCounter<T>::i_l+1;
			indexType cols[gLmerCounter<T>::i_l];
			
			for(indexType lmer=0;lmer<gLmerCounter<T>::i_nlmers;lmer++){
				indexType b=lmer;
				for(indexType i=0;i<gLmerCounter<T>::i_l;i++){
					cols[gLmerCounter<T>::i_l-i-1]=b & 0x00000003;
					b >>= 2;
				}
				T ptot=0;
				for(indexType p=0;p<nstep;p++){
					T prob=1;
					for(indexType i=0;i<gLmerCounter<T>::i_l;i++){
						prob*=apads[4*(p+i)+cols[i]];
					}
					ptot+=prob;
				}
				ptot/=(nstep * i_nmasks);
				indexType mlmer = lmer * i_nmasks;
				if(gLmerCounter<T>::i_doubleStrand){
                    for(indexType mask=0;mask<i_nmasks;mask++){
                        gkmers[base + i_gkmerDSGkmer[i_mlmerGkmer[mlmer++]]]+=ptot;
                    }
                }else{
                    for(indexType mask=0;mask<i_nmasks;mask++){
                        gkmers[base + i_mlmerGkmer[mlmer++]]+=ptot;
                    }                   
                }
				
			}
		}catch(...){
			te.CaptureException();
		}
	}
	te.Rethrow();
	//replace_if(gkmers.begin(), gkmers.end(),isSmall<T>(nrows),0.0);
    return move(gMDense<T>(nrows,PWMs.size(),gkmers.data()));
}

template<typename T>
gMDense<T> gGappedKmerCounter<T>::getPWM(const gMDense<T> & W,indexType pwmlength,const vector<double> & pwm0,double tol,indexType maxiter) const{
    nlopt::opt myopt(nlopt::LN_COBYLA,pwmlength*3);
    arma::Col<double> w=arma::conv_to<arma::Mat<double>>::from(armadillo::convert<gMDense,arma::Mat,T>(W)).col(0);
    pwmdataG data(gLmerCounter<T>::i_l,gLmerCounter<T>::i_nlmers,i_nDSgkmers,i_ngkmers,i_nmasks,pwmlength,w,pwm0,i_gkmerDSGkmer,i_mlmerGkmer,myopt);
//     data.i_outer.Start([=]{ std::cout << "\toptimizing...";});
    
    myopt.set_max_objective(pwmfuncG, &data);
    myopt.set_stopval(1-tol);
    myopt.set_ftol_abs(tol);
    myopt.set_xtol_abs(tol);
//     myopt.set_maxtime(120);
    myopt.set_maxeval(maxiter);
    myopt.set_lower_bounds(0.0);
    myopt.set_upper_bounds(1.0);
    myopt.add_inequality_mconstraint(pwmConstraintG, NULL, vector<double>(pwmlength,1e-8));
   
    double minf;
    arma::Mat<double> M0(pwm0.data(),4,data.i_ncols);
    arma::Mat<double> M1=M0.rows(data.R);
    vector<double> v(M1.begin(),M1.end());
//     data.i_inner.Start([=]{ std::cout << "\titeration: " << data.ncall+1;});
    try{
        myopt.optimize(v, minf);
    }catch(nlopt::roundoff_limited & e){
        warning("Roundoff limited: resuts could be inacccurate");
    }catch(nlopt::forced_stop & e){
        warning("user interruption: results cuold be still valid");
    }catch(std::runtime_error &e ){
        throw gException("something bad");
    }
//     data.i_inner.End([=]{ std::cout << "\titeration: " << data.ncall;});
    arma::Mat<double> M2(v.data(),3,data.i_ncols);
    arma::Mat<double> M(4,data.i_ncols);
    M.rows(data.R)=M2;
    M.row(3)=1-arma::sum(M2);
    arma::Col<indexType> bad=find(M.row(3)<0);
    if(bad.n_elem>0){
        for(indexType i=0;i<bad.n_elem;i++){
            arma::Col<double> c=M.col(bad(i));
            M.col(bad(i))=(c-c(3))/sum(c-c(3));
        }
        warning("Negative entries in pwm due to approximation in constraints: corrected");        
    }
    vector<T> ret(4*pwmlength);
    copy(M.memptr(),M.memptr()+4*pwmlength,ret.begin());
    //cout << myopt.get_xtol_abs() << endl;
    cout << myopt.get_ftol_abs() << endl;
//     data.i_outer.End([=]{ std::cout << " done";});
    return move(armadillo::convert<arma::Mat,gMDense,T>( arma::conv_to<arma::Mat<T>>::from(M)));
}


template<typename T>
pair< vector<string>,vector<string> > gGappedKmerCounter<T>::getStrings() const{
    pair< vector<string>,vector<string> > ret;    
    if(gLmerCounter<T>::i_doubleStrand){
        vector<string> falpha={"A","C","G","T"};
        vector<string> ralpha={"T","G","C","A"};                
        ret.first=vector<string>(i_nDSgkmers);
        ret.second=vector<string>(i_nDSgkmers);
#pragma omp parallel for
        for(indexType gkmer=0;gkmer<i_nDSgkmers;gkmer++){
            indexType gk=i_DSgkmers[gkmer];
            indexType kmer=gk/i_nmasks;
            indexType mask=i_masks[gk % i_nmasks];
            for(indexType i=0;i<gLmerCounter<T>::i_l;i++){
                if(mask & 1){
                    ret.first[gkmer]=falpha[kmer & 3] + ret.first[gkmer];
                    ret.second[gkmer]=ret.second[gkmer] + ralpha[kmer & 3];
                    kmer >>= 2;
                }else{
                    ret.first[gkmer]=string("N") + ret.first[gkmer];
                    ret.second[gkmer]=ret.second[gkmer] + string("N");
                }
                mask >>= 1;
            }
        }
    }else{
        vector<string> alpha={"A","C","G","T"};        
        ret.first=vector<string>(i_ngkmers);
        ret.second=vector<string>(i_ngkmers);
#pragma omp parallel for
        for(indexType gk=0;gk<i_ngkmers;gk++){
            indexType kmer=gk/i_nmasks;
            indexType mask=i_masks[gk % i_nmasks];
            for(indexType i=0;i<gLmerCounter<T>::i_l;i++){
                if(mask & 1){
                    ret.first[gk]=alpha[kmer & 3] + ret.first[gk];
                    kmer >>= 2;
                }else{
                    ret.first[gk]=string("N") + ret.first[gk];
                }
                mask >>= 1;
            }
        }
    }
    return ret;
}

template class geco::methods::kmers::gGappedKmerCounter<float>;
template class geco::methods::kmers::gGappedKmerCounter<double>;
template class geco::methods::kmers::gGappedKmerCounter<indexType>;
template class geco::methods::kmers::gGappedKmerCounter<unsigned short>;




template<typename T>
gLmerEstimator<T>::gLmerEstimator(geco::methods::indexType l,geco::methods::indexType k,bool doubleStrand):gGappedKmerCounter<T>(l,k,doubleStrand),i_mmcoeffs(gGappedKmerCounter<T>::i_nkmers),i_mm(gGappedKmerCounter<T>::i_nkmers){
    auto prime_time = std::chrono::high_resolution_clock::now();    
    std::chrono::duration<double> dur;
    
    //MAPPA W dipendente da nmismatch
//     cout << "building W map...";cout.flush();
    vector<T> tmpmmcoeffs(gGappedKmerCounter<T>::i_k+1);
    prime_time = std::chrono::high_resolution_clock::now();
    {
        int64_t kml = (int64_t) k - (int64_t) l;
        for(indexType nmm=0;nmm<=gGappedKmerCounter<T>::i_k;nmm++){
            double s=0;
            for(indexType n=0;n<=(k-nmm);n++){
                s+=choose(l,n) * pow((double)3,(double)n);
            }
            tmpmmcoeffs[nmm] = s * choose(kml,nmm) / ((double) pow((double)4,(double)l) * choose(l,k) * choose(k,nmm) );
        }
    }
//     dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "\t(" << dur.count() << ")" << endl;


    //MAPPA2 W dipendente da nmismatch
//     cout << "building W map2...";cout.flush();
    prime_time = std::chrono::high_resolution_clock::now();
    {
        for(indexType x=0;x<gGappedKmerCounter<T>::i_nkmers;x++){
            indexType c=~x;
            indexType counts=0;
            for(indexType i=0;i<gGappedKmerCounter<T>::i_k;i++){
                counts+=((c & 3)!=3);
                c>>=2;
            }
            i_mmcoeffs[x]=tmpmmcoeffs[counts];
            i_mm[x]=counts;
        }
    }
//     dur=std::chrono::high_resolution_clock::now()-prime_time;
//     cout << "\t(" << dur.count() << ")" << endl;
    if(gLmerCounter<T>::i_doubleStrand){    
        //i_gkmersNDSgkmers; //tells how many gkmers contribute ti a DSgkmers
        for(indexType gkmer=0;gkmer<gGappedKmerCounter<T>::i_ngkmers;gkmer++){
            i_gkmersNDSgkmers[gGappedKmerCounter<T>::i_gkmerDSGkmer[gkmer]]++;
        }
    }
    
    
    

#ifdef GECO_HAS_CUDA
    i_cudaObject=nullptr;
    if(useCuda()){
        i_cudaObject=new gCudaCounter<T>(*this);
    }else{
        i_cudaObject=nullptr;
    }
#endif
}

#ifdef GECO_HAS_CUDA
template<typename T>
gLmerEstimator<T>::~gLmerEstimator(){
    if(useCuda()){
        if(i_cudaObject!=nullptr){
            delete ((gCudaCounter<T> *) i_cudaObject);
        }
    }
}
#endif

template<typename T>
geco::methods::gMDense<T> gLmerEstimator<T>::countDense(const std::vector<string> & sequences,bool strict) const{
    indexType i_mask=pow(4,gGappedKmerCounter<T>::i_l)-1;
    indexType nrows=gGappedKmerCounter<T>::i_doubleStrand?gGappedKmerCounter<T>::i_nDSlmers:gGappedKmerCounter<T>::i_nlmers;    
    vector<T> lmers(sequences.size() * nrows);
    
    
#ifdef GECO_HAS_CUDA
    const gCudaCounter<T> & ccounter = *((gCudaCounter<T> *) i_cudaObject);
    if(useCuda()){
        gMThreadException te;
//#pragma omp parallel for
        for(indexType s=0;s<sequences.size();s++){
            if(te.exceptionOccurred()) continue;
            try{
                indexType base=s*nrows;
                bool nnzero=ccounter.estimatelmers(seqCode(sequences[s],s,strict),lmers,base,gGappedKmerCounter<T>::i_doubleStrand,strict);
                if(!nnzero) warning(string("gLmerEstimator: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
                checkInterrupt();
            }catch(...){
                te.CaptureException();
            }
        }
        te.Rethrow();
        return move(gMDense<T>(nrows,sequences.size(),lmers.data()));
    }
#endif

    T eps=std::numeric_limits<T>::epsilon();
    gMThreadException te;
#pragma omp parallel for
    for(indexType s=0;s<sequences.size();s++){
        if(te.exceptionOccurred()) continue;
        try{
            indexType base=s*nrows;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                auto bb=u->second.begin();
                if(u->second.size()>=gLmerCounter<T>::i_l){
                    nnzero=true;
                    indexType lmer=*bb++;
                    //bb++;
                    for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    map<indexType,indexType> gkmers;
                    while(bb!=u->second.end()){
                        lmer <<= 2;
                        lmer|=*bb++;
                        
                        indexType mlmer=(lmer & i_mask) * gGappedKmerCounter<T>::i_nmasks;
                        for(indexType mask=0;mask<gGappedKmerCounter<T>::i_nmasks;mask++){
                            gkmers[gGappedKmerCounter<T>::i_mlmerGkmer[mlmer++]]++;
                        }
                    }
                    
                    if(gLmerCounter<T>::i_doubleStrand){
                        for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                            T tot=0;                            
                            for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                                tot+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                            }
                            if(tot>eps){
                                lmers[base+gGappedKmerCounter<T>::i_lmerDSLmer[lmer]]=tot;
                            }
                            
                            
                        }
                    }else{
                        for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                            T tot=0;
                            indexType gkmer_kmer,lmer_kmer;
                            for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                                gkmer_kmer = gkmer->first / gGappedKmerCounter<T>::i_nmasks;
                                lmer_kmer = gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks;
                                tot+=i_mmcoeffs[gkmer_kmer ^ lmer_kmer] * gkmer->second;
                            }
                            if(tot>eps){
                                lmers[base+lmer]=tot;
                            }
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gLmerEstimator: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return move(gMDense<T>(nrows,sequences.size(),lmers.data()));
}

template<typename T>
geco::methods::gMSparse<T> gLmerEstimator<T>::countSparse(const std::vector<string> & sequences,bool strict) const{
    indexType i_mask=pow(4,gGappedKmerCounter<T>::i_l)-1;
    indexType nrows=gGappedKmerCounter<T>::i_doubleStrand?gGappedKmerCounter<T>::i_nDSlmers:gGappedKmerCounter<T>::i_nlmers;
    vector<map<indexType,T>> lmers(sequences.size());
    
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        const gCudaCounter<T> & ccounter = *((gCudaCounter<T> *) i_cudaObject);
        
        gMThreadException te;
#pragma omp parallel for
        for(indexType s=0;s<sequences.size();s++){
            if(te.exceptionOccurred()) continue;
            try{
                vector<T> tmp(nrows);
                bool nnzero=ccounter.estimatelmers(seqCode(sequences[s],s,strict),tmp,0,gGappedKmerCounter<T>::i_doubleStrand);
                if(!nnzero) warning(string("gLmerEstimator: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
                for(indexType i=0;i<nrows;i++){
                    if(tmp[i]>0){
                        lmers[s][i]=tmp[i];
                    }
                }
                checkInterrupt();
            }catch(...){
                te.CaptureException();
            }
        }
        te.Rethrow();
        return move(gMSparse<T>(lmers,nrows));
    }
#endif
    T eps=std::numeric_limits<T>::epsilon();
    
    gMThreadException te;
#pragma omp parallel for
    for(indexType s=0;s<sequences.size();s++){
        if(te.exceptionOccurred()) continue;
        try{
            map<indexType,T> tmers;
            indexType base=s*nrows;
            map<size_t,list<unsigned char>> nseq = seqCode(sequences[s],s,strict);
            bool nnzero=false;
            for(auto u=nseq.begin();u!=nseq.end();u++){
                auto bb=u->second.begin();            
                if(u->second.size()>=gLmerCounter<T>::i_l){
                    nnzero=true;
                    indexType lmer=*bb++;
                    for(indexType i=1;i<gLmerCounter<T>::i_l-1;i++){
                        lmer <<= 2;
                        lmer|=*bb++;
                    }
                    
                    map<indexType,indexType> gkmers;
                    while(bb!=u->second.end()){
                        lmer <<= 2;
                        lmer|=*bb++;
                        
                        indexType mlmer=(lmer & i_mask) * gGappedKmerCounter<T>::i_nmasks;
                        for(indexType mask=0;mask<gGappedKmerCounter<T>::i_nmasks;mask++){
                            gkmers[gGappedKmerCounter<T>::i_mlmerGkmer[mlmer++]]++;
                        }
                    }
                    
                    if(gGappedKmerCounter<T>::i_doubleStrand){
                        vector<T> tmplmers(gGappedKmerCounter<T>::i_nDSlmers,0);
                        for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                            indexType ss = lmer * gGappedKmerCounter<T>::i_nmasks;
                            for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                                tmplmers[gGappedKmerCounter<T>::i_lmerDSLmer[lmer]]+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[ss + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                            }
                        }
                        for(indexType dslmer=0;dslmer<gGappedKmerCounter<T>::i_nDSlmers;dslmer++){
                            if(tmplmers[dslmer]>eps){
                                tmers[dslmer]=tmplmers[dslmer];
                            }
                        }
                    }else{
                        for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                            indexType ss=lmer * gGappedKmerCounter<T>::i_nmasks;
                            T tot=0;
                            for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                                tot+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[ss + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                            }
                            if(tot>eps){
                                tmers[lmer]=tot;
                            }
                        }
                    }
                }
            }
            if(!nnzero) warning(string("gLmerEstimator: sequence ")+sequences[s]+" ("+to_string(s+1)+") "+" had zero counts");
            lmers[s]=tmers;
            checkInterrupt();
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return move(gMSparse<T>(lmers,nrows));
}

template<typename T>
gMDense<T> gLmerEstimator<T>::countDense(const gMDense<T> & gkmerCounts) const{
    indexType nrows=gGappedKmerCounter<T>::i_doubleStrand?gGappedKmerCounter<T>::i_nDSlmers:gGappedKmerCounter<T>::i_nlmers;    
    vector<T> lmers(gkmerCounts.i_ncols * nrows);
    arma::Mat<T> gkm=arma::square(armadillo::convert<gMDense,arma::Mat,T>(gkmerCounts));
    T eps=std::numeric_limits<T>::epsilon();
    
#ifdef GECO_HAS_CUDA
    const gCudaCounter<T> & ccounter = *((gCudaCounter<T> *) i_cudaObject);
    if(useCuda()){
        gMThreadException te;
#pragma omp parallel for
        for(indexType s=0;s<gkmerCounts.i_ncols;s++){
            if(te.exceptionOccurred()) continue;
// #pragma omp critical
//             {
//                 cout << s << endl;
//             }
        
            try{
                indexType base=s*nrows;

                if(gLmerCounter<T>::i_doubleStrand){
                    map<indexType,T> gkmers;
                    for(indexType gkmer=0;gkmer<gGappedKmerCounter<T>::i_ngkmers;gkmer++){
                        indexType dsgkmer=gGappedKmerCounter<T>::i_gkmerDSGkmer[gkmer];
                        T val=gkm(dsgkmer,s);
                        if(val>eps){
                            gkmers[gkmer]=val / i_gkmersNDSgkmers.at(dsgkmer);
                        }
                    }
                    ccounter.estimatelmers(gkmers,lmers,base,gGappedKmerCounter<T>::i_doubleStrand);
                }else{
                    map<indexType,T> gkmers;
                    for(indexType gkmer=0;gkmer<gGappedKmerCounter<T>::i_ngkmers;gkmer++){
                        T val=gkm(gkmer,s);
                        if(val>eps){
                            gkmers[gkmer]=val;
                        }
                    }
                    ccounter.estimatelmers(gkmers,lmers,base,gGappedKmerCounter<T>::i_doubleStrand);
                }
            }catch(...){
                te.CaptureException();
            }
        }
        te.Rethrow();
        return move(gMDense<T>(nrows,gkmerCounts.i_ncols,lmers.data()));
    }
#endif
    

    gMThreadException te;

#pragma omp parallel for
    for(indexType s=0;s<gkmerCounts.i_ncols;s++){
        if(te.exceptionOccurred()) continue;

// #pragma omp critical                                    
//         {
//             cout << s << endl;
//         }
        
        try{
            indexType base=s*nrows;

            if(gLmerCounter<T>::i_doubleStrand){
                map<indexType,T> gkmers;
                for(indexType gkmer=0;gkmer<gGappedKmerCounter<T>::i_ngkmers;gkmer++){
                    indexType dsgkmer=gGappedKmerCounter<T>::i_gkmerDSGkmer[gkmer];
                    T val=gkm(dsgkmer,s);
                    if(val>eps){
                        gkmers[gkmer]=val / i_gkmersNDSgkmers.at(dsgkmer);
                    }
                }
                
                indexType cc=0;
                for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                    indexType ss=lmer * gGappedKmerCounter<T>::i_nmasks;
                    for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                        lmers[base+gGappedKmerCounter<T>::i_lmerDSLmer[lmer]]+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[ss + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                    }
                    if(cc % 1000 ==0){
                        cc=0;
                        if(te.exceptionOccurred()){
                            break;
                        }
                        checkInterrupt();
                    }
                    cc++;
                }
                for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nDSlmers;lmer++){
                    if(lmers[lmer]<=eps){
                        lmers[lmer]=0;
                    }
                }
                
            }else{
                map<indexType,T> gkmers;
                for(indexType gkmer=0;gkmer<gGappedKmerCounter<T>::i_ngkmers;gkmer++){
                    T val=gkm(gkmer,s);
                    if(val>eps){
                        gkmers[gkmer]+=val;
                    }
                }
              
                indexType cc=0;                
                for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                    indexType ss=lmer * gGappedKmerCounter<T>::i_nmasks;
                    T tot=0;
                    for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                        tot+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[ss + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                    }
                    if(tot>eps){
                        lmers[base+lmer]=tot;
                    }
                    if(cc % 1000 ==0){
                        cc=0;
                        if(te.exceptionOccurred()){
                            break;
                        }
                        checkInterrupt();
                    }
                    cc++;
                }
            }
        }catch(...){
            te.CaptureException();
        }
    }
    te.Rethrow();
    return move(gMDense<T>(nrows,gkmerCounts.i_ncols,lmers.data()));
}



template<typename T>
gMSparse<T> gLmerEstimator<T>::calculateG(){
    cout << "CALCULATE G" << endl;
    indexType nrows=gGappedKmerCounter<T>::i_doubleStrand?gGappedKmerCounter<T>::i_nDSlmers:gGappedKmerCounter<T>::i_nlmers;
    vector<map<indexType,T>> lmers(gGappedKmerCounter<T>::i_nlmers);
    if(gGappedKmerCounter<T>::i_doubleStrand){
#pragma omp parallel for
        for(indexType rlmer=0;rlmer<gGappedKmerCounter<T>::i_nlmers;rlmer++){
            cout << "dsrlmer=" <<  rlmer;cout.flush();
            auto prime_time = std::chrono::high_resolution_clock::now();    
            map<indexType,indexType> gkmers;
            indexType mlmer=rlmer  * gGappedKmerCounter<T>::i_nmasks;
            for(indexType mask=0;mask<gGappedKmerCounter<T>::i_nmasks;mask++){
                gkmers[gGappedKmerCounter<T>::i_mlmerGkmer[mlmer++]]++;
            }
            T eps=10*std::numeric_limits<T>::epsilon();
            for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                T tot=0;
                for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                    tot+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                }
                if(tot>eps){
                    lmers[rlmer][gGappedKmerCounter<T>::i_lmerDSLmer[lmer]]=tot;
                }
            }
            std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "\t(" << dur.count() << ")" << endl;
        }
    }else{
//#pragma omp parallel for
        for(indexType rlmer=0;rlmer<gGappedKmerCounter<T>::i_nlmers;rlmer++){
            cout << "ssrlmer=" <<  rlmer;cout.flush();
            auto prime_time = std::chrono::high_resolution_clock::now();    
            
            map<indexType,indexType> gkmers;
            indexType mlmer=rlmer  * gGappedKmerCounter<T>::i_nmasks;
            for(indexType mask=0;mask<gGappedKmerCounter<T>::i_nmasks;mask++){
                gkmers[gGappedKmerCounter<T>::i_mlmerGkmer[mlmer++]]++;
            }
            T eps=std::numeric_limits<T>::epsilon();
            vector<vector<map<indexType,T>> > lmers2(gGappedKmerCounter<T>::i_k+1,vector<map<indexType,T>>(gGappedKmerCounter<T>::i_nlmers));
            for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
                T tot=0;
                for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
                    lmers2[i_mm[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)]][rlmer][lmer]++;
                    //tot+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                    //lmers[rlmer][lmer]+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
                }
                
//                 if(abs(tot)>eps){
// //                     cout << tot << "\t" << eps << "\tOK" << endl;                
//                     lmers[rlmer][lmer]=tot;
//                     
//                 }else{
// //                     cout << tot << "\t" << eps << endl;                
//                 }
            }
            
//             for(indexType lmer=0;lmer<gGappedKmerCounter<T>::i_nlmers;lmer++){
//                 T tot=0;
//                 for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
//                     tot+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
//                     //lmers[rlmer][lmer]+=i_mmcoeffs[(gkmer->first / gGappedKmerCounter<T>::i_nmasks) ^ (gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer->first % gGappedKmerCounter<T>::i_nmasks] / gGappedKmerCounter<T>::i_nmasks)] * gkmer->second;
//                 }
//                 
//                 if(abs(tot)>eps){
// //                     cout << tot << "\t" << eps << "\tOK" << endl;                
//                     lmers[rlmer][lmer]=tot;
//                     
//                 }else{
// //                     cout << tot << "\t" << eps << endl;                
//                 }
//             }
            std::chrono::duration<double> dur=std::chrono::high_resolution_clock::now()-prime_time;
            cout << "\t(" << dur.count() << ")" << endl;
        }
    }
    
    cout << "DONE" << endl;
    return move(gMSparse<T>(lmers,nrows));
}

template class geco::methods::kmers::gLmerEstimator<float>;
template class geco::methods::kmers::gLmerEstimator<double>;




template<typename T>
void hellinger_impl(gMDense<T> & counts){
#pragma omp parallel for
    for(indexType i=0;i<counts.i_ncols;i++){
        T s=accumulate((const T *)counts.i_values+i*counts.i_nrows,(const T *)counts.i_values+counts.i_nrows*(1+i),(T)0);
        transform((const T *)counts.i_values+i*counts.i_nrows,(const T *)counts.i_values+(i+1)*counts.i_nrows,(T *)counts.i_values+i*counts.i_nrows,[s](T n) -> T { return sqrt(n/s);});
    }
}

template<typename T>
void hellinger_impl(gMSparse<T> & counts){
#pragma omp parallel for
    for(indexType i=0;i<counts.i_ncols;i++){
        T s=accumulate((const T *)counts.i_cscValues+counts.i_cscColPtr[i],(const T *)counts.i_cscValues+counts.i_cscColPtr[i+1],(T)0);
        transform((const T *)counts.i_cscValues+counts.i_cscColPtr[i],(const T *)counts.i_cscValues+counts.i_cscColPtr[i+1],(T *)counts.i_cscValues+counts.i_cscColPtr[i],[s](T n) -> T { return sqrt(n/s);});
    }
}

template<template<typename> class D,typename T>
D<T> & geco::methods::kmers::toHellinger(D<T> & counts,bool use_cuda){
#ifdef GECO_HAS_CUDA
    if(use_cuda){
            hellinger_impl<T>(counts);
    }
#endif
    hellinger_impl<T>(counts);
    return counts;
}

template gMSparse<float> & geco::methods::kmers::toHellinger<gMSparse,float>(gMSparse<float> & counts,bool use_cuda);
template gMSparse<double> & geco::methods::kmers::toHellinger<gMSparse,double>(gMSparse<double> & counts,bool use_cuda);
template gMDense<float> & geco::methods::kmers::toHellinger<gMDense,float>(gMDense<float> & counts,bool use_cuda);
template gMDense<double> & geco::methods::kmers::toHellinger<gMDense,double>(gMDense<double> & counts,bool use_cuda);




