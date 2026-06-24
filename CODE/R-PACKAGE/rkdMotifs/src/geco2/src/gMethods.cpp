#include <list>

#include "gMethods_armadillo.h"
#include <omp.h>

using namespace geco::methods;
using namespace geco::methods::armadillo;

#ifdef GECO_HAS_CUDA
#include "gMethods_cuda.h"
using namespace geco::methods::cuda;
#endif

using namespace std;

class gMethodsConfigData{
public:
    gMethodsConfigData(){
        i_seedValue = 1;
        i_ncores=0;
        i_verbosity=0;
        arma::arma_rng::set_seed(i_seedValue); 
    }
    ~gMethodsConfigData(){
        if(i_ncores>0){
            omp_set_num_threads(i_ncores);
        }   
    }   

#ifdef GECO_HAS_CUDA
    bool                                        i_use_cuda=true;
#else
    bool                                        i_use_cuda=false;            
#endif
    int                                         i_ncores=0;
    unsigned int                                i_verbosity=1;
    gMethodInterruptFunction    i_interruptCallback=nullptr;
    gMethodWarningFunction      i_warningCallback=nullptr;
    
    unsigned int                i_seedValue;
};

static gMethodsConfigData __methodsConfig;
static list<string> __warnings;

void geco::methods::setNCores(int ncores){
    if(omp_in_parallel()==false){
        if(ncores>0){
            if(__methodsConfig.i_ncores==0){
                __methodsConfig.i_ncores=omp_get_max_threads();
            }
            omp_set_num_threads(ncores);
        }else{
            omp_set_num_threads(omp_get_max_threads());
            __methodsConfig.i_ncores=0;
        }
    }else{
        throw gMethodException("Cannot set n cores while in a parallel region.");
    }
}

void geco::methods::setgMethodsSeed(unsigned int seed){
    __methodsConfig.i_seedValue=seed;
    arma::arma_rng::set_seed(seed);
}

unsigned int geco::methods::getgMehtodsSeed(){
    return __methodsConfig.i_seedValue;
}

void geco::methods::setInterruptCallback(gMethodInterruptFunction function){
    __methodsConfig.i_interruptCallback=function;
}

void geco::methods::setWarningCallback(gMethodWarningFunction function){
    __methodsConfig.i_warningCallback=function;
}

void geco::methods::setVerbosity(unsigned short verbosity){
    __methodsConfig.i_verbosity=verbosity;
}

void geco::methods::setCudaUsage(bool use_cuda){
#ifdef GECO_HAS_CUDA
    __methodsConfig.i_use_cuda=use_cuda;
#else
    cout << "Setting useCuda to true when cuda is not available" << endl; 
    __methodsConfig.i_use_cuda=false;
#endif
}



bool geco::methods::useCuda(){
    return __methodsConfig.i_use_cuda;
}

unsigned short geco::methods::getVerbosity(){
    return __methodsConfig.i_verbosity;
}

void warningflush(){
    if(omp_get_thread_num()==0){
        if(__methodsConfig.i_warningCallback!=nullptr){
#pragma omp critical               
            {
                for(auto omsg=__warnings.begin();omsg!=__warnings.end();omsg++){
                    __methodsConfig.i_warningCallback(*omsg);
                }
                __warnings.clear();
            }
        }else{
#pragma omp critical               
            {
                for(auto omsg=__warnings.begin();omsg!=__warnings.end();omsg++){
                    cout << "Warning: " << *omsg << endl;
                }
                __warnings.clear();
            }
        }
    }    
}


void geco::methods::warning(const std::string & message){
    if(omp_get_thread_num()==0){
        warningflush();
        if(__methodsConfig.i_warningCallback!=nullptr){
#pragma omp critical               
            {
                __methodsConfig.i_warningCallback(message);
            }
        }else{
#pragma omp critical               
            {
                cout << "Warning: " << message << endl;
            }
        }
    }else{
#pragma omp critical
        {
            __warnings.push_back(message);
        }
    }
}

void geco::methods::checkInterrupt(int tn){
    if(__methodsConfig.i_interruptCallback!=nullptr){
        if(omp_get_thread_num()==tn){
            try{
                __methodsConfig.i_interruptCallback();
            }catch(...){
                throw gMethodsInterrupt();
            }
        }
    }
}


gMThreadException::gMThreadException():Ptr(nullptr){
}

gMThreadException::~gMThreadException(){
    warningflush();
}

void gMThreadException::Rethrow(){
    if(this->Ptr) std::rethrow_exception(this->Ptr);
}

void gMThreadException::CaptureException(){
    //std::unique_lock<std::mutex> guard(this->Lock);
#pragma omp critical
    {
        this->Ptr = std::current_exception();
    }
}

bool gMThreadException::exceptionOccurred(){
    return this->Ptr!=nullptr;
}



template<typename T>
gMSparse<T> geco::methods::toSparse(const gMDense<T> & dmat){
    arma::SpMat<T> ret(armadillo::convert<gMDense,arma::Mat,T>(dmat));
    return armadillo::convert<arma::SpMat,gMSparse,T>(ret);
}
        
template<typename T>
gMDense<T> geco::methods::toDense(const gMSparse<T> & smat){
    arma::Mat<T> ret(armadillo::convert<gMSparse,arma::SpMat,T>(move(smat)));
    return armadillo::convert<arma::Mat,gMDense,T>(ret);
}

template<typename T>
T geco::methods::mse(const gMDense<T> & M1, const gMDense<T> & M2){
    if(M1.i_nrows!=M2.i_nrows || M1.i_ncols!=M2.i_ncols){
        throw gMethodException("MSE: dimensions of M1 and M2 must be the same");
    }
    arma::Mat<T> m1(armadillo::convert<gMDense,arma::Mat,T>(move(M1)));
    arma::Mat<T> m2(armadillo::convert<gMDense,arma::Mat,T>(move(M2)));
    return sqrt( arma::mean(arma::mean(square(m1-m2)))/ arma::mean(arma::mean(square(m1))) );
}


//FILL DENSE MATRICES
template<typename T>
class geco::methods::fill_method<gMDense,T>{
    gMDense<T> operator () (indexType nrows,indexType ncols,T value) const{
        arma::Mat<T> ret(nrows,ncols);
        ret.fill(value);
        return armadillo::convert<arma::Mat,gMDense,T>(ret);
    }
};

template<typename T>
std::vector<indexType>  geco::methods::find(const gMDense<T> & M1, T threshold, bool greater){
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        try{
            cuda::autoCudaHandles cudaHandles;
            cuda::cudaDense<T> cM1=cuda::convert<gMDense,cuda::cudaDense,T>(M1,false,cudaHandles);
            std::vector<cuda::cudaIndexType> ret=cuda::find(cM1,threshold,greater);
            return std::vector<indexType>(ret.begin(),ret.end());
        }catch(gMethodException & e){
            throw e;
        }catch(exception & e){
            throw e;
        }catch(...){
            throw gMethodException("unknown ecxeption");
        }
    }
#endif
    arma::Mat<T> aM1=armadillo::convert<gMDense,arma::Mat>(move(M1));
    arma::Mat<arma::uword> ret;
    if(greater){
        ret=arma::find(aM1>threshold);
    }else{
        ret=arma::find(aM1<threshold);
    }
    return std::vector<indexType>(ret.begin(),ret.end());
}

template<typename T>
gMDense<T> geco::methods::select_columns(gMDense<T> M1, const std::vector<indexType> columns){
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        try{
            cuda::autoCudaHandles cudaHandles;
            cuda::cudaDense<T> cM1=cuda::convert<gMDense,cuda::cudaDense,T>(move(M1),false,cudaHandles);
            std::vector<cudaIndexType> cColumns(columns.begin(),columns.end());
            cuda::select_columns_inplace<T>(cM1,cColumns);
            return cuda::convert<cudaDense,gMDense,T>(cM1,false,cudaHandles);
        }catch(gMethodException & e){
            throw e;
        }catch(exception & e){
            throw e;
        }catch(...){
            throw gMethodException("unknown ecxeption");
        }
    }
#endif
    arma::Mat<T> aM1=armadillo::convert<gMDense,arma::Mat>(move(M1));
    arma::Mat<T> ret=aM1.cols(arma::uvec(columns));
    return armadillo::convert<arma::Mat,gMDense,T>(move(ret));
}

template<typename T>
std::vector< T > geco::methods::randperm(T m, T n){
    arma::uvec tmp=arma::randperm(m,n);
    vector<T> scols(n);
    std::copy(tmp.begin(),tmp.end(),scols.begin());
    return scols;
}

//RANDOM DENSE MATRICES
template<typename T>
class geco::methods::initRandom_method<gMDense,T>{
    gMDense<T> operator () (indexType nrows,indexType ncols,T sparseness) const{
        arma::Mat<T> ret=armadillo::initRandom_dense(nrows,ncols,sparseness);
        return armadillo::convert<arma::Mat,gMDense,T>(ret);
    }
};

//RANDOM SPARSE MATRICES
template<typename T>
class geco::methods::initRandom_method<gMSparse,T>{
    gMSparse<T> operator () (indexType nrows,indexType ncols,T sparseness) const{
        arma::SpMat<T> ret=armadillo::initRandom_sparse(nrows,ncols,sparseness);
        return armadillo::convert<arma::SpMat,gMSparse,T>(ret);
    }
};




//DENSE DENSE specialization
template<typename T>
class geco::methods::multiply_method<gMDense,gMDense,T>{
public:
    gMDense<T> operator () (const gMDense<T> & M1, const gMDense<T> & M2, bool transposeA,bool transposeB) const{
#ifdef GECO_HAS_CUDA
        if(useCuda()){
            try{
                autoCudaHandles cudaHandles;
                cudaDense<T> C((transposeA)?(M1.i_ncols):(M1.i_nrows),(transposeB)?(M2.i_nrows):(M2.i_ncols));
                if(&M1==&M2){
                    cudaDense<T> TM=cuda::convert<gMDense,cudaDense,T>(M1,false,cudaHandles);
                    multiply_cuda<T>(TM,TM,C,transposeA,transposeB,1,0,cudaHandles);
                }else{
                    multiply_cuda<T>(cuda::convert<gMDense,cudaDense,T>(M1,false,cudaHandles),cuda::convert<gMDense,cudaDense,T>(M2,false,cudaHandles),C,transposeA,transposeB,1,0,cudaHandles);
                }
                return std::move(cuda::convert<cudaDense,gMDense,T>(std::move(C),false,cudaHandles));
            }catch(gMethodException & e){
                throw e;
            }catch(exception & e){
                throw e;
            }catch(...){
                throw gMethodException("unknown ecxeption");
            }
        }
#endif
        arma::Mat<T> R=armadillo::multiply_arma(armadillo::convert<gMDense,arma::Mat,T>(M1),armadillo::convert<gMDense,arma::Mat,T>(M2),transposeA,transposeB);
        return armadillo::convert<arma::Mat,gMDense,T>(R);
    }
};

//SPARSE DENSE specialization
template<typename T>
class geco::methods::multiply_method<gMSparse,gMDense,T>{
public:
    gMDense<T> operator ()(const gMSparse<T> & M1, const gMDense<T> & M2, bool transposeA, bool transposeB) const{
#ifdef GECO_HAS_CUDA
        if(useCuda()){
            try{
                autoCudaHandles cudaHandles;

                cudaDense<T> C((transposeA)?(M1.i_ncols):(M1.i_nrows),(transposeB)?(M2.i_nrows):(M2.i_ncols));
                multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(move(M1),false,cudaHandles),cuda::convert<gMDense,cudaDense,T>(move(M2),false,cudaHandles),C,transposeA,transposeB,1,0,cudaHandles);
                return std::move(cuda::convert<cudaDense,gMDense,T>(std::move(C),false,cudaHandles));
            }catch(gMethodException & e){
                throw e;
            }catch(exception & e){
                throw e;
            }catch(...){
                throw gMethodException("unknown ecxeption");
            }
        }
#endif
        return std::move(armadillo::convert<arma::Mat,gMDense,T>(armadillo::multiply_arma(armadillo::convert<gMSparse,arma::SpMat,T>(M1),armadillo::convert<gMDense,arma::Mat,T>(M2),transposeA,transposeB)));
    }
};

//DENSE SPARSE specialization
template<typename T>
class geco::methods::multiply_method<gMDense,gMSparse,T>{
public:
    gMSparse<T> operator ()(const gMDense<T> & M1, const gMSparse<T> & M2, bool transposeA, bool transposeB) const{
#ifdef GECO_HAS_CUDA
        if(useCuda()){
//             try{
//                 autoCublasHandle cublasHandle;
//                 autoCusparseHandle cusparseHandle;
//                 cudaDense<T> C((transposeA)?(M1.i_ncols):(M1.i_nrows),(transposeB)?(M2.i_nrows):(M2.i_ncols));
//                 multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(move(M1),false,cublasHandle,cusparseHandle),cuda::convert<gMDense,cudaDense,T>(move(M2),false,cublasHandle,cusparseHandle),C,transposeA,transposeB,1,0,cusparseHandle);
//                 return std::move(cuda::convert<cudaDense,gMDense,T>(std::move(C),false,cublasHandle,cusparseHandle));
//             }catch(gException & e){
//             }            
              cout << "warning dense matrix multiplication not implemented in cuda using armadillo instead" << endl;
              return std::move(armadillo::convert<arma::SpMat,gMSparse,T>(armadillo::multiply_arma(armadillo::convert<gMDense,arma::Mat,T>(M1),armadillo::convert<gMSparse,arma::SpMat,T>(M2),transposeA,transposeB)));
        }
#endif
        return move(armadillo::convert<arma::SpMat,gMSparse,T>(armadillo::multiply_arma(armadillo::convert<gMDense,arma::Mat,T>(M1),armadillo::convert<gMSparse,arma::SpMat,T>(M2),transposeA,transposeB)));
    }
};

//SPARSE SPARSE specialization
template<typename T>
class geco::methods::multiply_method<gMSparse,gMSparse,T>{
public:
    gMSparse<T> operator ()(const gMSparse<T> & M1, const gMSparse<T> & M2, bool transposeA, bool transposeB) const{
#ifdef GECO_HAS_CUDA
        if(useCuda()){
            try{
                autoCudaHandles cudaHandles;
                indexType Cr=(transposeA)?(M1.i_ncols):(M1.i_nrows);
                indexType Cc=(transposeB)?(M2.i_nrows):(M2.i_ncols);
                gMSparse<T> ret;
                if(&M1==&M2){
                    if(!transposeA && !transposeB){
                        cudaSparse<T> C(Cc,Cr,0);
                        cudaSparse<T> TM=move(cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles));
                        multiply_cuda<T>(TM,TM,C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),true,cudaHandles));
                    }else if(transposeA && transposeB){
                        cudaSparse<T> C(Cr,Cc,0);
                        cudaSparse<T> TM=move(cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles));
                        multiply_cuda<T>(TM,TM,C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),false,cudaHandles));
                    }else if(transposeA && !transposeB){
                        cudaSparse<T> C(Cc,Cr,0);
                        multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles),cuda::convert<gMSparse,cudaSparse,T>(M1,false,cudaHandles),C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),true,cudaHandles));
                    }else if(!transposeA && transposeB){
                        cudaSparse<T> C(Cc,Cr,0);
                        multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(M1,false,cudaHandles),cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles),C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),true,cudaHandles));
                    }
                }else{
                    if(!transposeA && !transposeB){
                        cudaSparse<T> C(Cc,Cr,0);
                        multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(M2,true,cudaHandles),cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles),C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),true,cudaHandles));
                    }else if(transposeA && transposeB){
                        cudaSparse<T> C(Cr,Cc,0);
                        multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles),cuda::convert<gMSparse,cudaSparse,T>(M2,true,cudaHandles),C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),false,cudaHandles));
                    }else if(transposeA && !transposeB){
                        cudaSparse<T> C(Cr,Cc,0);
                        multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(M1,true,cudaHandles),cuda::convert<gMSparse,cudaSparse,T>(M2,false,cudaHandles),C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),false,cudaHandles));
                    }else if(!transposeA && transposeB){
                        cudaSparse<T> C(Cr,Cc,0);
                        multiply_cuda<T>(cuda::convert<gMSparse,cudaSparse,T>(M1,false,cudaHandles),cuda::convert<gMSparse,cudaSparse,T>(M2,true,cudaHandles),C,1,cudaHandles);
                        ret=move(cuda::convert<cudaSparse,gMSparse,T>(move(C),false,cudaHandles));
                    }
                }
                return move(ret);
            }catch(gMethodException & e){
                throw e;
            }catch(exception & e){
                throw e;
            }catch(...){
                throw gMethodException("unknown ecxeption");
            }
        }
#endif
        return move(armadillo::convert<arma::SpMat,gMSparse,T>(armadillo::multiply_arma(armadillo::convert<gMSparse,arma::SpMat,T>(M1),armadillo::convert<gMSparse,arma::SpMat,T>(M2),transposeA,transposeB)));
    }
};


template gMDense<float> geco::methods::toDense<float>(const gMSparse<float> & smat);
template gMDense<double> geco::methods::toDense<double>(const gMSparse<double> & smat);

template gMSparse<float> geco::methods::toSparse<float>(const gMDense<float> & dmat);
template gMSparse<double> geco::methods::toSparse<double>(const gMDense<double> & dmat);

template  float geco::methods::mse<float>(const gMDense<float> & M1, const gMDense<float> & M2);
template double geco::methods::mse<double>(const gMDense<double> & M1, const gMDense<double> & M2);

template class geco::methods::fill_method<gMDense,float>;
template class geco::methods::fill_method<gMDense,double>;

template vector<indexType> geco::methods::find(const gMDense<float> & M1, float threshold, bool greater);
template vector<indexType> geco::methods::find(const gMDense<double> & M1, double threshold, bool greater);

template gMDense<float> geco::methods::select_columns(gMDense<float> M1, const std::vector<indexType> columns);
template gMDense<double> geco::methods::select_columns(gMDense<double> M1, const std::vector<indexType> columns);


template std::vector<indexType> geco::methods::randperm(indexType n, indexType m);

#ifdef GECO_HAS_CUDA
//cuda::cudaIndexType is the same than indexType. Should be better to keep intantiations separated but this require cuda
template std::vector<cuda::cudaIndexType> geco::methods::randperm(cuda::cudaIndexType n, cuda::cudaIndexType m);
#endif


template class geco::methods::initRandom_method<gMDense,float>;
template class geco::methods::initRandom_method<gMDense,double>;

template class geco::methods::initRandom_method<gMSparse,float>;
template class geco::methods::initRandom_method<gMSparse,double>;


template class geco::methods::multiply_method<gMDense,gMDense,float>;
template class geco::methods::multiply_method<gMDense,gMDense,double>;

template class geco::methods::multiply_method<gMSparse,gMDense,float>;
template class geco::methods::multiply_method<gMSparse,gMDense,double>;

template class geco::methods::multiply_method<gMDense,gMSparse,float>;
template class geco::methods::multiply_method<gMDense,gMSparse,double>;

template class geco::methods::multiply_method<gMSparse,gMSparse,float>;
template class geco::methods::multiply_method<gMSparse,gMSparse,double>;



