#include "gNMF.h"
#ifdef GECO_HAS_CUDA
#include "gNMF_cuda.h"
#endif

#include "gNMF_armadillo.h"
using namespace std;
using namespace geco;
using namespace geco::methods;
using namespace geco::methods::nmf;


//BEGIN initNNDSVD
template<template<typename> class D, typename T>
bool nndsvdNMFInit<D,T>::operator ()(const D<T> & V, gNMF<D,T> & M, indexType rank) const{
    gMethodsOut msg(1);
    if(useCuda()){
        msg.Start([=]{cout << "NNDSVD NMF initialization (using cuda) is not implemented yet, using armadillo instead";});
        msg.End([=]{});
    }
    arma::Mat<T> W,H;
    if(std::is_same<D<T>,gMDense<T>>::value){
        msg.Start([=]{cout << "NNDSVD NMF dense initialization (using armadillo)...";});
        armadillo::NNDSVD(armadillo::convert<D,arma::Mat>(V),rank,i_flag,baseNMFInit<D,T>::i_initW,baseNMFInit<D,T>::i_initH,W,H);
    }else{
        msg.Start([=]{cout << "NNDSVD NMF sparse initialization (using armadillo)...";});
        armadillo::NNDSVD(armadillo::convert<D,arma::SpMat>(V),rank,i_flag,baseNMFInit<D,T>::i_initW,baseNMFInit<D,T>::i_initH,W,H);
    }
    M.i_W=armadillo::convert<arma::Mat,gMDense>(move(W));
    M.i_H=armadillo::convert<arma::Mat,gMDense>(move(H));
    msg.End([=]{});
    return true;
}
template class geco::methods::nmf::nndsvdNMFInit<gMDense,float>;
template class geco::methods::nmf::nndsvdNMFInit<gMDense,double>;
template class geco::methods::nmf::nndsvdNMFInit<gMSparse,float>;
template class geco::methods::nmf::nndsvdNMFInit<gMSparse,double>;
//END initNNDSVD


//BEGIN PNMF
template<template<typename> class D, typename T>
void PNMF<D,T>::operator ()(const D<T> & V, T tolerance, indexType maxIterations){
    gMethodsOut msgs(1);
    if(std::is_same<D<T>,gMDense<T>>::value){
        msgs.Start([=]{
            cout << "PNMF" << ((i_ard)?("ARD"):("")) << " Dense specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF<D,T>::i_W.i_ncols <<  ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });
    }else{
        msgs.Start([=]{
            cout << "PNMF" << ((i_ard)?("ARD"):("")) << " Sparse specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF<D,T>::i_W.i_ncols << ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });            
    }

#ifdef GECO_HAS_CUDA
    if(useCuda()){
        try{
            cuda::autoCudaHandles cudaHandles;
            cuda::cudaDense<T> aW=cuda::convert<gMDense,cuda::cudaDense>(gNMF<D,T>::i_W,false,cudaHandles);
            cuda::cudaDense<T> sigma(1,aW.i_ncols);
            std::pair<T,indexType> ret;

            if(std::is_same<D<T>,gMDense<T>>::value){
                ret=cuda::PNMF(cuda::convert<D,cuda::cudaDense>(move(V),true,cudaHandles),aW,i_ard,sigma,tolerance,maxIterations, i_normW,cudaHandles);
            }else{
                ret=cuda::PNMF( cuda::convert<D,cuda::cudaSparse>(move(V),true,cudaHandles),aW,i_ard,sigma,tolerance,maxIterations, i_normW,cudaHandles);
            }
            gNMF<D,T>::i_W=cuda::convert<cuda::cudaDense,gMDense>(move(aW),false,cudaHandles);
            gNMF<D,T>::i_error=ret.first;
            gNMF<D,T>::i_iterations=ret.second;

            i_sigma=cuda::convert<cuda::cudaDense,gMDense>(move(sigma),false,cudaHandles);
            msgs.End([=]{
                cout << "\tTotal iterations: " << ret.second << "\t Total time: ";
            });
        }catch(gMethodException &e){
            cout << "PNMF failed with error: \"" << e.what() << "\"." << endl;
            cout << "Try setCudaUsage(false)" << endl;
        }
        return;
    }
#endif
    arma::Mat<T> aW=armadillo::convert<gMDense,arma::Mat,T>(move(gNMF<D,T>::i_W));
    arma::Mat<T> sigma(gNMF<D,T>::i_W.i_ncols,1);
    std::pair<T,indexType> ret;
    if(std::is_same<D<T>,gMDense<T>>::value){
        //arma::Mat<T> tV=arma::trans(armadillo::convert<D,arma::Mat>(move(V)));
        ret=armadillo::PNMF(armadillo::convert<D,arma::Mat>(V),aW,i_ard,sigma,tolerance,maxIterations,i_normW);
    }else{
        //arma::SpMat<T> tV=arma::trans(armadillo::convert<D,arma::SpMat>(move(V)));
        ret=armadillo::PNMF(armadillo::convert<D,arma::SpMat>(V),aW,i_ard,sigma,tolerance,maxIterations,i_normW);
    }
    gNMF<D,T>::i_W=armadillo::convert<arma::Mat,gMDense,T>(move(aW));
    gNMF<D,T>::i_error=ret.first;
    gNMF<D,T>::i_iterations=ret.second;
    i_sigma=armadillo::convert<arma::Mat,gMDense,T>(move(sigma));
    msgs.End([=]{
        cout << "\tTotal iterations: " << ret.second << "\t Total time: ";
    });    
    return;
}
template class geco::methods::nmf::PNMF<gMDense,float>;
template class geco::methods::nmf::PNMF<gMDense,double>;
template class geco::methods::nmf::PNMF<gMSparse,float>;
template class geco::methods::nmf::PNMF<gMSparse,double>;
//END PNMF

//BEGIN OPNMF
template<template<typename> class D, typename T>
void OPNMF<D,T>::operator ()(const D<T> & V, T tolerance, indexType maxIterations){
    gMethodsOut msgs(1);
    if(std::is_same<D<T>,gMDense<T>>::value){
        msgs.Start([=]{
            cout << "OPNMF Dense specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF<D,T>::i_W.i_ncols <<  ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });
    }else{
        msgs.Start([=]{
            cout << "OPNMF Sparse specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF<D,T>::i_W.i_ncols << ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });            
    }
        
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        try{
            cuda::autoCudaHandles cudaHandles;
            cuda::cudaDense<T> aW=cuda::convert<gMDense,cuda::cudaDense>(gNMF<D,T>::i_W,false,cudaHandles);
            std::pair<T,indexType> ret;

            cuda::cudaDense<T> Q(aW.i_nrows,aW.i_ncols,0);
            if(std::is_same<D<T>,gMDense<T>>::value){
                ret=cuda::OPNMF(cuda::convert<D,cuda::cudaDense>(move(V),true,cudaHandles),aW,i_nnew,Q,tolerance,maxIterations, i_normW,cudaHandles);
            }else{
                ret=cuda::OPNMF( cuda::convert<D,cuda::cudaSparse>(move(V),true,cudaHandles),aW,i_nnew,Q,tolerance,maxIterations, i_normW,cudaHandles);
            }
            gNMF<D,T>::i_W=cuda::convert<cuda::cudaDense,gMDense>(move(aW),false,cudaHandles);
            gNMF<D,T>::i_error=ret.first;
            gNMF<D,T>::i_iterations=ret.second;
            i_Q=cuda::convert<cuda::cudaDense,gMDense,T>(move(Q),false,cudaHandles);
            msgs.End([=]{cout << "\tTotal iterations: " << ret.second << "\t Total time: ";});    
        }catch(gMethodException &e){
            cout << "PNMF failed with error: \"" << e.what() << "\"." << endl;
            cout << "Try setting use_cuda to false" << endl;
        }
        return;
    }
#endif
    arma::Mat<T> aW=armadillo::convert<gMDense,arma::Mat,T>(move(gNMF<D,T>::i_W));
    arma::Mat<T> Q(aW.n_rows,aW.n_cols,arma::fill::zeros);
    std::pair<T,indexType> ret;
    if(std::is_same<D<T>,gMDense<T>>::value){
        ret=armadillo::OPNMF(armadillo::convert<D,arma::Mat>(V),aW,i_nnew,Q,tolerance,maxIterations,i_normW);
    }else{
        ret=armadillo::OPNMF(armadillo::convert<D,arma::SpMat>(V),aW,i_nnew,Q,tolerance,maxIterations,i_normW);
    }
    gNMF<D,T>::i_W=armadillo::convert<arma::Mat,gMDense,T>(move(aW));
    gNMF<D,T>::i_error=ret.first;
    gNMF<D,T>::i_iterations=ret.second;
    i_Q=armadillo::convert<arma::Mat,gMDense,T>(move(Q));
    msgs.End([=]{cout << "\tTotal iterations: " << ret.second << "\t Total time: ";});
    return;
}
template class geco::methods::nmf::OPNMF<gMDense,float>;
template class geco::methods::nmf::OPNMF<gMDense,double>;
template class geco::methods::nmf::OPNMF<gMSparse,float>;
template class geco::methods::nmf::OPNMF<gMSparse,double>;
//END PNMF


//BEGIN initNNDSVD2
template<template<typename> class D, typename T>
bool nndsvdNMFInit2<D,T>::operator ()(const dataFeeder<D,T> & feeder, gNMF2<D,T> & M, indexType rank) const{
    gMethodsOut msg(1);
    if(useCuda()){
        msg.Start([=]{cout << "NNDSVD NMF initialization (using cuda) is not implemented yet, using armadillo instead";});
        msg.End([=]{});
    }
    arma::Mat<T> W,H;
    if(std::is_same<D<T>,gMDense<T>>::value){
        msg.Start([=]{cout << "NNDSVD NMF dense initialization (using armadillo)...";});
        armadillo::NNDSVD(armadillo::convert<D,arma::Mat>(feeder(0)),rank,i_flag,baseNMFInit2<D,T>::i_initW,baseNMFInit2<D,T>::i_initH,W,H);
    }else{
        msg.Start([=]{cout << "NNDSVD NMF sparse initialization (using armadillo)...";});
        armadillo::NNDSVD(armadillo::convert<D,arma::SpMat>(feeder(0)),rank,i_flag,baseNMFInit2<D,T>::i_initW,baseNMFInit2<D,T>::i_initH,W,H);
    }
    M.i_W=armadillo::convert<arma::Mat,gMDense>(move(W));
    M.i_H=armadillo::convert<arma::Mat,gMDense>(move(H));
    msg.End([=]{});
    return true;
}
template class geco::methods::nmf::nndsvdNMFInit2<gMDense,float>;
template class geco::methods::nmf::nndsvdNMFInit2<gMDense,double>;
template class geco::methods::nmf::nndsvdNMFInit2<gMSparse,float>;
template class geco::methods::nmf::nndsvdNMFInit2<gMSparse,double>;
//END initNNDSVD2


//BEGIN PNMF2
template<template<typename> class D, typename T>
void PNMF2<D,T>::operator ()(const dataFeeder<D,T> & feeder, T tolerance, indexType maxIterations){
    gMethodsOut msgs(1);
    if(std::is_same<D<T>,gMDense<T>>::value){
        msgs.Start([=]{
            cout << "PNMF" << ((i_ard)?("ARD"):("")) << " Dense specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF2<D,T>::i_W.i_ncols <<  ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });
    }else{
        msgs.Start([=]{
            cout << "PNMF" << ((i_ard)?("ARD"):("")) << " Sparse specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF2<D,T>::i_W.i_ncols << ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });            
    }

#ifdef GECO_HAS_CUDA

    if(useCuda()){
        try{
            cuda::autoCudaHandles cudaHandles;
            cuda::cudaDense<T> aW=cuda::convert<gMDense,cuda::cudaDense>(gNMF2<D,T>::i_W,false,cudaHandles);
            cuda::cudaDense<T> sigma(1,aW.i_ncols);
            std::pair<T,indexType> ret;

            if(std::is_same<D<T>,gMDense<T>>::value){
                ret=cuda::PNMF2<D,cuda::cudaDense,T>(feeder,aW,i_ard,sigma,tolerance,maxIterations, i_normW,cudaHandles);
            }else{
                ret=cuda::PNMF2<D,cuda::cudaSparse,T>(feeder,aW,i_ard,sigma,tolerance,maxIterations, i_normW,cudaHandles);
            }
            gNMF2<D,T>::i_W=cuda::convert<cuda::cudaDense,gMDense>(move(aW),false,cudaHandles);
            gNMF2<D,T>::i_error=ret.first;
            gNMF2<D,T>::i_iterations=ret.second;

            i_sigma=cuda::convert<cuda::cudaDense,gMDense>(move(sigma),false,cudaHandles);
            msgs.End([=]{
                cout << "\tTotal iterations: " << ret.second << "\t Total time: ";
            });
        }catch(gMethodException &e){
            cout << "PNMF failed with error: \"" << e.what() << "\"." << endl;
            cout << "Try setCudaUsage(false)" << endl;
        }
        return;
    }
#endif
    arma::Mat<T> aW=armadillo::convert<gMDense,arma::Mat,T>(move(gNMF2<D,T>::i_W));
    arma::Mat<T> sigma(gNMF2<D,T>::i_W.i_ncols,1);
    std::pair<T,indexType> ret;
    if(std::is_same<D<T>,gMDense<T>>::value){
        //arma::Mat<T> tV=arma::trans(armadillo::convert<D,arma::Mat>(move(V)));
        ret=armadillo::PNMF2<D,arma::Mat,T>(feeder,aW,i_ard,sigma,tolerance,maxIterations,i_normW);
    }else{
        //arma::SpMat<T> tV=arma::trans(armadillo::convert<D,arma::SpMat>(move(V)));
        ret=armadillo::PNMF2<D,arma::SpMat,T>(feeder,aW,i_ard,sigma,tolerance,maxIterations,i_normW);
    }
    gNMF2<D,T>::i_W=armadillo::convert<arma::Mat,gMDense,T>(move(aW));
    gNMF2<D,T>::i_error=ret.first;
    gNMF2<D,T>::i_iterations=ret.second;
    i_sigma=armadillo::convert<arma::Mat,gMDense,T>(move(sigma));
    msgs.End([=]{
        cout << "\tTotal iterations: " << ret.second << "\t Total time: ";
    });    
    return;
}
template class geco::methods::nmf::PNMF2<gMDense,float>;
template class geco::methods::nmf::PNMF2<gMDense,double>;
template class geco::methods::nmf::PNMF2<gMSparse,float>;
template class geco::methods::nmf::PNMF2<gMSparse,double>;
//END PNMF2


//BEGIN OPNMF2
template<template<typename> class D, typename T>
void OPNMF2<D,T>::operator ()(const dataFeeder<D,T> & feeder, T tolerance, indexType maxIterations){
    gMethodsOut msgs(1);
    if(std::is_same<D<T>,gMDense<T>>::value){
        msgs.Start([=]{
            cout << "OPNMF Dense specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF2<D,T>::i_W.i_ncols <<  ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });
    }else{
        msgs.Start([=]{
            cout << "OPNMF Sparse specialization with " << ( (typeid(T)==typeid(float) )?("float"):("double") ) << " values and rank=" << gNMF2<D,T>::i_W.i_ncols << ((useCuda())?(" (using cuda)..."):(" (using armadillo)...")) << endl;
        });            
    }
        
#ifdef GECO_HAS_CUDA
    if(useCuda()){
        try{
            cuda::autoCudaHandles cudaHandles;
            cuda::cudaDense<T> aW=cuda::convert<gMDense,cuda::cudaDense>(gNMF2<D,T>::i_W,false,cudaHandles);
            std::pair<T,indexType> ret;

            cuda::cudaDense<T> Q(aW.i_nrows,aW.i_ncols,0);
            if(std::is_same<D<T>,gMDense<T>>::value){
                ret=cuda::OPNMF2<D,cuda::cudaDense,T>(feeder,aW,i_nnew,Q,tolerance,maxIterations, i_normW,cudaHandles);
            }else{
                ret=cuda::OPNMF2<D,cuda::cudaSparse,T>(feeder,aW,i_nnew,Q,tolerance,maxIterations, i_normW,cudaHandles);
            }
            gNMF2<D,T>::i_W=cuda::convert<cuda::cudaDense,gMDense>(move(aW),false,cudaHandles);
            gNMF2<D,T>::i_error=ret.first;
            gNMF2<D,T>::i_iterations=ret.second;
            i_Q=cuda::convert<cuda::cudaDense,gMDense,T>(move(Q),false,cudaHandles);
            msgs.End([=]{cout << "\tTotal iterations: " << ret.second << "\t Total time: ";});    
        }catch(gMethodException &e){
            cout << "PNMF failed with error: \"" << e.what() << "\"." << endl;
            cout << "Try setting use_cuda to false" << endl;
        }
        return;
    }
#endif
    arma::Mat<T> aW=armadillo::convert<gMDense,arma::Mat,T>(move(gNMF2<D,T>::i_W));
    arma::Mat<T> Q(aW.n_rows,aW.n_cols,arma::fill::zeros);
    std::pair<T,indexType> ret;
    if(std::is_same<D<T>,gMDense<T>>::value){
        ret=armadillo::OPNMF2<D,arma::Mat,T>(feeder,aW,i_nnew,Q,tolerance,maxIterations,i_normW);
    }else{
        ret=armadillo::OPNMF2<D,arma::SpMat,T>(feeder,aW,i_nnew,Q,tolerance,maxIterations,i_normW);
    }
    gNMF2<D,T>::i_W=armadillo::convert<arma::Mat,gMDense,T>(move(aW));
    gNMF2<D,T>::i_error=ret.first;
    gNMF2<D,T>::i_iterations=ret.second;
    i_Q=armadillo::convert<arma::Mat,gMDense,T>(move(Q));
    msgs.End([=]{cout << "\tTotal iterations: " << ret.second << "\t Total time: ";});
    return;
}
template class geco::methods::nmf::OPNMF2<gMDense,float>;
template class geco::methods::nmf::OPNMF2<gMDense,double>;
template class geco::methods::nmf::OPNMF2<gMSparse,float>;
template class geco::methods::nmf::OPNMF2<gMSparse,double>;
//END PNMF2

