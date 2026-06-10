#pragma once
#include "gMethods.h"
#include "gPWM.h"

#include <list>
#include <map>

namespace geco{
    namespace methods{
        namespace kmers{

            template<typename T>
            class gLmerCounter{
           protected:
                bool i_doubleStrand;
                indexType i_l;            
                indexType i_nlmers;
                indexType i_nDSlmers;
                std::vector<indexType> i_lmerDSLmer;
                std::vector<indexType> i_DSlmers;
                
                //std::map<size_t,std::list<unsigned char>> seqCode(const std::string  & seq,indexType sID,bool strict) const;

            public:
                gLmerCounter(indexType l,bool doubleStrand=false);
                gMDense<T> countDense(const std::vector<std::string> & sequences,bool strict) const;
                gMSparse<T> countSparse(const std::vector<std::string> & sequences,bool strict) const;
                
                gMDense<T> countWinDense(const std::string & sequence,indexType win,bool strict) const;
                gMSparse<T> countWinSparse(const std::string & sequence,indexType win,bool strict) const;
                
                gMDense<T> countPWM(const geco::kmers::gPWMSet & PWMs,indexType padSize) const;
                gMDense<T> getPWM(const gMDense<T> & W,indexType pwmlength,const std::vector<double> & pwm0,double tol,indexType maxiter) const;
                gMDense<T> countPWM2(const geco::kmers::gPWMSet & PWMs) const;
                std::pair< std::vector<std::string>,std::vector<std::string> >getStrings() const;
            };
            
            template<typename T>
            class gGappedKmerCounter:public gLmerCounter<T>{
            protected:
                indexType i_k;
                indexType i_nmasks;
                indexType i_nkmers;
                indexType i_ngkmers;
                indexType i_nmlmers;
                indexType i_nDSgkmers;
                
                std::vector<indexType> i_masks;
                std::vector<indexType> i_DSgkmers;

                std::vector<indexType> i_mlmerGkmer;
                std::vector<indexType> i_gkmerDSGkmer;
                
            public:
                gGappedKmerCounter(indexType l,indexType k,bool doubleStrand=false);
                gMDense<T> countDense(const std::vector<std::string> & sequences,bool strict) const;
                gMSparse<T> countSparse(const std::vector<std::string> & sequences,bool strict) const;
                
                gMDense<T> countWinDense(const std::string & sequence,indexType win,bool strict) const;
                gMSparse<T> countWinSparse(const std::string & sequence,indexType win,bool strict) const;
                
                gMDense<T> countPWM(const geco::kmers::gPWMSet & PWMs,indexType padSize) const;
                gMDense<T> getPWM(const gMDense<T> & W,indexType pwmlength,const std::vector<double> & pwm0,double tol,indexType maxiter) const;
                std::pair< std::vector<std::string>,std::vector<std::string> >getStrings() const;
            };
        
#ifdef GECO_HAS_CUDA
            template<typename T> class gCudaCounter;
#endif
        
            template<typename T>
            class gLmerEstimator:public gGappedKmerCounter<T>{
#ifdef GECO_HAS_CUDA
                friend class gCudaCounter<T>;
#endif
            protected:
                std::vector<T> i_mmcoeffs;
                std::map<indexType,indexType> i_gkmersNDSgkmers; //tells how many DSgkmers a gkmers contribute to
                std::vector<indexType> i_mm;
#ifdef GECO_HAS_CUDA
                void * i_cudaObject;
#endif
            
            
            public:
                gLmerEstimator(indexType l,indexType k,bool doubleStrand=false);
#ifdef GECO_HAS_CUDA
                ~gLmerEstimator();
#endif
            
                gMDense<T> countDense(const std::vector<std::string> & sequences,bool strict) const;
                gMSparse<T> countSparse(const std::vector<std::string> & sequences,bool strict) const;
                
                gMDense<T> countDense(const gMDense<T> & gkmers) const;
                
                gMSparse<T> calculateG();
            };
            
            template<template<typename> class D,typename T>
            D<T> & toHellinger(D<T> & counts,bool use_cuda=false);
            
            template<typename T>
            std::vector<T> getPWM(const gMDense<T> & W,indexType l,indexType pwmlength,const std::vector<double> & pwm0,double tol=1e-4,indexType maxiter=10000);
        }
    }
}
