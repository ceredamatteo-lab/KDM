#pragma once
#include "gMethods.h"
#include <cstdarg>
#include <boost/any.hpp>
//#include <any>

namespace geco{
    namespace methods{
        namespace nmf{
            
//BEGIN  NMF methods            
            template< template<typename> class D,typename T>
            class gNMF{
            public:
                gMDense<T> i_W;
                gMDense<T> i_H;
                
                bool i_init=false;
                
                T i_error=0;
                indexType i_iterations=0;

                gNMF(){

                }
                
                gNMF(gNMF & other){
                    i_W=other.i_W;
                    i_H=other.i_H;
                    i_init=other.i_init;
                    i_error=other.i_error;
                    i_iterations=other.i_iterations;
                }
                
                gNMF & operator = (gNMF other){
                    swap(*this,other);
                    return *this;
                }
                
                gNMF(gNMF && other){
                    swap(*this,other);
                }
                
//                 gNMF & operator = (gNMF && other){
//                     std::cout << "move assigned" << std::endl;
//                     swap(*this,other);
//                 }
                
                friend void swap(gNMF<D,T> & a,gNMF<D,T> & b){
                    using std::swap;
                    swap(a.i_W,b.i_W);
                    swap(a.i_H,b.i_H);
                    swap(a.i_init,b.i_init);
                    swap(a.i_error,b.i_error);
                    swap(a.i_iterations,b.i_iterations);
                }
                
                void operator ()(const D<T> & V, T tolerance, indexType maxIterations){
                    
                }
            };

            template<template<typename> class D, typename T>
            class PNMF:public gNMF<D,T>{
            public:
                gMDense<T> i_sigma;                
                bool i_ard;
                bool i_normW;
                PNMF(bool ard,bool normalizeW):i_ard(ard),i_normW(normalizeW){
                }
                void operator ()(const D<T> & V, T tolerance, indexType maxIterations);
            };
            
            template<template<typename> class D, typename T>
            class OPNMF:public gNMF<D,T>{
            public:
                indexType i_nnew;
                bool i_normW;
                gMDense<T> i_Q;
                OPNMF(indexType nnew,bool normalizeW):i_nnew(nnew),i_normW(normalizeW){
                }
                void operator ()(const D<T> & V, T tolerance, indexType maxIterations);
            };
//END  NMF methods
            
//BEGIN NMF initialization methods
            
            template<template<typename> class D,typename T>
            class baseNMFInit{
            public:
                bool i_initW;
                bool i_initH;
                baseNMFInit(bool initW,bool initH):i_initW(initW),i_initH(initH){
                }
            };
            
            template<template<typename> class D,typename T>
            class customNMFInit:public baseNMFInit<D,T>{
                gMDense<T> i_W0;
                gMDense<T> i_H0;
            public:
                customNMFInit(gMDense<T>  W0, gMDense<T>  H0):baseNMFInit<D,T>(W0.i_nrows>0,H0.i_nrows>0){
                    i_W0=W0;
                    i_H0=H0;
                }
                bool operator ()(const D<T> & V, gNMF<D,T> & M, indexType rank) const{
                    gMethodsOut msg(1);
                    msg.Start([=]{std::cout << "Custom NMF initialization...";});
                    if(baseNMFInit<D,T>::i_initW){
                        if(i_W0.i_ncols!=rank|| i_W0.i_nrows!=V.i_nrows){
                            throw gMethodException("customNMFInit: W mismatched dimensions");
                        }
                        M.i_W=i_W0;
                    }
                    if(baseNMFInit<D,T>::i_initH){
                        if(i_H0.i_ncols!=V.i_ncols || i_H0.i_nrows!=rank){
                            throw gMethodException("customNMFInit: H mismatched dimensions");
                        }
                        M.i_H=i_H0;
                    }
                    msg.End([=]{});
                    return true;
                }
            };
            
            template<template<typename> class D, typename T>
            class randomNMFInit:public baseNMFInit<D,T>{
            public:
                randomNMFInit(bool initW,bool initH):baseNMFInit<D,T>(initW,initH){
                }
                bool operator ()(const D<T> & V, gNMF<D,T> & M, indexType rank) const{
                    gMethodsOut msg(1);
                    msg.Start([=]{std::cout << "Random NMF initialization...";});
                    if(baseNMFInit<D,T>::i_initW){
                        M.i_W=initRandom<gMDense,T>(V.i_nrows,rank,1);
                    }
                    if(baseNMFInit<D,T>::i_initH){
                        M.i_H=initRandom<gMDense,T>(rank,V.i_ncols,1);
                    }
                    msg.End([=]{});
                    return true;
                }
            };
           
            template<template<typename> class D, typename T>
            class nndsvdNMFInit:public baseNMFInit<D,T>{
            public:
                unsigned short i_flag;
                
                nndsvdNMFInit(bool initW=true,bool initH=true,unsigned short flag=0):baseNMFInit<D,T>(initW,initH),i_flag(flag){
                }
                
                bool operator ()(const D<T> & V, gNMF<D,T> & M, indexType rank) const;
            };
//END NMF initialization methods

            template<template<typename> class D, template<template<typename> class,typename> class I,template<template<typename> class, typename> class M,typename T>
            M<D,T> NMF(const D<T> & V, const I<D,T> & initMethod, M<D,T> method, indexType rank,T tolerance, indexType maxIterations){
                method.i_init=initMethod(V,method,rank);
                method(V,tolerance,maxIterations);
                return method;
            }
            
            
//BEGIN  NMF2 methods                        
            template<template<typename> class D,typename T>
            class dataFeeder{
            public:
                boost::any i_data;
                D<T> (*i_datafunction)(const boost::any & data,indexType);
                dataFeeder(const boost::any & data,D<T> (*datafunction)(const boost::any & ,indexType)):i_data(data),i_datafunction(datafunction){
                }
                
                D<T> operator()(indexType n) const{
                    return i_datafunction(i_data,n);
                }
            };


            template< template<typename> class D,typename T>
            class gNMF2{
            public:
                gMDense<T> i_W;
                gMDense<T> i_H;
                
                bool i_init=false;
                
                T i_error=0;
                indexType i_iterations=0;
                

                gNMF2(){

                }
                
                gNMF2(gNMF2 & other){
                    i_W=other.i_W;
                    i_H=other.i_H;
                    i_init=other.i_init;
                    i_error=other.i_error;
                    i_iterations=other.i_iterations;
                }
                
                gNMF2 & operator = (gNMF2 other){
                    swap(*this,other);
                    return *this;
                }
                
                gNMF2(gNMF2 && other){
                    swap(*this,other);
                }
                
//                 gNMF & operator = (gNMF && other){
//                     std::cout << "move assigned" << std::endl;
//                     swap(*this,other);
//                 }
                
                friend void swap(gNMF2<D,T> & a,gNMF2<D,T> & b){
                    using std::swap;
                    swap(a.i_W,b.i_W);
                    swap(a.i_H,b.i_H);
                    swap(a.i_init,b.i_init);
                    swap(a.i_error,b.i_error);
                    swap(a.i_iterations,b.i_iterations);
                }
                
                void operator ()(const dataFeeder<D,T> & feeder, T tolerance, indexType maxIterations){
                    
                }
            };
            
            template<template<typename> class D, typename T>
            class PNMF2:public gNMF2<D,T>{
            public:
                gMDense<T> i_sigma;                
                bool i_ard;
                bool i_normW;
                PNMF2(bool ard,bool normalizeW):i_ard(ard),i_normW(normalizeW){
                }
                void operator ()(const dataFeeder<D,T> & feeder, T tolerance, indexType maxIterations);
            };
            
            template<template<typename> class D, typename T>
            class OPNMF2:public gNMF2<D,T>{
            public:
                indexType i_nnew;
                bool i_normW;
                gMDense<T> i_Q;
                OPNMF2(indexType nnew,bool normalizeW):i_nnew(nnew),i_normW(normalizeW){
                }
                void operator ()(const dataFeeder<D,T> & feeder, T tolerance, indexType maxIterations);
            };
//END  NMF2 methods                                    
            
//BEGIN NMF2 initialization methods
            
            template<template<typename> class D,typename T>
            class baseNMFInit2{
            public:
                bool i_initW;
                bool i_initH;
                baseNMFInit2(bool initW,bool initH):i_initW(initW),i_initH(initH){
                }
            };
            
            template<template<typename> class D,typename T>
            class customNMFInit2:public baseNMFInit2<D,T>{
                gMDense<T> i_W0;
                gMDense<T> i_H0;
            public:
                customNMFInit2(gMDense<T>  W0, gMDense<T>  H0):baseNMFInit2<D,T>(W0.i_nrows>0,H0.i_nrows>0){
                    i_W0=W0;
                    i_H0=H0;
                }
                bool operator ()(const dataFeeder<D,T> & feeder, gNMF2<D,T> & M, indexType rank) const{
                    gMethodsOut msg(1);
                    msg.Start([=]{std::cout << "Custom NMF initialization...";});
                    if(baseNMFInit2<D,T>::i_initW){
                        if(i_W0.i_ncols!=rank|| i_W0.i_nrows!=feeder(1).i_nrows){
                            throw gMethodException("customNMFInit: W mismatched dimensions");
                        }
                        M.i_W=i_W0;
                    }
                    if(baseNMFInit2<D,T>::i_initH){
                        if(i_H0.i_ncols!=feeder(1).i_ncols || i_H0.i_nrows!=rank){
                            throw gMethodException("customNMFInit: H mismatched dimensions");
                        }
                        M.i_H=i_H0;
                    }
                    msg.End([=]{});
                    return true;
                }
            };
            
            template<template<typename> class D, typename T>
            class randomNMFInit2:public baseNMFInit2<D,T>{
            public:
                randomNMFInit2(bool initW,bool initH):baseNMFInit2<D,T>(initW,initH){
                }
                bool operator ()(const dataFeeder<D,T> & feeder, gNMF2<D,T> & M, indexType rank) const{
                    gMethodsOut msg(1);
                    msg.Start([=]{std::cout << "Random NMF initialization...";});
                    if(baseNMFInit2<D,T>::i_initW){
                        M.i_W=initRandom<gMDense,T>(feeder(1).i_nrows,rank,1);
                    }
                    if(baseNMFInit2<D,T>::i_initH){
                        M.i_H=initRandom<gMDense,T>(rank,feeder(1).i_ncols,1);
                    }
                    msg.End([=]{});
                    return true;
                }
            };
           
            template<template<typename> class D, typename T>
            class nndsvdNMFInit2:public baseNMFInit2<D,T>{
            public:
                unsigned short i_flag;
                
                nndsvdNMFInit2(bool initW=true,bool initH=true,unsigned short flag=0):baseNMFInit2<D,T>(initW,initH),i_flag(flag){
                }
                
                bool operator ()(const dataFeeder<D,T> & feeder, gNMF2<D,T> & M, indexType rank) const;
            };
//END NMF2 initialization methods

            template<template<typename> class D, template<template<typename> class,typename> class I,template<template<typename> class, typename> class M,typename T>
            M<D,T> NMF2(const dataFeeder<D,T> & feeder, const I<D,T> & initMethod, M<D,T> method, indexType rank,T tolerance, indexType maxIterations){
                method.i_init=initMethod(feeder,method,rank);
                method(feeder,tolerance,maxIterations);
                return method;
            }
            
        }
    }
}
