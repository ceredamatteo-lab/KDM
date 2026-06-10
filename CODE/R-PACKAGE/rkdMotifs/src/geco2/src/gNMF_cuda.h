#pragma once
#include "gNMF.h"
#include "gMethods_cuda.h"

namespace geco{
    namespace methods{
        namespace cuda{
            
            template<typename T>
            T nmf_change(const cudaDense<T> & W,const cudaDense<T> & oldW);
            
            template<typename T>
            cudaDense<T> & get_sigma(const cudaDense<T> & dsigma,cudaDense<T> & sigma);
            
            template<typename T>
            cudaDense<T> & set_dgma(const cudaDense<T> & sigma, cudaDense<T> & dsigma);

            template<typename T>
            cudaDense<T> & set_dgma2(const cudaDense<T> & sigma, cudaDense<T> & dsigma);
            
            template <typename T>
            void printout(const cudaDense<T> & M,const autoCudaHandles & cudaHandles);
            
            //ATTENTION V must be transposed (for efficiency)
            template<template<typename> class M,typename T>
            std::pair<T,indexType> PNMF(const M<T> & Vt,cudaDense<T> & W, bool ard, cudaDense<T> & sigma, T tolerance, indexType maxIter, bool normalise, autoCudaHandles & cudaHandles){
                gMethodsOut inner(2);
                T change=1;
                indexType iterations=0;
                std::vector<cudaIndexType> G;

                cudaDense<T> A(Vt.i_nrows,W.i_ncols);
                cudaDense<T> N(W.i_nrows,W.i_ncols);
                cudaDense<T> D(W.i_nrows,W.i_ncols);
                cudaDense<T> C(W.i_ncols,W.i_ncols);
                
                cudaDense<T> WtW(W.i_ncols,W.i_ncols);
                cudaDense<T> E(W.i_ncols,1);
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    multiply_cuda<T>(W,W,WtW,true,false,1,0,cudaHandles); //WtW=W'W
                    if(ard){
                        get_diag(WtW,sigma);

                        G=find(sigma,std::numeric_limits<T>::epsilon(),true);
                        if(G.size()<sigma.i_ncols){
                            select_columns_inplace(sigma,G);
                            select_columns_inplace(W,G);
                            A.updateNcols(W.i_ncols);
                            N.updateNcols(W.i_ncols);
                            D.updateNcols(W.i_ncols);
                            C.updateNrows(W.i_ncols);
                            C.updateNcols(W.i_ncols);
                            WtW.updateNrows(W.i_ncols);
                            WtW.updateNcols(W.i_ncols);
                            E.updateNrows(W.i_ncols);
                            multiply_cuda<T>(W,W,WtW,true,false,1,0,cudaHandles); //WtW=W'W MAYBE BETTER TO HAVE select_inplace(vecrows,veccols)....
                        }
                    }
                    
                    //n:
                    multiply_cuda<T>(Vt,W,A,false,false,1,0,cudaHandles); //A=V'W
                    multiply_cuda<T>(Vt,A,N,true,false,1,0,cudaHandles); //N=VA=VV'W
                    
                    //d:
                    multiply_cuda<T>(W,N,C,true,false,1,0,cudaHandles); //C= W'N = W'VV'W
                    multiply_cuda<T>(W,C,D,false,false,1,0,cudaHandles); //D= WC=WW'N = WW'VV'W
                    multiply_cuda<T>(N,WtW,D,false,false,1,1,cudaHandles); //D = D + N * WtW = WW'VV'W + VV'W * W'W
                    
                    if(ard){
                        multiply_cuda<T>(W,set_dgma(sigma,C),D,false,false,1,1,cudaHandles); //D = D + W * C where C<-diag(1/S)
                    }
                    elementwise_add_constant(D,std::numeric_limits<T>::epsilon(),D);
                    
                    // W <- W * (N / D)
                    elementwise_divide<T>(N,D,N); //N <- N / D
                    elementwise_multiply<T>(W,N,N); // N <- W * N <- W * N/D
                    
                    //W <- W / sqrt(max(eigenv(W'W)
                    syrk<T>(N,C,true,true,1.0,0.0,cudaHandles);
                    syevd(C,E,true,cudaHandles);
                    elementwise_multiply_constant(N,1/::sqrt(cuda::elementwise_max(E)),N);
                    change=nmf_change(N,W);
                    W=N;
                    
                    iterations++;
                    checkInterrupt();                    
                    
                    if(iterations % 10==0){
                        if(ard){
                            inner.End([=]{ std::cout << "\tard rank: " << sigma.i_ncols << "\tchange: " << change;});
                        }else{
                            inner.End([=]{ std::cout << "\tchange: " << change;});
                        }
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                
                if(ard){
                    inner.End([=]{ std::cout << "\tard rank: " << sigma.i_ncols << "\tchange: " << change;});
                }else{
                    inner.End([=]{ std::cout << "\tchange: " << change;});
                }

                if(iterations>=maxIter){
                    warning(std::string("PNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    normalise_inplace(W);
                }

                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);
                return std::make_pair(0,iterations);
            }
            
            template<template<typename> class M,typename T>
            std::pair<T,indexType> OPNMF(const M<T> & Vt,cudaDense<T> & W, cudaIndexType p, cudaDense<T> & Q, T tolerance, cudaIndexType maxIter, bool normalise, autoCudaHandles & cudaHandles){
                gMethodsOut inner(2);
                T change=1;
                indexType iterations=0;
                std::vector<cudaIndexType> G;

                cudaDense<T> A(p,W.i_ncols);
                cudaDense<T> N(W.i_nrows,W.i_ncols);
                cudaDense<T> D(W.i_nrows,W.i_ncols);
                cudaDense<T> C(W.i_ncols,W.i_ncols);
                
                cudaDense<T> WtW(W.i_ncols,W.i_ncols);
                cudaDense<T> E(W.i_ncols,1);
                
                //prime_time = std::chrono::high_resolution_clock::now();
                //second_time = std::chrono::high_resolution_clock::now();
                
                //printout(Q,cudaHandles);
                //printout(Q,cudaHandles);
                //second_time = std::chrono::high_resolution_clock::now();
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    multiply_cuda<T>(W,W,WtW,true,false,1,0,cudaHandles); //WtW=W'W
                    
                    M<T> pvt=select_rows(Vt,randperm(Vt.i_nrows,p));
                    //Q:
                    multiply_cuda<T>(pvt,W,A,false,false,1,0,cudaHandles); //A=pV'W
                    multiply_cuda<T>(pvt,A,Q,true,false,1,1,cudaHandles); //Q = Q + pVA = Q + pVpV'W
                    
                    //d:
                    multiply_cuda<T>(W,Q,C,true,false,1,0,cudaHandles); //C= W'Q = W'pVpV'W
                    multiply_cuda<T>(W,C,D,false,false,1,0,cudaHandles); //D= WC=WW'Q = WW'pVpV'W
                    multiply_cuda<T>(Q,WtW,D,false,false,1,1,cudaHandles); //D = D + Q * WtW = WW'pVpV'W + pVpV'W * W'W
                    elementwise_add_constant(D,std::numeric_limits<T>::epsilon(),D);
                    
                    // W <- W * (N / D)
                    elementwise_divide<T>(Q,D,N); //N <- N / D
                    elementwise_multiply<T>(W,N,N); // N <- W * N <- W * N/D
                    
                    //W <- W / sqrt(max(eigenv(W'W)
                    syrk<T>(N,C,true,true,1.0,0.0,cudaHandles);
                    syevd(C,E,true,cudaHandles);
                    elementwise_multiply_constant(N,1/::sqrt(cuda::elementwise_max(E)),N);
                    
                    change=nmf_change(N,W);
                    
                    W=N;
                    
                    iterations++;
                    checkInterrupt();  
                    
                    if(iterations % 10==0){
                        inner.End([=]{std::cout << "\tchange=" << change;});
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                inner.End([=]{std::cout << "\tchange=" << change;});
                
                if(iterations>=maxIter){
                    warning(std::string("OPNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    normalise_inplace(W);
                }

                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);

                return std::make_pair(0,iterations);
            }

            //ATTENTION V must be transposed (for efficiency)
            template<template<typename> class B,template<typename> class M,typename T>
            std::pair<T,indexType> PNMF2(const nmf::dataFeeder<B,T> & feeder,cudaDense<T> & W, bool ard, cudaDense<T> & sigma, T tolerance, indexType maxIter, bool normalise, autoCudaHandles & cudaHandles){
                gMethodsOut inner(2);
                T change=1;
                indexType iterations=0;
                std::vector<cudaIndexType> G;
                
                M<T> Vt=cuda::convert<B,M,T>(feeder(0),true,cudaHandles);

                cudaDense<T> A(Vt.i_nrows,W.i_ncols);
                cudaDense<T> N(W.i_nrows,W.i_ncols);
                cudaDense<T> D(W.i_nrows,W.i_ncols);
                cudaDense<T> C(W.i_ncols,W.i_ncols);
                
                cudaDense<T> WtW(W.i_ncols,W.i_ncols);
                cudaDense<T> E(W.i_ncols,1);
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    multiply_cuda<T>(W,W,WtW,true,false,1,0,cudaHandles); //WtW=W'W
                    if(ard){
                        get_diag(WtW,sigma);

                        G=find(sigma,std::numeric_limits<T>::epsilon(),true);
                        if(G.size()<sigma.i_ncols){
                            select_columns_inplace(sigma,G);
                            select_columns_inplace(W,G);
                            A.updateNcols(W.i_ncols);
                            N.updateNcols(W.i_ncols);
                            D.updateNcols(W.i_ncols);
                            C.updateNrows(W.i_ncols);
                            C.updateNcols(W.i_ncols);
                            WtW.updateNrows(W.i_ncols);
                            WtW.updateNcols(W.i_ncols);
                            E.updateNrows(W.i_ncols);
                            multiply_cuda<T>(W,W,WtW,true,false,1,0,cudaHandles); //WtW=W'W MAYBE BETTER TO HAVE select_inplace(vecrows,veccols)....
                        }
                    }
                    
                    //n:
                    multiply_cuda<T>(Vt,W,A,false,false,1,0,cudaHandles); //A=V'W
                    multiply_cuda<T>(Vt,A,N,true,false,1,0,cudaHandles); //N=VA=VV'W
                    
                    //d:
                    multiply_cuda<T>(W,N,C,true,false,1,0,cudaHandles); //C= W'N = W'VV'W
                    multiply_cuda<T>(W,C,D,false,false,1,0,cudaHandles); //D= WC=WW'N = WW'VV'W
                    multiply_cuda<T>(N,WtW,D,false,false,1,1,cudaHandles); //D = D + N * WtW = WW'VV'W + VV'W * W'W
                    
                    if(ard){
                        multiply_cuda<T>(W,set_dgma(sigma,C),D,false,false,1,1,cudaHandles); //D = D + W * C where C<-diag(1/S)
                    }
                    elementwise_add_constant(D,std::numeric_limits<T>::epsilon(),D);
                    
                    // W <- W * (N / D)
                    elementwise_divide<T>(N,D,N); //N <- N / D
                    elementwise_multiply<T>(W,N,N); // N <- W * N <- W * N/D
                    
                    //W <- W / sqrt(max(eigenv(W'W)
                    syrk<T>(N,C,true,true,1.0,0.0,cudaHandles);
                    syevd(C,E,true,cudaHandles);
                    elementwise_multiply_constant(N,1/::sqrt(cuda::elementwise_max(E)),N);
                    
                    change=nmf_change(N,W);
                    
                    W=N;
                    
                    iterations++;
                    checkInterrupt();                    
                    
                    if(iterations % 10==0){
                        if(ard){
                            inner.End([=]{ std::cout << "\tard rank: " << sigma.i_ncols << "\tchange: " << change;});
                        }else{
                            inner.End([=]{ std::cout << "\tchange: " << change;});
                        }
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                
                if(ard){
                    inner.End([=]{ std::cout << "\tard rank: " << sigma.i_ncols << "\tchange: " << change;});
                }else{
                    inner.End([=]{ std::cout << "\tchange: " << change;});
                }

                if(iterations>=maxIter){
                    warning(std::string("PNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    normalise_inplace(W);
                }

                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);
                return std::make_pair(0,iterations);
            }
            
            
            template<template<typename> class B,template<typename> class M,typename T>
            std::pair<T,indexType> OPNMF2(const nmf::dataFeeder<B,T> & feeder,cudaDense<T> & W, cudaIndexType p, cudaDense<T> & Q, T tolerance, cudaIndexType maxIter, bool normalise, autoCudaHandles & cudaHandles){
                gMethodsOut inner(2);
                T change=1;
                indexType iterations=0;
                std::vector<cudaIndexType> G;

                cudaDense<T> A(p,W.i_ncols);
                cudaDense<T> N(W.i_nrows,W.i_ncols);
                cudaDense<T> D(W.i_nrows,W.i_ncols);
                cudaDense<T> C(W.i_ncols,W.i_ncols);
                
                cudaDense<T> WtW(W.i_ncols,W.i_ncols);
                cudaDense<T> E(W.i_ncols,1);
                
                //prime_time = std::chrono::high_resolution_clock::now();
                //second_time = std::chrono::high_resolution_clock::now();
                
                //printout(Q,cudaHandles);
                //printout(Q,cudaHandles);
                //second_time = std::chrono::high_resolution_clock::now();
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    multiply_cuda<T>(W,W,WtW,true,false,1,0,cudaHandles); //WtW=W'W
                    
                    M<T> pvt=cuda::convert<B,M,T>(feeder(p),true,cudaHandles);
                    //Q:
                    multiply_cuda<T>(pvt,W,A,false,false,1,0,cudaHandles); //A=pV'W
                    multiply_cuda<T>(pvt,A,Q,true,false,1,1,cudaHandles); //Q = Q + pVA = Q + pVpV'W
                    
                    //d:
                    multiply_cuda<T>(W,Q,C,true,false,1,0,cudaHandles); //C= W'Q = W'pVpV'W
                    multiply_cuda<T>(W,C,D,false,false,1,0,cudaHandles); //D= WC=WW'Q = WW'pVpV'W
                    multiply_cuda<T>(Q,WtW,D,false,false,1,1,cudaHandles); //D = D + Q * WtW = WW'pVpV'W + pVpV'W * W'W
                    elementwise_add_constant(D,std::numeric_limits<T>::epsilon(),D);
                    
                    // W <- W * (N / D)
                    elementwise_divide<T>(Q,D,N); //N <- N / D
                    elementwise_multiply<T>(W,N,N); // N <- W * N <- W * N/D
                    
                    //W <- W / sqrt(max(eigenv(W'W)
                    syrk<T>(N,C,true,true,1.0,0.0,cudaHandles);
                    syevd(C,E,true,cudaHandles);
                    elementwise_multiply_constant(N,1/::sqrt(cuda::elementwise_max(E)),N);
                    
                    change=nmf_change(N,W);
                    
                    W=N;
                    
                    iterations++;
                    checkInterrupt();  
                    
                    if(iterations % 10==0){
                        inner.End([=]{std::cout << "\tchange=" << change;});
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                inner.End([=]{std::cout << "\tchange=" << change;});
                
                if(iterations>=maxIter){
                    warning(std::string("OPNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    normalise_inplace(W);
                }

                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);

                return std::make_pair(0,iterations);
            }
            
            
        }
    }
}


