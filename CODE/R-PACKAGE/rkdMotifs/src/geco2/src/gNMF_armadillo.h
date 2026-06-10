#pragma once
#include "gMethods_armadillo.h"
#include <chrono>

template<template<typename> class  M,typename T>
class svd_impl{
public:
    bool operator ()(arma::Mat<T> & U, arma::Col<T> & s, arma::Mat<T> & V, const M<T> & X, geco::methods::indexType k){
//         if(std::is_same<M<T>,arma::SpMat<T>>::value){
//             return arma::svds(U,s,V,X,k);
//         }else{
//             return arma::svd_econ(U,s,V,X);
//         }
        return false;
    }
};

template<typename T>
class svd_impl<arma::SpMat,T>{
public:
    bool operator ()(arma::Mat<T> & U, arma::Col<T> & s, arma::Mat<T> & V, const arma::SpMat<T> & X,geco::methods::indexType k){
        return arma::svds(U,s,V,X,k);
    }
};

template<typename T>
class svd_impl<arma::Mat,T>{
public:
    bool operator ()(arma::Mat<T> & U, arma::Col<T> & s, arma::Mat<T> & V, const arma::Mat<T> & X,geco::methods::indexType k){
        return arma::svd_econ(U,s,V,X);
    }
};

template<template<typename> class  M,typename T>
bool svd(arma::Mat<T> & U, arma::Col<T> & s, arma::Mat<T> & V, const M<T> & X, geco::methods::indexType k){
    svd_impl<M,T> method;
    return method(U,s,V,X,k);
}

template<typename T>
T calcChange(const arma::Mat<T> & W, const arma::Mat<T> & oldW){
    return pow(arma::norm(oldW-W,"fro")/arma::norm(oldW,"fro"),2);
}

namespace geco{
    namespace methods{
        namespace armadillo{
            
            template<template<typename> class M,typename T>
            void NNDSVD(M<T> A,indexType rank, unsigned short flag, bool initW, bool initH, arma::Mat<T> & W, arma::Mat<T> & H){
                // This function implements the NNDSVD algorithm described in [1] for
                // initialization of Nonnegative Matrix Factorization Algorithms.
                //
                // [W,H] = nndsvd(A,k,flag);
                //
                // INPUT
                // ------------
                //
                // A    : the input nonnegative m x n matrix A
                // k    : the rank of the computed factors W,H
                // flag : indicates the variant of the NNDSVD Algorithm
                //        flag = 0 --> NNDSVD
                //        flag = 1 --> NNDSVDa
                //        flag = 2 --> NNDSVDar
                //
                // OUTPUT
                // -------------
                //   
                // W   : nonnegative m x k matrix
                // H   : nonnegative k x n matrix
                //
                // 
                // References:
                // 
                // [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
                //     start for nonnegative matrix factorization, Pattern Recognition,
                //     Elsevier
                //
                // This code is kindly provided by the authors for research porpuses.
                // - Efstratios Gallopoulos (stratis@ceid.upatras.gr)
                // - Christos Boutsidis (boutsc@cs.rpi.edu)     
                //
                // For any problems or questions please send an email to boutsc@cs.rpi.edu 
                //--------------------------------------------------------------------------                
                
                auto prime_time = std::chrono::high_resolution_clock::now();
                auto second_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> dur;
                
                //%size of the input matrix
                //[m,n] = size(A);
                arma::uword m=A.n_rows;
                arma::uword n=A.n_cols;
                
                //%the matrices of the factorization
                //W = zeros(m,k);
                //H = zeros(k,n);
                if(initW){
                    W = arma::zeros<arma::Mat<T>>(m,rank);
                }
                
                if(initH){
                    H = arma::zeros<arma::Mat<T>>(rank,n);
                }

                //% 1st SVD --> partial SVD rank-k to the input matrix A. 
                //[U,S,V] = svds(A,k);

                bool success=false;
                arma::Mat<T> U;
                arma::Col<T> s;
                arma::Mat<T> V;

                success=svd(U,s,V,A,rank);
                
                if(!success) throw gMethodException("NNDSVD: svd failed");
                
                //%choose the first singular triplet to be nonnegative
                //W(:,1)     =  sqrt(S(1,1)) * abs(U(:,1) );         
                //H(1,:)     =  sqrt(S(1,1)) * abs(V(:,1)'); 
                if(initW){
                    W.col(0) = sqrt(s(0)) * abs(U.col(0));
                }
                
                if(initH){
                    H.row(0) = sqrt(s(0)) * abs(trans(V.col(0)));
                }

                //% 2nd SVD for the other factors (see table 1 in our paper)
                //for i=2:k
                //    uu = U(:,i); vv = V(:,i);
                //    uup = pos(uu); uun = neg(uu) ;
                //    vvp = pos(vv); vvn = neg(vv);
                //    n_uup = norm(uup);
                //    n_vvp = norm(vvp) ;
                //    n_uun = norm(uun) ;
                //    n_vvn = norm(vvn) ;
                //    termp = n_uup*n_vvp; termn = n_uun*n_vvn;
                //    if (termp >= termn)
                //        W(:,i) = sqrt(S(i,i)*termp)*uup/n_uup; 
                //        H(i,:) = sqrt(S(i,i)*termp)*vvp'/n_vvp;
                //    else
                //        W(:,i) = sqrt(S(i,i)*termn)*uun/n_uun; 
                //        H(i,:) = sqrt(S(i,i)*termn)*vvn'/n_vvn;
                //    end
                //end
                //%------------------------------------------------------------
                for(arma::uword i=1;i<rank;i++){
                    arma::Col<T> uu = U.col(i); 
                    arma::Col<T> uup = (uu>=0) % uu;
                    arma::Col<T> uun = (uu<0) % (-uu);
                    T n_uup = norm(uup);
                    T n_uun = norm(uun);                    
                    
                    arma::Col<T> vv = V.col(i);
                    arma::Col<T> vvp = (vv>=0) % vv;
                    arma::Col<T> vvn = (vv<0) % (-vv);
                    T n_vvp = norm(vvp);
                    T n_vvn = norm(vvn);
                    
                    T termp = n_uup*n_vvp; 
                    T termn = n_uun*n_vvn;
                    
                    if (termp >= termn){
                        if(initW){
                            W.col(i) = sqrt(s(i)*termp)*uup/n_uup; 
                        }
                        if(initH){
                            H.row(i) = sqrt(s(i)*termp)*trans(vvp)/n_vvp;
                        }
                    }else{
                        if(initW){
                            W.col(i) = sqrt(s(i)*termn)*uun/n_uun; 
                        }
                        
                        if(initH){
                            H.row(i) = sqrt(s(i)*termn)*trans(vvn)/n_vvn;
                        }
                    }
                }
                
                //%actually these numbers are zeros
                //W((W<0.0000000001))=0;
                //H((H<0.0000000001))=0;
                if(initW){
                    W.elem(arma::find(W<0.0000000001)).fill(0);
                }
                if(initH){
                    H.elem(arma::find(H<0.0000000001)).fill(0);
                }
                


                //% NNDSVDa: fill in the zero elements with the average 
                //if flag==1
                //ind1      =  find(W==0) ;
                //ind2      =  find(H==0) ;
                //average   =  mean(A(:)) ; 
                //W( ind1 ) =  average    ; 
                //H( ind2 ) =  average    ;
                //
                //% NNDSVDar: fill in the zero elements with random values in the space [0:average/100]
                //elseif flag==2
                //ind1      =  find(W==0) ;
                //ind2      =  find(H==0) ;
                //n1        =  numel(ind1);
                //n2        =  numel(ind2);
                //
                //average   =  mean(A(:))       ;
                //W( ind1 ) =  (average*rand(n1,1)./100)  ; 
                //H( ind2 ) =  (average*rand(n2,1)./100)  ;   
                //end                
                if(flag!=0){
                    T average = mean(mean(A));
                    if(initW){
                        arma::uvec ind=arma::find(W==0);
                        if(flag==1){
                            W.elem(ind) += average;
                        }else if(flag==2){
                            W(ind) +=  (average*arma::randu<arma::Col<T>>(ind.n_elem,1)/100); 
                        }
                    }
                    if(initH){
                        arma::uvec ind=arma::find(H==0);
                        if(flag==1){
                            H.elem(ind) += average;
                        }else if(flag==2){
                            H.elem(ind) +=  (average*arma::randu<arma::Col<T>>(ind.n_elem,1)/100);   
                        }
                    }
                }
            }

            template<template<typename> class M,typename T>
            std::pair<T,indexType> PNMF(const M<T> & V,arma::Mat<T> & W, bool ard, arma::Mat<T> & sigma,T tolerance,indexType maxIter,bool normalise){
                gMethodsOut inner(2);
                
                /*
                // check_step = 100;
                // 
                // W = Winit;
                // 
                // XX = X * X';
                // for iter=1:max_iter
                //     W_old = W;
                //     if mod(iter,check_step)==0
                //         fprintf('iter=% 5d ', iter);
                //     end
                //     sigma = sum(bsxfun(@times, W, W));
                //     
                //     W = W .* ((XX*W)*diag(sigma) ./ ( (W*(W'*XX*W) + XX*W*(W'*W))*diag(sigma) + W + eps));
                //     W = W ./ norm(W);
                //     
                //     diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
                //     if diffW<tol
                //         fprintf('converged after %d steps.\n', iter);
                //         break;
                //     end
                //     
                //     if mod(iter,check_step)==0
                //         fprintf('diff=%.10f, ', diffW);
                //         fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
                //         fprintf('\n');
                //         fprintf('max sigma: %f    min sigma: %f\n', max(sigma), min(sigma));
                //     end
                // end
                // 
                // wnorm = [];
                // rinit = size(W,2);
                // for i = 1:rinit
                //     wnorm = [wnorm norm(W(:,i))];
                // end
                */
                
                arma::Mat<T> tW=trans(W);
                arma::Mat<T> prevW=W;
                arma::Mat<T> n,d,WtW;
                arma::Col<T> eigval;
                
                T change=10000000;
                indexType iterations=0;
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    WtW =  tW * W ;
                    
                    if(ard){
                        sigma=arma::Row<T>(WtW.diag());
                        arma::uvec G=arma::find(sigma>arma::Datum<T>::eps);
                        if(G.n_elem<sigma.n_elem){
                            sigma=sigma.elem(G);
                            W=W.cols(G);
                            tW=tW.rows(G);
                            prevW=prevW.cols(G);
                            WtW=WtW(G,G);
                        }
                    }
                    
                    n = V * (trans(V) * W);

                    d = (( W * ( tW * n ) ) + ( n * WtW ));
                    if(ard){
                        d += W * arma::diagmat(1/sigma);
                    }
                    d += arma::Datum<T>::eps;
                    
                    W = W % (n/d);
                    
                    eigval=arma::eig_sym(arma::trans(W) * W);
                    W = W / sqrt(max(eigval));
                    //W = W / arma::norm(W);  
                    
                    change=calcChange(W,prevW);
                    tW=arma::trans(W);
                    prevW=W;
                    iterations++;
                    checkInterrupt();
                    if(iterations % 10==0){
                        if(ard){
                            inner.End([=]{ std::cout << "\tard rank: " << sigma.n_cols << "\tchange: " << change;});
                        }else{
                            inner.End([=]{ std::cout << "\tchange: " << change;});
                        }
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                
                if(ard){
                    inner.End([=]{ std::cout << "\tard rank: " << sigma.n_cols << "\tchange: " << change;});
                }else{
                    inner.End([=]{ std::cout << "\tchange: " << change;});
                }
                
                if(iterations>=maxIter){
                    warning(std::string("PNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    W=arma::normalise(W);
                }
                
                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);
                return std::make_pair(0,iterations);
            }
            
            
            template<template<typename> class M,typename T>
            std::pair<T,indexType> OPNMF(const M<T> & V,arma::Mat<T> & W, indexType p, arma::Mat<T> & Q,T tolerance,indexType maxIter,bool normalise){
                gMethodsOut inner(2);
                /*
                // check_step = 100;
                // 
                // W = Winit;
                // 
                // XX = X * X';
                // for iter=1:max_iter
                //     W_old = W;
                //     if mod(iter,check_step)==0
                //         fprintf('iter=% 5d ', iter);
                //     end
                //     sigma = sum(bsxfun(@times, W, W));
                //     
                //     W = W .* ((XX*W)*diag(sigma) ./ ( (W*(W'*XX*W) + XX*W*(W'*W))*diag(sigma) + W + eps));
                //     W = W ./ norm(W);
                //     
                //     diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
                //     if diffW<tol
                //         fprintf('converged after %d steps.\n', iter);
                //         break;
                //     end
                //     
                //     if mod(iter,check_step)==0
                //         fprintf('diff=%.10f, ', diffW);
                //         fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
                //         fprintf('\n');
                //         fprintf('max sigma: %f    min sigma: %f\n', max(sigma), min(sigma));
                //     end
                // end
                // 
                // wnorm = [];
                // rinit = size(W,2);
                // for i = 1:rinit
                //     wnorm = [wnorm norm(W(:,i))];
                // end
                */
                

                arma::Mat<T> tW=trans(W);
                arma::Mat<T> prevW=W;
                arma::Mat<T> d,WtW;
                arma::Col<T> eigval;
                
                T change=10000000;
                indexType iterations=0;
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    
                    WtW =  tW * W ;
                    M<T> pV=V.cols(arma::randperm(V.n_cols,p));
                    Q += pV * (trans(pV) * W);

                    d = (( W * ( tW * Q ) ) + ( Q * WtW ));
                    d += arma::Datum<T>::eps;
                    
                    W = W % (Q/d);
                    
                    eigval=arma::eig_sym(arma::trans(W) * W);
                    W = W / sqrt(max(eigval));
                    //W = W / arma::norm(W);  
                    
                    
                    change=calcChange(W,prevW);
                    tW=arma::trans(W);
                    prevW=W;
                    iterations++;
                    checkInterrupt();
                    if(iterations % 10==0){
                        inner.End([=]{ std::cout << "\tchange: " << change;});
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                inner.End([=]{ std::cout << "\tchange: " << change;});
                
                if(iterations>=maxIter){
                    warning(std::string("OPNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    W=arma::normalise(W);
                }

                
                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);
                return std::make_pair(0,iterations);
            }

            template<template<typename> class B,template<typename> class M, typename T>
            std::pair<T,indexType> PNMF2(const nmf::dataFeeder<B,T> & feeder,arma::Mat<T> & W, bool ard, arma::Mat<T> & sigma,T tolerance,indexType maxIter,bool normalise){
                gMethodsOut inner(2);
                
                /*
                // check_step = 100;
                // 
                // W = Winit;
                // 
                // XX = X * X';
                // for iter=1:max_iter
                //     W_old = W;
                //     if mod(iter,check_step)==0
                //         fprintf('iter=% 5d ', iter);
                //     end
                //     sigma = sum(bsxfun(@times, W, W));
                //     
                //     W = W .* ((XX*W)*diag(sigma) ./ ( (W*(W'*XX*W) + XX*W*(W'*W))*diag(sigma) + W + eps));
                //     W = W ./ norm(W);
                //     
                //     diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
                //     if diffW<tol
                //         fprintf('converged after %d steps.\n', iter);
                //         break;
                //     end
                //     
                //     if mod(iter,check_step)==0
                //         fprintf('diff=%.10f, ', diffW);
                //         fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
                //         fprintf('\n');
                //         fprintf('max sigma: %f    min sigma: %f\n', max(sigma), min(sigma));
                //     end
                // end
                // 
                // wnorm = [];
                // rinit = size(W,2);
                // for i = 1:rinit
                //     wnorm = [wnorm norm(W(:,i))];
                // end
                */
                
                M<T> V = armadillo::convert<B,M,T>(feeder(0));
                arma::Mat<T> tW=trans(W);
                arma::Mat<T> prevW=W;
                arma::Mat<T> n,d,WtW;
                arma::Col<T> eigval;
                
                T change=10000000;
                indexType iterations=0;
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    WtW =  tW * W ;
                    
                    if(ard){
                        sigma=arma::Row<T>(WtW.diag());
                        arma::uvec G=arma::find(sigma>arma::Datum<T>::eps);
                        if(G.n_elem<sigma.n_elem){
                            sigma=sigma.elem(G);
                            W=W.cols(G);
                            tW=tW.rows(G);
                            prevW=prevW.cols(G);
                            WtW=WtW(G,G);
                        }
                    }
                    
                    n = V * (trans(V) * W);

                    d = (( W * ( tW * n ) ) + ( n * WtW ));
                    if(ard){
                        d += W * arma::diagmat(1/sigma);
                    }
                    d += arma::Datum<T>::eps;
                    
                    W = W % (n/d);
                    
                    eigval=arma::eig_sym(arma::trans(W) * W);
                    W = W / sqrt(max(eigval));
                    //W = W / arma::norm(W);  
                    
                    change=calcChange(W,prevW);
                    tW=arma::trans(W);
                    prevW=W;
                    iterations++;
                    checkInterrupt();
                    if(iterations % 10==0){
                        if(ard){
                            inner.End([=]{ std::cout << "\tard rank: " << sigma.n_cols << "\tchange: " << change;});
                        }else{
                            inner.End([=]{ std::cout << "\tchange: " << change;});
                        }
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                
                if(ard){
                    inner.End([=]{ std::cout << "\tard rank: " << sigma.n_cols << "\tchange: " << change;});
                }else{
                    inner.End([=]{ std::cout << "\tchange: " << change;});
                }
                
                if(iterations>=maxIter){
                    warning(std::string("PNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    W=arma::normalise(W);
                }
                
                //return std::make_pair(norm(V-(W * (tW * V)),"fro"),iterations);
                return std::make_pair(0,iterations);
            }


            template<template<typename> class B,template<typename> class M, typename T>
            std::pair<T,indexType> OPNMF2(const nmf::dataFeeder<B,T> & feeder,arma::Mat<T> & W, indexType p, arma::Mat<T> & Q,T tolerance,indexType maxIter,bool normalise){
                gMethodsOut inner(2);
                /*
                // check_step = 100;
                // 
                // W = Winit;
                // 
                // XX = X * X';
                // for iter=1:max_iter
                //     W_old = W;
                //     if mod(iter,check_step)==0
                //         fprintf('iter=% 5d ', iter);
                //     end
                //     sigma = sum(bsxfun(@times, W, W));
                //     
                //     W = W .* ((XX*W)*diag(sigma) ./ ( (W*(W'*XX*W) + XX*W*(W'*W))*diag(sigma) + W + eps));
                //     W = W ./ norm(W);
                //     
                //     diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
                //     if diffW<tol
                //         fprintf('converged after %d steps.\n', iter);
                //         break;
                //     end
                //     
                //     if mod(iter,check_step)==0
                //         fprintf('diff=%.10f, ', diffW);
                //         fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
                //         fprintf('\n');
                //         fprintf('max sigma: %f    min sigma: %f\n', max(sigma), min(sigma));
                //     end
                // end
                // 
                // wnorm = [];
                // rinit = size(W,2);
                // for i = 1:rinit
                //     wnorm = [wnorm norm(W(:,i))];
                // end
                */
                

                arma::Mat<T> tW=trans(W);
                arma::Mat<T> prevW=W;
                arma::Mat<T> d,WtW;
                arma::Col<T> eigval;
                
                T change=10000000;
                indexType iterations=0;
                
                inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                while((change>=tolerance) && (iterations<maxIter)){
                    
                    WtW =  tW * W ;
                    //M<T> pV=V.cols(arma::randperm(V.n_cols,p));
                    M<T> pV=convert<B,M,T>(feeder(p));
                    Q += pV * (arma::trans(pV) * W);

                    d = (( W * ( tW * Q ) ) + ( Q * WtW ));
                    d += arma::Datum<T>::eps;
                    
                    W = W % (Q/d);
                    
                    eigval=arma::eig_sym(arma::trans(W) * W);
                    W = W / sqrt(max(eigval));
                    //W = W / arma::norm(W);  
                    
                    
                    change=calcChange(W,prevW);
                    tW=arma::trans(W);
                    prevW=W;
                    iterations++;
                    checkInterrupt();
                    if(iterations % 10==0){
                        inner.End([=]{ std::cout << "\tchange: " << change;});
                        inner.Start([=]{std::cout << "\titerations [" << (iterations+1) << "-" << (iterations+10) << "]...";});
                    }
                }
                inner.End([=]{ std::cout << "\tchange: " << change;});
                
                if(iterations>=maxIter){
                    warning(std::string("OPNMF: Warning! maxIter reached before convergence"));
                }
                
                if(normalise){
                    W=arma::normalise(W);
                }

                return std::make_pair(0,iterations);
            }

            
            
        }
    }
}



