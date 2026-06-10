#pragma once
#include "gMethods.h"

#if defined(GMETHODS_USE_32BITS_INDEX)
#define ARMA_32BIT_WORD
#else
#define ARMA_64BIT_WORD
#endif
#define ARMA_DONT_USE_WRAPPER
#define ARMA_EXTRA_DEBUG
#define ARMA_NO_DEBUG
//#define ARMA_EXTRA_DEBUG
#define ARMA_USE_OPENMP
#define ARMA_USE_CXX11
#include "armadillo/armadillo"

namespace geco{
    namespace methods{
        namespace armadillo{
            
//BEGIN armadillo conversion methods
            template<template<typename> class A, template<typename> class B, typename T>
            class convert_method{
            public:
                B<T> operator()(A<T>  && toConvert) const;
            };
            
            template<template<typename> class A, template<typename> class B, typename T>
            class convert_method_const_ref{
            public:
                const B<T> operator()(const A<T> & toConvert) const;
            };
            
            
            template<typename T>
            class convert_method<arma::Mat,gMDense,T>{
            public:
                gMDense<T> operator()(arma::Mat<T> && mat){
                    gMDense<T> ret;
                    if(mat.n_alloc==0 && mat.mem_state==0){
                        ret.i_values=gMMemory<T>(mat.mem,mat.n_elem);
                    }else{
                        ret.i_values.i_ptr=mat.mem;
                    }
                    ret.i_values.i_nelem=mat.n_elem;
                    ret.i_values.i_owns=true;
                    ret.i_nrows=mat.n_rows;
                    ret.i_ncols=mat.n_cols;
                    arma::access::rw(mat.mem)=nullptr;
                    arma::access::rw(mat.mem_state)=0;
                    arma::access::rw(mat.n_alloc)=0;
                    arma::access::rw(mat.n_cols)=0;
                    arma::access::rw(mat.n_rows)=0;
                    arma::access::rw(mat.n_elem)=0;
                    return std::move(ret);
                }
            };
            
            template<typename T>
            class convert_method<gMDense,arma::Mat,T>{
            public:
                arma::Mat<T> operator()(gMDense<T> && mat) const{
                    //std::cout << "move convert" << std::endl;
                    if(mat.i_values.i_nelem<arma::arma_config::mat_prealloc){
                        arma::Mat<T> ret((T*)mat.i_values,mat.i_nrows,mat.i_ncols,true);
                        return std::move(ret);
                    }else{
                        arma::Mat<T> ret;
                        arma::access::rw(ret.mem)=mat.i_values.i_ptr;
                        arma::access::rw(ret.mem_state)=0;
                        arma::access::rw(ret.n_alloc)=mat.i_values.i_nelem;
                        arma::access::rw(ret.n_cols)=mat.i_ncols;
                        arma::access::rw(ret.n_elem)=mat.i_values.i_nelem;
                        arma::access::rw(ret.n_rows)=mat.i_nrows;
                        mat.i_values.i_ptr=nullptr;
                        mat.i_values.i_nelem=0;
                        mat.i_values.i_owns=false;
                        return std::move(ret);
                    }
                }
            };
            
            template<typename T>
            class convert_method<arma::SpMat,gMSparse,T>{
            public:
                gMSparse<T> operator()(arma::SpMat<T>  && mat) const{
//                     cout << "to Sparse..";cout.flush();
                    gMSparse<T> ret;
                    ret.i_cscValues.release((T *) ret.i_cscValues);
                    ret.i_cscValues.i_ptr=mat.values;
                    ret.i_cscValues.i_nelem=mat.n_nonzero + 1;
                    ret.i_cscValues.i_owns=true;
                    
                    ret.i_cscColPtr.release((indexType *) ret.i_cscColPtr);
                    ret.i_cscColPtr.i_ptr=(indexType *) mat.col_ptrs;
                    ret.i_cscColPtr.i_nelem=mat.n_cols+2;
                    ret.i_cscColPtr.i_owns=true;
                    
                    ret.i_cscRowIndex.release((indexType *) ret.i_cscRowIndex);
                    ret.i_cscRowIndex.i_ptr=(indexType *) mat.row_indices;
                    ret.i_cscRowIndex.i_nelem=mat.n_nonzero + 1;
                    ret.i_cscRowIndex.i_owns=true;

                    ret.i_nrows=mat.n_rows;
                    ret.i_ncols=mat.n_cols;
                    ret.i_nnz=mat.n_nonzero;
                    
                    arma::access::rw(mat.values)=nullptr;
                    arma::access::rw(mat.col_ptrs)=nullptr;
                    arma::access::rw(mat.row_indices)=nullptr;
                    arma::access::rw(mat.n_nonzero)=0;
                    arma::access::rw(mat.n_cols)=0;
                    arma::access::rw(mat.n_rows)=0;
                    arma::access::rw(mat.n_elem)=0;
                    
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
            };
            
            template<typename T>
            class convert_method<gMSparse,arma::SpMat,T>{
            public:
                arma::SpMat<T> operator()(gMSparse<T> && mat) const{
//                     cout << "to SpMat..";cout.flush();
                    arma::SpMat<T> ret;
                    arma::memory::release(ret.values);
                    arma::memory::release(ret.col_ptrs);
                    arma::memory::release(ret.row_indices);
                    arma::access::rw(ret.values)=mat.i_cscValues.i_ptr;
                    arma::access::rw(ret.col_ptrs)=(arma::uword *) mat.i_cscColPtr.i_ptr;
                    arma::access::rw(ret.row_indices)=(arma::uword *) mat.i_cscRowIndex.i_ptr;
                    arma::access::rw(ret.n_nonzero)=mat.i_nnz;
                    arma::access::rw(ret.n_cols)=mat.i_ncols;
                    arma::access::rw(ret.n_rows)=mat.i_nrows;
                    arma::access::rw(ret.n_elem)=mat.i_nrows * mat.i_ncols;
                    arma::access::rw(ret.n_cols)=mat.i_ncols;
                    arma::access::rw(ret.n_rows)=mat.i_nrows;

                    mat.i_cscValues.i_ptr=0;
                    mat.i_cscValues.i_nelem=0;
                    mat.i_cscValues.i_owns=false;
                    
                    mat.i_cscColPtr.i_ptr=0;
                    mat.i_cscColPtr.i_nelem=0;
                    mat.i_cscColPtr.i_owns=false;
                    
                    mat.i_cscRowIndex.i_ptr=0;
                    mat.i_cscRowIndex.i_nelem=0;
                    mat.i_cscRowIndex.i_owns=false;
                    
                    mat.i_nrows=0;
                    mat.i_ncols=0;
                    mat.i_nnz=0;
                    
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
            };
            

            template<typename T>
            class convert_method_const_ref<arma::Mat,gMDense,T>{
            public:
                const gMDense<T> operator()(const arma::Mat<T> & mat){
                    gMDense<T> ret;
                    if(mat.n_alloc==0 && mat.mem_state==0){
                        ret.i_values=gMMemory<T>(mat.mem,mat.n_elem);
                    }else{
                        ret.i_values.i_ptr=mat.mem;
                    }
                    ret.i_values.i_nelem=mat.n_elem;
                    ret.i_values.i_owns=true;
                    ret.i_nrows=mat.n_rows;
                    ret.i_ncols=mat.n_cols;
                    arma::access::rw(mat.mem)=nullptr;
                    arma::access::rw(mat.mem_state)=0;
                    arma::access::rw(mat.n_alloc)=0;
                    arma::access::rw(mat.n_cols)=0;
                    arma::access::rw(mat.n_rows)=0;
                    arma::access::rw(mat.n_elem)=0;
                    return std::move(ret);
                }
            };
            
            template<typename T>
            class convert_method_const_ref<gMDense,arma::Mat,T>{
            public:
                const arma::Mat<T> operator()(const gMDense<T> & mat) const{
                    //std::cout << "const ref convert" << std::endl;
                    arma::Mat<T> ret;
                    arma::access::rw(ret.mem)=mat.i_values.i_ptr;
                    arma::access::rw(ret.mem_state)=2;
                    arma::access::rw(ret.n_alloc)=0;
                    arma::access::rw(ret.n_cols)=mat.i_ncols;
                    arma::access::rw(ret.n_elem)=mat.i_values.i_nelem;
                    arma::access::rw(ret.n_rows)=mat.i_nrows;
//                         mat.i_values.i_ptr=nullptr;
//                         mat.i_values.i_nelem=0;
//                         mat.i_values.i_owns=false;
                    return std::move(ret);
                }
            };
            
            template<typename T>
            class convert_method_const_ref<arma::SpMat,gMSparse,T>{
            public:
                const gMSparse<T> operator()(const arma::SpMat<T> & mat) const{
//                     cout << "to Sparse..";cout.flush();
                    gMSparse<T> ret;
                    ret.i_cscValues.release((T *) ret.i_cscValues);
                    ret.i_cscValues.i_ptr=mat.values;
                    ret.i_cscValues.i_nelem=mat.n_nonzero + 1;
                    ret.i_cscValues.i_owns=true;
                    
                    ret.i_cscColPtr.release((indexType *) ret.i_cscColPtr);
                    ret.i_cscColPtr.i_ptr=(indexType *) mat.col_ptrs;
                    ret.i_cscColPtr.i_nelem=mat.n_cols+2;
                    ret.i_cscColPtr.i_owns=true;
                    
                    ret.i_cscRowIndex.release((indexType *) ret.i_cscRowIndex);
                    ret.i_cscRowIndex.i_ptr=(indexType *) mat.row_indices;
                    ret.i_cscRowIndex.i_nelem=mat.n_nonzero + 1;
                    ret.i_cscRowIndex.i_owns=true;

                    ret.i_nrows=mat.n_rows;
                    ret.i_ncols=mat.n_cols;
                    ret.i_nnz=mat.n_nonzero;
                    
                    arma::access::rw(mat.values)=nullptr;
                    arma::access::rw(mat.col_ptrs)=nullptr;
                    arma::access::rw(mat.row_indices)=nullptr;
                    arma::access::rw(mat.n_nonzero)=0;
                    arma::access::rw(mat.n_cols)=0;
                    arma::access::rw(mat.n_rows)=0;
                    arma::access::rw(mat.n_elem)=0;
                    
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
            };
            
            template<typename T>
            class convert_method_const_ref<gMSparse,arma::SpMat,T>{
            public:
                const arma::SpMat<T> operator()(const gMSparse<T> & mat) const{
//                     cout << "to SpMat..";cout.flush();
                    arma::SpMat<T> ret;
                    arma::memory::release(ret.values);
                    arma::memory::release(ret.col_ptrs);
                    arma::memory::release(ret.row_indices);
                    arma::access::rw(ret.values)=mat.i_cscValues.i_ptr;
                    arma::access::rw(ret.col_ptrs)=(arma::uword *) mat.i_cscColPtr.i_ptr;
                    arma::access::rw(ret.row_indices)=(arma::uword *) mat.i_cscRowIndex.i_ptr;
                    arma::access::rw(ret.n_nonzero)=mat.i_nnz;
                    arma::access::rw(ret.n_cols)=mat.i_ncols;
                    arma::access::rw(ret.n_rows)=mat.i_nrows;
                    arma::access::rw(ret.n_elem)=mat.i_nrows * mat.i_ncols;
                    arma::access::rw(ret.n_cols)=mat.i_ncols;
                    arma::access::rw(ret.n_rows)=mat.i_nrows;
                    arma::access::rw(ret.owns)=false;

//                     mat.i_cscValues.i_ptr=0;
//                     mat.i_cscValues.i_nelem=0;
//                     mat.i_cscValues.i_owns=false;
//                     
//                     mat.i_cscColPtr.i_ptr=0;
//                     mat.i_cscColPtr.i_nelem=0;
//                     mat.i_cscColPtr.i_owns=false;
//                     
//                     mat.i_cscRowIndex.i_ptr=0;
//                     mat.i_cscRowIndex.i_nelem=0;
//                     mat.i_cscRowIndex.i_owns=false;
//                     
//                     mat.i_nrows=0;
//                     mat.i_ncols=0;
//                     mat.i_nnz=0;
                    
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
            };

            
            
            template<template<typename> class A,template<typename> class B,typename T>
            B<T> convert(A<T> && mat){
//                 cout << "About to call convert_method" << endl;
                convert_method<A,B,T> method;
                return method(std::move(mat));
            }

            template<template<typename> class A,template<typename> class B,typename T>
            const B<T> convert(const A<T> & mat){
//                 cout << "About to call convert_method" << endl;
                convert_method_const_ref<A,B,T> method;
                return method(mat);
            }
//END armadillo conversion methods


//BEGIN armadillo methods
            template<typename T> 
            arma::Mat<T> initRandom_dense(indexType nrows, indexType ncols,T sparseness){
                return arma::randu<arma::Mat<T>>(nrows,ncols);
            }
            
            template<typename T> 
            arma::SpMat<T> initRandom_sparse(indexType nrows, indexType ncols,T sparseness){
                return arma::sprandu<arma::SpMat<T>>(nrows,ncols,sparseness);
            }

            template<typename T>
            arma::Mat<T> multiply_arma(const arma::Mat<T> & M1,const arma::Mat<T> & M2,bool transposeA,bool transposeB){
                if(&M1==&M2){
                    if(transposeA){
                        return arma::trans(M1) * M1;
                    }else{
                        return M1 * arma::trans(M1);
                    }
                }else{
                    if(transposeA){
                        if(transposeB){
                            return arma::trans(M1) * arma::trans(M2);
                        }else{
                            return arma::trans(M1) * M2;
                        }
                    }else{
                        if(transposeB){
                            return M1 * arma::trans(M2);
                        }else{
                            return M1 * M2;
                        }
                    }
                }
            }

            
            template<typename T>
            arma::SpMat<T> multiply_arma(const arma::Mat<T> & M1,const arma::SpMat<T> & M2,bool transposeA,bool transposeB){
                if(transposeA){
                    if(transposeB){
                        return arma::SpMat<T>(trans(M1) * arma::trans(M2));
                    }else{
                        return arma::SpMat<T>((arma::trans(M1) * M2));
                    }
                }else{
                    if(transposeB){
                        return arma::SpMat<T>(M1 * arma::trans(M2));
                    }else{
                        return arma::SpMat<T>(M1 * M2);
                    }
                }
            }
            
            template<typename T>
            arma::Mat<T> multiply_arma(const arma::SpMat<T> & M1,const arma::Mat<T> & M2,bool transposeA,bool transposeB){
                if(transposeA){
                    if(transposeB){
                        return trans(M1) * trans(M2);
                    }else{
                        return trans(M1) * M2;
                    }
                }else{
                    if(transposeB){
                        return M1 * trans(M2);
                    }else{
                        return M1 * M2;
                    }
                }
            }

            template<typename T>
            arma::SpMat<T> multiply_arma(const arma::SpMat<T> & M1,const arma::SpMat<T> & M2,bool transposeA,bool transposeB){
                if(&M1==&M2){
                    if(transposeA){
                        return trans(M1) * M1;
                    }else{
                        return M1 * trans(M1);
                    }
                }else{
                    if(transposeA){
                        if(transposeB){
                            return trans(M1) * trans(M2);
                        }else{
                            return trans(M1) * M2;
                        }
                    }else{
                        if(transposeB){
                            return M1 * trans(M2);
                        }else{
                            return M1 * M2;
                        }
                    }
                }
            }

            template<typename T>
            arma::Mat<T> multiply_transpose_arma(const arma::Mat<T> & M, bool transpose){
                return multiply_arma(M,M,transpose,!transpose);
            }

            template<typename T>
            arma::SpMat<T> multiply_transpose_arma(const arma::SpMat<T> & M, bool transpose){
                return multiply_arma(M,M,transpose,!transpose);
            }
//END armadillo methods            
            
        }
    }
}



