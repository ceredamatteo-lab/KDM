#pragma once
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "gMethods.h"
#include <sstream>
#include <cuda_runtime_api.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#define THRUST_IGNORE_CUB_VERSION_CHECK

namespace geco{
	namespace methods{
        namespace cuda{
            
            typedef uint32_t cudaIndexType;
            //typedef indexType cudaIndexType;

//BEGIN cuda checks
            inline void CHECK_CUDART_ERROR(cudaError_t error){
                if(error!=cudaSuccess){
                    std::stringstream msg;
                    msg << "CUDA ERROR: " << cudaGetErrorString(error) << "\t" << cudaGetErrorName (error); 
                    throw gMethodException(msg.str());
                }
            }

            inline void CHECK_CUSPARSE_ERROR(cusparseStatus_t status) {
                if (status != CUSPARSE_STATUS_SUCCESS) {
                    std::stringstream msg;	
                    msg << "CUSPARSE ERROR: " << cusparseGetErrorString(status) << "\t" << cusparseGetErrorName(status);
                    throw gMethodException(msg.str());
                }
            }

            inline void CHECK_CUBLAS_ERROR(cublasStatus_t status) {
                if (status != CUBLAS_STATUS_SUCCESS) {
                    std::stringstream msg;
                    switch(status){
                        case CUBLAS_STATUS_NOT_INITIALIZED: msg << "The cuBLAS library was not initialized. This is usually caused by the lack of a prior cublasCreate() call, an error in the CUDA Runtime API called by the cuBLAS routine, or an error in the hardware setup.";break;
                        case CUBLAS_STATUS_ALLOC_FAILED: msg << "Resource allocation failed inside the cuBLAS library. This is usually caused by a cudaMalloc() failure.";break;
                        case CUBLAS_STATUS_INVALID_VALUE: msg << "An unsupported value or parameter was passed to the function (a negative vector size, for example).";break;
                        case CUBLAS_STATUS_ARCH_MISMATCH: msg << "The function requires a feature absent from the device architecture; usually caused by the lack of support for double precision.";break;
                        case CUBLAS_STATUS_MAPPING_ERROR: msg << "An access to GPU memory space failed, which is usually caused by a failure to bind a texture.";
                        case CUBLAS_STATUS_EXECUTION_FAILED: msg << "The GPU program failed to execute. This is often caused by a launch failure of the kernel on the GPU, which can be caused by multiple reasons.";break;
                        case CUBLAS_STATUS_INTERNAL_ERROR: msg << "An internal cuBLAS operation failed. This error is usually caused by a cudaMemcpyAsync() failure.";break;
                        case CUBLAS_STATUS_NOT_SUPPORTED: msg << "The functionnality requested is not supported";break;
                        case CUBLAS_STATUS_LICENSE_ERROR: msg << "The functionnality requested requires some license and an error was detected when trying to check the current licensing. This error can happen if the license is not present or is expired or if the environment variable NVIDIA_LICENSE_FILE is not set properly";break;
                        default: break;
                    }
                    throw gMethodException(msg.str());
                }
            }

            inline void CHECK_CUSOLVER_ERROR(cusolverStatus_t status,int devInfo){
                if (status != CUSOLVER_STATUS_SUCCESS) {
                    std::stringstream msg;
                    switch(status){
                        case CUSOLVER_STATUS_NOT_INITIALIZED: msg << "The cuSolver library was not initialized. This is usually caused by the lack of a prior call, an error in the CUDA Runtime API called by the cuSolver routine, or an error in the hardware setup.";break;
                        case CUSOLVER_STATUS_ALLOC_FAILED: msg << "Resource allocation failed inside the cuSolver library. This is usually caused by a cudaMalloc() failure.";break;
                        case CUSOLVER_STATUS_INVALID_VALUE: msg << "An unsupported value or parameter was passed to the function (a negative vector size, for example).";break;
                        case CUSOLVER_STATUS_ARCH_MISMATCH: msg << "The function requires a feature absent from the device architecture; usually caused by the lack of support for atomic operations or double precision.";break;
                        case CUSOLVER_STATUS_EXECUTION_FAILED: msg << "he GPU program failed to execute. This is often caused by a launch failure of the kernel on the GPU, which can be caused by multiple reasons.";break;
                        case CUSOLVER_STATUS_INTERNAL_ERROR: msg << "An internal cuSolver operation failed. This error is usually caused by a cudaMemcpyAsync() failure.";break;
                        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED: msg << "The matrix type is not supported by this function. This is usually caused by passing an invalid matrix descriptor to the function.";break;
                        default: break;
                    }
                    msg << " devinfo=" << devInfo;
                    throw gMethodException(msg.str());
                }
            }
//END cuda checks

//BEGIN templated cuda calls
            template<typename F1, typename F2, typename ... Args>
            auto invokeOne(F1 f1, F2 f2, Args && ... args) -> decltype(f1(args...)){
                return f1(args...);
            }

            template<typename F1, typename F2, typename ... Args>
            auto invokeOne(F1 f1, F2 f2, Args && ... args) -> decltype(f2(args...)){
                return f2(args...);
            }

            template<typename T>
            cublasStatus_t cublasTgeam(
                cublasHandle_t handle,                           
                cublasOperation_t transa, 
                cublasOperation_t transb,
                int m, 
                int n,
                const T *alpha, /* host or device pointer */ 
                const T *A, 
                int lda,
                const T *beta , /* host or device pointer */ 
                const T *B, 
                int ldb,
                T *C, 
                int ldc
            ){
                return invokeOne(cublasSgeam,cublasDgeam,handle,transa,transb,m,n,alpha,A,lda,beta,B,ldb,C,ldc);
            }

            template<typename T>
            cublasStatus_t cublasTgemm(
                cublasHandle_t handle,
                cublasOperation_t transa, 
                cublasOperation_t transb,
                int m, 
                int n,
                int k,
                const T * alpha,
                const T * A, 
                int lda,
                const T * B, 
                int ldb,
                const T * beta,
                T * C, 
                int ldc
            ){
                return invokeOne(cublasSgemm,cublasDgemm,handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
            }

            template<typename T>
            cublasStatus_t cublasTsyrk(
                cublasHandle_t      handle,
                cublasFillMode_t    uplo, 
                cublasOperation_t   trans,
                int                 n, 
                int                 k,
                const T *           alpha,
                const T *           A, 
                int                 lda,
                const T *           beta,
                T *                 C, 
                int                 ldc
            ){
                return invokeOne(cublasSsyrk,cublasDsyrk,handle,uplo,trans,n,k,alpha,A,lda,beta,C,ldc); 
            }
            
            template<typename T>
            cusolverStatus_t cusolverDnTsyevd_bufferSize(
                cusolverDnHandle_t  handle,
                cusolverEigMode_t   jobz,
                cublasFillMode_t    uplo,
                int                 n,
                const T *           A,
                int                 lda,
                const T *           W,
                int *               lwork
            ){
                return invokeOne(cusolverDnSsyevd_bufferSize,cusolverDnDsyevd_bufferSize,handle,jobz,uplo,n,A,lda,W,lwork);
            }
            
            template<typename T>
            cusolverStatus_t cusolverDnTsyevd(
                cusolverDnHandle_t  handle,
                cusolverEigMode_t   jobz,
                cublasFillMode_t    uplo,
                int                 n,
                T *                 A,
                int                 lda,
                T *                 W,
                T *                 work,
                int                 lwork,
                int *               devInfo
            ){
                return invokeOne(cusolverDnSsyevd,cusolverDnDsyevd,handle,jobz,uplo,n,A,lda,W,work,lwork,devInfo);
            }
//END templated cuda calls
            
//BEGIN autoHandles
            template<typename T>
            class cudaXHandle{
            private:
                T handle=nullptr;
                
                template<typename F1, typename F2, typename F3>
                auto cudaCreateXHandle(F1 f1, F2 f2, F3 f3, T & handle) -> decltype(CHECK_CUBLAS_ERROR(f1(&handle))){
                    CHECK_CUBLAS_ERROR(f1(&handle));
                }
                
                template<typename F1, typename F2, typename F3>
                auto cudaCreateXHandle(F1 f1, F2 f2, F3 f3, T & handle) -> decltype(CHECK_CUSPARSE_ERROR(f2(&handle))){
                    CHECK_CUSPARSE_ERROR(f2(&handle));
                }
                
                template<typename F1, typename F2, typename F3>
                auto cudaCreateXHandle(F1 f1, F2 f2, F3 f3, T & handle) -> decltype(CHECK_CUSOLVER_ERROR(f3(&handle),0)){
                    CHECK_CUSOLVER_ERROR(f3(&handle),0);
                }
                
                
                template<typename F1, typename F2, typename F3>
                auto cudaDestroyXHandle(F1 f1, F2 f2, F3 f3, T & handle) -> decltype(f1(handle)){
                    return f1(handle);
                }
                
                template<typename F1, typename F2, typename F3>
                auto cudaDestroyXHandle(F1 f1, F2 f2, F3 f3, T & handle) -> decltype(f2(handle)){
                    return f2(handle);
                }
                
                template<typename F1, typename F2, typename F3>
                auto cudaDestroyXHandle(F1 f1, F2 f2, F3 f3, T & handle) -> decltype(f3(handle)){
                    return f3(handle);
                }
                
            public:
                cudaXHandle(){
                }
                
                cudaXHandle(const cudaXHandle & other):handle(other.handle){
                }
                
                ~cudaXHandle(){
                    // no check, destructors must not throw
                    if(handle!=nullptr){
                        cudaDestroyXHandle(cublasDestroy,cusparseDestroy,cusolverDnDestroy,handle);
                    }
                }
                
                operator T (){
                    if(handle==nullptr){
                        cudaCreateXHandle(cublasCreate,cusparseCreate,cusolverDnCreate,handle);
                    }
                    return handle;
                }
                
            };
            
//             typedef cudaXHandle<cublasHandle_t> autoCublasHandle;
//             typedef cudaXHandle<cusparseHandle_t> autoCusparseHandle;
//             typedef cudaXHandle<cusolverDnHandle_t> autoCusolverHandle;
            
            class autoCudaHandles{
            private:
                cudaXHandle<cublasHandle_t> i_cublasHandle;
                cudaXHandle<cusparseHandle_t> i_cusparseHandle;
                cudaXHandle<cusolverDnHandle_t> i_cusolverHandle;
            public:
                operator cublasHandle_t () {
                    return i_cublasHandle;
                }
                
                operator cusparseHandle_t () {
                    return i_cusparseHandle;
                }
                
                operator cusolverDnHandle_t () {
                    return i_cusolverHandle;
                }
            };
//END autoHandles
            
//BEGIN cuda matrix classes
            template<typename T>
            class cudaMemory{
                T * i_ptr=nullptr;
            public:
                cudaMemory(){
                    i_ptr=nullptr;
                }
                
                cudaMemory(indexType n_elem,bool init0=false):cudaMemory(){
                    if(n_elem>0){
                        //cout << "CUDA ALLOC" << endl;
                        CHECK_CUDART_ERROR(cudaMalloc((void **)&i_ptr, n_elem * sizeof(T)));
                        if(init0){
                            CHECK_CUDART_ERROR(cudaMemset((void *)i_ptr, 0,n_elem * sizeof(T)));
                        }
                        CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                    }                        
                }
                
                template<typename D>
                cudaMemory(const D * hostPtr,indexType n_elem):cudaMemory(){
                    copyFrom(hostPtr,n_elem);
                }
                
                cudaMemory(cudaMemory<T> && other):cudaMemory(){
                    swap(*this,other);
                }
                
                cudaMemory<T> & operator = (cudaMemory<T> && other){
                    swap(*this,other);
                    return *this;
                }
                
                
                ~cudaMemory(){
                    if(i_ptr!=nullptr){
                        // no check, destructors must not throw
                        cudaFree(i_ptr);
                        //CHECK_CUDART_ERROR(cudaFree(i_ptr));
                    }
                }
                
                operator T * (){
                    return i_ptr;
                }
                
                operator const T * () const {
                    return i_ptr;
                }
                
                template<typename D>
                void copyFrom(const D * hostPtr,indexType n_elem){
                    if(n_elem>0 && hostPtr!=nullptr){
                        if(i_ptr==nullptr){
                            CHECK_CUDART_ERROR(cudaMalloc((void **)&i_ptr, n_elem * sizeof(T)));
                            CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                        }
                        if(typeid(D)==typeid(T)){
                            CHECK_CUDART_ERROR(cudaMemcpy(i_ptr,hostPtr,n_elem *sizeof(T),cudaMemcpyHostToDevice));
                            CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                        }else{
                            T * tmpHostPtr = new T[n_elem];
//apparently std::copy is faster
#pragma omp parallel for
                            for(indexType i=0;i<n_elem;i++){
                                tmpHostPtr[i]=hostPtr[i];
                            }
//                            std::copy(hostPtr,hostPtr+n_elem,tmpHostPtr);
                            CHECK_CUDART_ERROR(cudaMemcpy(i_ptr,tmpHostPtr,n_elem *sizeof(T),cudaMemcpyHostToDevice));
                            CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                            delete [] tmpHostPtr;
                        }
                    }
                }
                
                friend void swap(cudaMemory<T> & a,cudaMemory<T> & b){
                    std::swap(a.i_ptr,b.i_ptr);
                }
                
                
            };
            
            template<typename T> 
            class cudaDense{
            private:
                void updateDescr(){
                    if(i_descr!=nullptr){
                        cusparseDestroyDnMat(i_descr);
                        cudaDeviceSynchronize();                        
                        i_descr=nullptr;
                    }
                    CHECK_CUSPARSE_ERROR(cusparseCreateDnMat(
                        &i_descr,
                        i_nrows,
                        i_ncols,
                        i_nrows,
                        i_values,
                        (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                        CUSPARSE_ORDER_COL
                    ));
                    cudaDeviceSynchronize();
                }
                
            public:
                cudaIndexType i_nrows;
                cudaIndexType i_ncols;
                cusparseDnMatDescr_t i_descr;                
                cudaMemory<T> i_values;
                
                cudaDense():i_nrows(0),i_ncols(0),i_descr(nullptr){
                }
                
                cudaDense(indexType nrows,indexType ncols):i_nrows(nrows),i_ncols(ncols),i_descr(nullptr),i_values(nrows*ncols){
                    updateDescr();
                }
                
                cudaDense(indexType nrows,indexType ncols, T val):i_nrows(nrows),i_ncols(ncols),i_descr(nullptr),i_values(nrows*ncols,true){
                    updateDescr();
                    //cuda::elementwise_add(*this,val,*this);
                }

                
                //copy
                cudaDense(const cudaDense & other):cudaDense(other.i_nrows,other.i_ncols){
                    CHECK_CUDART_ERROR(cudaMemcpy(i_values, other.i_values, i_nrows * i_ncols * sizeof(T), cudaMemcpyDeviceToDevice));
                    cudaDeviceSynchronize();
                }

                cudaDense & operator = (cudaDense other){
                    swap(*this,other);
                    return *this;
                }

                
                //move
                cudaDense(cudaDense && other):cudaDense<T>(){
                    swap(*this,other);
                }

//                 cudaDense & operator = (cudaDense && other){
//                     swap(*this,other);
//                     return *this;
//                 }
                
                //destructor
                ~cudaDense(){
                    //no check, destructors cannot throw
                    cusparseDestroyDnMat(i_descr);
                    cudaDeviceSynchronize();
                    //CHECK_CUSPARSE_ERROR(cusparseDestroyDnMat(i_descr));
                }
                
                void updateNrows(cudaIndexType nrows){
                    i_nrows=nrows;
                    updateDescr();
                }
                
                void updateNcols(cudaIndexType ncols){
                    i_ncols=ncols;
                    updateDescr();
                }
                


                //inplace operations
                cudaDense & inplace_transpose(autoCudaHandles & cudaHandles){
                    cudaMemory<T> tmp_values(i_nrows * i_ncols);
                    CHECK_CUDART_ERROR(cudaMemcpy((T*) tmp_values, i_values, i_nrows * i_ncols * sizeof(T), cudaMemcpyHostToDevice));
                    T alpha=1,beta=0;
                    CHECK_CUBLAS_ERROR(cublasTgeam<T>(
                        cudaHandles,
                        CUBLAS_OP_T, 
                        CUBLAS_OP_N,
                        i_nrows, 
                        i_ncols,
                        &alpha,
                        (T*) tmp_values, 
                        i_nrows,
                        &beta,
                        nullptr, 
                        0,
                        (T*) i_values,
                        i_nrows
                    ));
                    cudaDeviceSynchronize();
                    std::swap(i_nrows,i_ncols);
                    return *this;
                }

                
//                 //convert to armadillo
//                 operator arma::Mat<T> () const{
//                     arma::Mat<T> ret(i_nrows,i_ncols);
//                     CHECK_CUDART_ERROR(cudaMemcpy(ret.memptr(), i_values, i_nrows * i_ncols * sizeof(T), cudaMemcpyDeviceToHost));
//                     return ret;
//                 }
                
                friend void swap(cudaDense & a,cudaDense & b){
                    using std::swap;
                    swap(a.i_nrows,b.i_nrows);
                    swap(a.i_ncols,b.i_ncols);
                    swap(a.i_values,b.i_values);
                    swap(a.i_descr,b.i_descr);
                }

                
                
                friend std::ostream & operator << (std::ostream & out, const cudaDense<T> & mat){
                    T * values= new T[mat.i_nrows*mat.i_ncols];
                    CHECK_CUDART_ERROR(cudaMemcpy(values, mat.i_values, mat.i_nrows * mat.i_ncols * sizeof(T), cudaMemcpyDeviceToHost));
                    
                    for(indexType row=0;row<mat.i_nrows;row++){
                        for(indexType col=0; col<mat.i_ncols;col++){
                            out << "\t" << values[row + col*mat.i_nrows];
                        }
                        out << std::endl;
                    }
                    delete [] values;
                    return out;
                }
                
            };
            
            template<typename T>
            class cudaSparse{
            public:
                void updateDescr(){
                    if(i_descr!=nullptr){
                        cusparseDestroySpMat(i_descr);
                        cudaDeviceSynchronize();
                        i_descr=nullptr;
                    }
                    CHECK_CUSPARSE_ERROR(cusparseCreateCsr(
                        &i_descr,
                        i_nrows,
                        i_ncols,
                        i_nnz,
                        i_csrRowPtr,
                        i_csrColInd,
                        i_csrValues,
                        CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO,
                        (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F)
                    ));
                    CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                }
            public:
                cudaIndexType i_nrows;
                cudaIndexType i_ncols;
                cudaIndexType i_nnz;
                
                cudaMemory<T> i_csrValues;
                cudaMemory<cudaIndexType> i_csrRowPtr;
                cudaMemory<cudaIndexType> i_csrColInd;
                
                cusparseSpMatDescr_t i_descr;                
                
                cudaSparse():i_nrows(0),i_ncols(0),i_nnz(0),i_csrValues(),i_csrRowPtr(),i_csrColInd(),i_descr(nullptr){
                }

                cudaSparse(indexType nrows,indexType ncols,indexType nnz):i_nrows(nrows),i_ncols(ncols),i_nnz(nnz),i_csrValues(nnz),i_csrRowPtr(nrows+1),i_csrColInd(nnz),i_descr(nullptr){
                    CHECK_CUDART_ERROR(cudaMemset(i_csrRowPtr,0,(i_nrows+1) * sizeof(cudaIndexType)));
                    CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                    updateDescr();
                }
                
                //copy
                cudaSparse(const cudaSparse & other):cudaSparse(other.i_nrows,other.i_ncols,other.i_nnz){
                    CHECK_CUDART_ERROR(cudaMemcpy(i_csrValues, other.i_csrValues, i_nnz * sizeof(T), cudaMemcpyDeviceToDevice));
                    CHECK_CUDART_ERROR(cudaMemcpy(i_csrRowPtr, other.i_csrRowPtr, (i_nrows + 1) * sizeof(cudaIndexType), cudaMemcpyDeviceToDevice));
                    CHECK_CUDART_ERROR(cudaMemcpy(i_csrColInd, other.i_csrColInd, i_nnz * sizeof(cudaIndexType), cudaMemcpyDeviceToDevice));
                    cudaDeviceSynchronize();
                }
                
                //move
                cudaSparse(cudaSparse && other):cudaSparse<T>(){
                    swap(*this,other);
                }
                
                cudaSparse & operator = (cudaSparse && other){
                    swap(*this,other);
                    return *this;
                }
                
                //destructor
                ~cudaSparse(){
                    cusparseDestroySpMat(i_descr);
                    cudaDeviceSynchronize();
                };

                
                //inplace operations
                cudaSparse & inplace_transpose();
                
                friend void swap(cudaSparse & a,cudaSparse & b){
                    using std::swap;
                    swap(a.i_nrows,b.i_nrows);
                    swap(a.i_ncols,b.i_ncols);
                    swap(a.i_csrValues,b.i_csrValues);
                    swap(a.i_csrRowPtr,b.i_csrRowPtr);
                    swap(a.i_csrColInd,b.i_csrColInd);
                    swap(a.i_nnz,b.i_nnz);
                    swap(a.i_descr,b.i_descr);
                }
                
            };
            
//END cuda matrix classes
            
//BEGIN cuda conversion methods
            template<template<typename> class A, template<typename> class B, typename T>
            class convert_method{
            public:
                B<T> operator()(A<T> toConvert, bool transpose,autoCudaHandles & cudaHandles) const;
                //B<T> operator()(const A<T> & toConvert, bool transpose,autoCudaHandles & cudaHandles) const;
            };
            
            template<typename T>
            class convert_method<cudaDense,gMDense,T>{
            public:
                gMDense<T> operator()(cudaDense<T> mat, bool transpose,autoCudaHandles & cudaHandles){
//                     cout << "to Dense..";cout.flush();
                    if(transpose){
                        mat.inplace_transpose(cudaHandles);
                    }
                    gMDense<T> ret(mat.i_nrows,mat.i_ncols);
                    CHECK_CUDART_ERROR(cudaMemcpy(( T*)ret.i_values, mat.i_values, mat.i_nrows * mat.i_ncols * sizeof(T), cudaMemcpyDeviceToHost));
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
                
//                 gMDense<T> operator()(const cudaDense<T> & mat, bool transpose,autoCudaHandles & cudaHandles){
// //                     cout << "to Dense..";cout.flush();
//                     if(transpose){
//                         mat.inplace_transpose(cudaHandles);
//                     }
//                     gMDense<T> ret(mat.i_nrows,mat.i_ncols);
//                     CHECK_CUDART_ERROR(cudaMemcpy(( T*)ret.i_values, mat.i_values, mat.i_nrows * mat.i_ncols * sizeof(T), cudaMemcpyDeviceToHost));
// //                     cout << "returning" << endl;
//                     return std::move(ret);
//                 }
                
            };
            
            template<typename T>
            class convert_method<gMDense,cudaDense,T>{
            public:
                cudaDense<T> operator()(gMDense<T> mat,bool transpose,autoCudaHandles & cudaHandles) const{
                    cudaDense<T> ret((transpose)?(mat.i_ncols):(mat.i_nrows),(transpose)?(mat.i_nrows):(mat.i_ncols));
                    if(transpose){
                        cudaMemory<T> tmp((T *) mat.i_values,mat.i_nrows * mat.i_ncols);
                        T alpha=1,beta=0;
                        CHECK_CUBLAS_ERROR(cublasTgeam<T>(
                            cudaHandles,
                            CUBLAS_OP_T, 
                            CUBLAS_OP_N,
                            ret.i_nrows, 
                            ret.i_ncols,
                            &alpha,
                            tmp, 
                            ret.i_ncols,
                            &beta,
                            ret.i_values,
                            ret.i_nrows,
                            ret.i_values,
                            ret.i_nrows
                        ));
                    }else{
                        ret.i_values.copyFrom((T *) mat.i_values,mat.i_nrows * mat.i_ncols);
                    }
                    CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                    return std::move(ret);
                }
                
//                 cudaDense<T> operator()(const gMDense<T> & mat,bool transpose,autoCudaHandles & cudaHandles) const{
//                     cudaDense<T> ret((transpose)?(mat.i_ncols):(mat.i_nrows),(transpose)?(mat.i_nrows):(mat.i_ncols));
//                     if(transpose){
//                         cudaMemory<T> tmp((const T *) mat.i_values,mat.i_nrows * mat.i_ncols);
//                         T alpha=1,beta=0;
//                         CHECK_CUBLAS_ERROR(cublasTgeam<T>(
//                             cudaHandles,
//                             CUBLAS_OP_T, 
//                             CUBLAS_OP_N,
//                             ret.i_nrows, 
//                             ret.i_ncols,
//                             &alpha,
//                             tmp, 
//                             ret.i_ncols,
//                             &beta,
//                             ret.i_values,
//                             ret.i_nrows,
//                             ret.i_values,
//                             ret.i_nrows
//                         ));
//                     }else{
//                         ret.i_values.copyFrom((const T *) mat.i_values,mat.i_nrows * mat.i_ncols);
//                     }
//                     CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                     return std::move(ret);
//                 }
                
            };
            
            template<typename T>
            class convert_method<cudaSparse,gMSparse,T>{
            public:
                gMSparse<T> operator()(cudaSparse<T> mat,bool transpose,autoCudaHandles & cudaHandles) const{
//                     cout << "to Sparse..";cout.flush();
                    
                    gMSparse<T> ret(transpose?mat.i_ncols:mat.i_nrows,transpose?mat.i_nrows:mat.i_ncols,mat.i_nnz);
                    //ret.reserve((transpose)?(i_ncols):(i_nrows),(transpose)?(i_nrows):(i_ncols),(indexType) i_nnz);
                    if(transpose){
                        CHECK_CUDART_ERROR(cudaMemcpy((void *) (T *) ret.i_cscValues,(void *) mat.i_csrValues,mat.i_nnz*sizeof(T), cudaMemcpyDeviceToHost));
                        
                        if(typeid(cudaIndexType)==typeid(indexType)){
                            CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscColPtr,mat.i_csrRowPtr,(mat.i_nrows+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                            CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscRowIndex,mat.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                        }else{
                            cudaIndexType * tmp = new cudaIndexType[mat.i_nrows+1];
                            CHECK_CUDART_ERROR(cudaMemcpy(tmp,mat.i_csrRowPtr,(mat.i_nrows+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                            std::copy(tmp,tmp+mat.i_nrows+1,(indexType *) ret.i_cscColPtr);
                            delete [] tmp ;
                            tmp = new cudaIndexType[mat.i_nnz];
                            CHECK_CUDART_ERROR(cudaMemcpy(tmp,mat.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                            std::copy(tmp,tmp+mat.i_nnz,(indexType *) ret.i_cscRowIndex);
                            delete [] tmp;
                        }
                    }else{
                        cudaSparse<T> A(mat.i_ncols,mat.i_nrows,mat.i_nnz);

                        size_t bufferSize;
                        CHECK_CUSPARSE_ERROR( cusparseCsr2cscEx2_bufferSize(
                            cudaHandles,
                            mat.i_nrows,
                            mat.i_ncols,
                            mat.i_nnz,
                            mat.i_csrValues,
                            (const int *) (const cudaIndexType *) mat.i_csrRowPtr,
                            (const int *)  (const cudaIndexType *) mat.i_csrColInd,
                            A.i_csrValues,
                            (int *)  (cudaIndexType *) A.i_csrRowPtr,
                            (int *)  (cudaIndexType *) A.i_csrColInd,
                            (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                            CUSPARSE_ACTION_NUMERIC,
                            CUSPARSE_INDEX_BASE_ZERO,
                            CUSPARSE_CSR2CSC_ALG1,
                            &bufferSize
                        ));
                        CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                        
                        cudaMemory<unsigned char> buffer(bufferSize);
                        
                        CHECK_CUSPARSE_ERROR( cusparseCsr2cscEx2(
                            cudaHandles,
                            mat.i_nrows,
                            mat.i_ncols,
                            mat.i_nnz,
                            mat.i_csrValues,
                            (const int *)  (const cudaIndexType *) mat.i_csrRowPtr,
                            (const int *)  (const cudaIndexType *) mat.i_csrColInd,
                            A.i_csrValues,
                            (int *)  (cudaIndexType *) A.i_csrRowPtr,
                            (int *)  (cudaIndexType *) A.i_csrColInd,
                            (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                            CUSPARSE_ACTION_NUMERIC,
                            CUSPARSE_INDEX_BASE_ZERO,
                            CUSPARSE_CSR2CSC_ALG1,
                            (void *) buffer
                        ));
                        CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                        
                        CHECK_CUDART_ERROR(cudaMemcpy((T *) ret.i_cscValues,A.i_csrValues,mat.i_nnz*sizeof(T), cudaMemcpyDeviceToHost));

                        if(typeid(cudaIndexType)==typeid(indexType)){
                            CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscColPtr,A.i_csrRowPtr,(mat.i_ncols+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                            CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscRowIndex,A.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                        }else{
                            cudaIndexType * tmp = new cudaIndexType[mat.i_ncols+1];
                            CHECK_CUDART_ERROR(cudaMemcpy(tmp,A.i_csrRowPtr,(mat.i_ncols+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                            std::copy(tmp,tmp+mat.i_ncols+1,(indexType *) ret.i_cscColPtr);
                            delete [] tmp ;
                            tmp = new cudaIndexType[mat.i_nnz];
                            CHECK_CUDART_ERROR(cudaMemcpy(tmp,A.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
                            std::copy(tmp,tmp+mat.i_nnz,(indexType *) ret.i_cscRowIndex);
                            delete [] tmp;
                        }
                    }
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
                
//                 gMSparse<T> operator()(const cudaSparse<T> & mat,bool transpose,autoCudaHandles & cudaHandles) const{
// //                     cout << "to Sparse..";cout.flush();
//                     
//                     gMSparse<T> ret(transpose?mat.i_ncols:mat.i_nrows,transpose?mat.i_nrows:mat.i_ncols,mat.i_nnz);
//                     //ret.reserve((transpose)?(i_ncols):(i_nrows),(transpose)?(i_nrows):(i_ncols),(indexType) i_nnz);
//                     if(transpose){
//                         CHECK_CUDART_ERROR(cudaMemcpy((void *) (T *) ret.i_cscValues,(void *) mat.i_csrValues,mat.i_nnz*sizeof(T), cudaMemcpyDeviceToHost));
//                         
//                         if(typeid(cudaIndexType)==typeid(indexType)){
//                             CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscColPtr,mat.i_csrRowPtr,(mat.i_nrows+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                             CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscRowIndex,mat.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                         }else{
//                             cudaIndexType * tmp = new cudaIndexType[mat.i_nrows+1];
//                             CHECK_CUDART_ERROR(cudaMemcpy(tmp,mat.i_csrRowPtr,(mat.i_nrows+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                             std::copy(tmp,tmp+mat.i_nrows+1,(indexType *) ret.i_cscColPtr);
//                             delete [] tmp ;
//                             tmp = new cudaIndexType[mat.i_nnz];
//                             CHECK_CUDART_ERROR(cudaMemcpy(tmp,mat.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                             std::copy(tmp,tmp+mat.i_nnz,(indexType *) ret.i_cscRowIndex);
//                             delete [] tmp;
//                         }
//                     }else{
//                         cudaSparse<T> A(mat.i_ncols,mat.i_nrows,mat.i_nnz);
// 
//                         size_t bufferSize;
//                         CHECK_CUSPARSE_ERROR( cusparseCsr2cscEx2_bufferSize(
//                             cudaHandles,
//                             mat.i_nrows,
//                             mat.i_ncols,
//                             mat.i_nnz,
//                             mat.i_csrValues,
//                             (const int *) (const cudaIndexType *) mat.i_csrRowPtr,
//                             (const int *)  (const cudaIndexType *) mat.i_csrColInd,
//                             A.i_csrValues,
//                             (int *)  (cudaIndexType *) A.i_csrRowPtr,
//                             (int *)  (cudaIndexType *) A.i_csrColInd,
//                             (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
//                             CUSPARSE_ACTION_NUMERIC,
//                             CUSPARSE_INDEX_BASE_ZERO,
//                             CUSPARSE_CSR2CSC_ALG1,
//                             &bufferSize
//                         ));
//                         CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                         
//                         cudaMemory<unsigned char> buffer(bufferSize);
//                         
//                         CHECK_CUSPARSE_ERROR( cusparseCsr2cscEx2(
//                             cudaHandles,
//                             mat.i_nrows,
//                             mat.i_ncols,
//                             mat.i_nnz,
//                             mat.i_csrValues,
//                             (const int *)  (const cudaIndexType *) mat.i_csrRowPtr,
//                             (const int *)  (const cudaIndexType *) mat.i_csrColInd,
//                             A.i_csrValues,
//                             (int *)  (cudaIndexType *) A.i_csrRowPtr,
//                             (int *)  (cudaIndexType *) A.i_csrColInd,
//                             (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
//                             CUSPARSE_ACTION_NUMERIC,
//                             CUSPARSE_INDEX_BASE_ZERO,
//                             CUSPARSE_CSR2CSC_ALG1,
//                             (void *) buffer
//                         ));
//                         CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                         
//                         CHECK_CUDART_ERROR(cudaMemcpy((T *) ret.i_cscValues,A.i_csrValues,mat.i_nnz*sizeof(T), cudaMemcpyDeviceToHost));
// 
//                         if(typeid(cudaIndexType)==typeid(indexType)){
//                             CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscColPtr,A.i_csrRowPtr,(mat.i_ncols+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                             CHECK_CUDART_ERROR(cudaMemcpy((indexType *) ret.i_cscRowIndex,A.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                         }else{
//                             cudaIndexType * tmp = new cudaIndexType[mat.i_ncols+1];
//                             CHECK_CUDART_ERROR(cudaMemcpy(tmp,A.i_csrRowPtr,(mat.i_ncols+1)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                             std::copy(tmp,tmp+mat.i_ncols+1,(indexType *) ret.i_cscColPtr);
//                             delete [] tmp ;
//                             tmp = new cudaIndexType[mat.i_nnz];
//                             CHECK_CUDART_ERROR(cudaMemcpy(tmp,A.i_csrColInd,(mat.i_nnz)*sizeof(cudaIndexType), cudaMemcpyDeviceToHost)); 
//                             std::copy(tmp,tmp+mat.i_nnz,(indexType *) ret.i_cscRowIndex);
//                             delete [] tmp;
//                         }
//                     }
// //                     cout << "returning" << endl;
//                     return std::move(ret);
//                 }
                
            };
            
            template<typename T>
            class convert_method<gMSparse,cudaSparse,T>{
            public:
                cudaSparse<T> operator()(gMSparse<T> mat,bool transpose,autoCudaHandles & cudaHandles) const{
//                     cout << "to SpMat..";cout.flush();
                    cudaSparse<T> ret( (transpose)?(mat.i_ncols):(mat.i_nrows),(transpose)?(mat.i_nrows):(mat.i_ncols),mat.i_nnz);
                    if(!transpose){                     
                        cudaSparse<T> csrTrans(mat.i_ncols,mat.i_nrows,mat.i_nnz);
                        csrTrans.i_csrValues.copyFrom((T*) mat.i_cscValues,mat.i_nnz);
                        csrTrans.i_csrRowPtr.copyFrom((indexType *) mat.i_cscColPtr,mat.i_ncols+1);
                        csrTrans.i_csrColInd.copyFrom((indexType *) mat.i_cscRowIndex,mat.i_nnz);
                     
                        size_t bufferSize;
                        CHECK_CUSPARSE_ERROR( cusparseCsr2cscEx2_bufferSize(
                            cudaHandles,
                            ret.i_ncols,
                            ret.i_nrows,
                            ret.i_nnz,
                            csrTrans.i_csrValues,
                            (const int *) (const cudaIndexType *) csrTrans.i_csrRowPtr,
                            (const int *) (const cudaIndexType *) csrTrans.i_csrColInd,
                            ret.i_csrValues,
                            (int *) (cudaIndexType *) ret.i_csrRowPtr,
                            (int *) (cudaIndexType *) ret.i_csrColInd,
                            (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                            CUSPARSE_ACTION_NUMERIC,
                            CUSPARSE_INDEX_BASE_ZERO,
                            CUSPARSE_CSR2CSC_ALG1,
                            &bufferSize
                        ));
                        CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                        
                        cudaMemory<unsigned char> buffer(bufferSize);
                    
                        CHECK_CUSPARSE_ERROR(cusparseCsr2cscEx2(
                            cudaHandles,
                            ret.i_ncols,
                            ret.i_nrows,
                            ret.i_nnz,
                            csrTrans.i_csrValues,
                            (const int *) (const cudaIndexType *) csrTrans.i_csrRowPtr,
                            (const int *) (const cudaIndexType *) csrTrans.i_csrColInd,
                            ret.i_csrValues,
                            (int *) (cudaIndexType *) ret.i_csrRowPtr,
                            (int *) (cudaIndexType *) ret.i_csrColInd,
                            (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                            CUSPARSE_ACTION_NUMERIC,
                            CUSPARSE_INDEX_BASE_ZERO,
                            CUSPARSE_CSR2CSC_ALG1,
                            (void *) buffer
                        ));
                        CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                    }else{
                        ret.i_csrValues.copyFrom((T*) mat.i_cscValues,mat.i_nnz);
                        ret.i_csrRowPtr.copyFrom((indexType *) mat.i_cscColPtr,mat.i_ncols+1);
                        ret.i_csrColInd.copyFrom((indexType *) mat.i_cscRowIndex,mat.i_nnz);
                    }
//                     cout << "returning" << endl;
                    return std::move(ret);
                }
                
//                 cudaSparse<T> operator()(const gMSparse<T> & mat,bool transpose,autoCudaHandles & cudaHandles) const{
// //                     cout << "to SpMat..";cout.flush();
//                     cudaSparse<T> ret( (transpose)?(mat.i_ncols):(mat.i_nrows),(transpose)?(mat.i_nrows):(mat.i_ncols),mat.i_nnz);
//                     if(!transpose){                     
//                         cudaSparse<T> csrTrans(mat.i_ncols,mat.i_nrows,mat.i_nnz);
//                         csrTrans.i_csrValues.copyFrom((const T*) mat.i_cscValues,mat.i_nnz);
//                         csrTrans.i_csrRowPtr.copyFrom((const indexType *) mat.i_cscColPtr,mat.i_ncols+1);
//                         csrTrans.i_csrColInd.copyFrom((const indexType *) mat.i_cscRowIndex,mat.i_nnz);
//                      
//                         size_t bufferSize;
//                         CHECK_CUSPARSE_ERROR( cusparseCsr2cscEx2_bufferSize(
//                             cudaHandles,
//                             ret.i_ncols,
//                             ret.i_nrows,
//                             ret.i_nnz,
//                             csrTrans.i_csrValues,
//                             (const int *) (const cudaIndexType *) csrTrans.i_csrRowPtr,
//                             (const int *) (const cudaIndexType *) csrTrans.i_csrColInd,
//                             ret.i_csrValues,
//                             (int *) (cudaIndexType *) ret.i_csrRowPtr,
//                             (int *) (cudaIndexType *) ret.i_csrColInd,
//                             (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
//                             CUSPARSE_ACTION_NUMERIC,
//                             CUSPARSE_INDEX_BASE_ZERO,
//                             CUSPARSE_CSR2CSC_ALG1,
//                             &bufferSize
//                         ));
//                         CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                         
//                         cudaMemory<unsigned char> buffer(bufferSize);
//                     
//                         CHECK_CUSPARSE_ERROR(cusparseCsr2cscEx2(
//                             cudaHandles,
//                             ret.i_ncols,
//                             ret.i_nrows,
//                             ret.i_nnz,
//                             csrTrans.i_csrValues,
//                             (const int *) (const cudaIndexType *) csrTrans.i_csrRowPtr,
//                             (const int *) (const cudaIndexType *) csrTrans.i_csrColInd,
//                             ret.i_csrValues,
//                             (int *) (cudaIndexType *) ret.i_csrRowPtr,
//                             (int *) (cudaIndexType *) ret.i_csrColInd,
//                             (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
//                             CUSPARSE_ACTION_NUMERIC,
//                             CUSPARSE_INDEX_BASE_ZERO,
//                             CUSPARSE_CSR2CSC_ALG1,
//                             (void *) buffer
//                         ));
//                         CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                     }else{
//                         ret.i_csrValues.copyFrom((const T*) mat.i_cscValues,mat.i_nnz);
//                         ret.i_csrRowPtr.copyFrom((const indexType *) mat.i_cscColPtr,mat.i_ncols+1);
//                         ret.i_csrColInd.copyFrom((const indexType *) mat.i_cscRowIndex,mat.i_nnz);
//                     }
// //                     cout << "returning" << endl;
//                     return std::move(ret);
//                 }
                
            };
            
            template<template<typename> class A,template<typename> class B,typename T>
            B<T> convert(A<T> mat,bool transpose,autoCudaHandles & cudaHandles){
//                 cout << "About to call convert_method" << endl;
                convert_method<A,B,T> method;
                return method(std::move(mat),transpose,cudaHandles);
            }
//END cuda conversion methods
            
//BEGIN cuda methods
            template<typename T>
            cudaDense<T> & multiply_cuda(const cudaDense<T> & A,const cudaDense<T> & B,cudaDense<T> & C,bool transposeA,bool transposeB,T coefA, T coefC, autoCudaHandles & cudaHandles){
                CHECK_CUBLAS_ERROR(cublasTgemm(
                    cudaHandles,
                    (transposeA)?(CUBLAS_OP_T):(CUBLAS_OP_N), 
                    (transposeB)?(CUBLAS_OP_T):(CUBLAS_OP_N),
                    C.i_nrows, 
                    C.i_ncols, 
                    (transposeA)?(A.i_nrows):(A.i_ncols),
                    &coefA,
                    (const T *) A.i_values, 
                    A.i_nrows,
                    (const T *)B.i_values, 
                    B.i_nrows,
                    &coefC,
                    (T*)C.i_values, 
                    C.i_nrows            
                ));
                return C;
            }
            
            template<typename T>
            cudaDense<T> & multiply_cuda(const cudaSparse<T> & A,const cudaDense<T> & B, cudaDense<T> & C,bool transposeA,bool transposeB,T coefA,T coefC, autoCudaHandles & cudaHandles){
            // 	auto prime_time = std::chrono::high_resolution_clock::now();
                size_t bufferSize;
                
                CHECK_CUSPARSE_ERROR(cusparseSpMM_bufferSize(
                    cudaHandles,
                    (transposeA)?(CUSPARSE_OPERATION_TRANSPOSE):(CUSPARSE_OPERATION_NON_TRANSPOSE),
                    (transposeB)?(CUSPARSE_OPERATION_TRANSPOSE):(CUSPARSE_OPERATION_NON_TRANSPOSE),
                    (void*)&coefA,
                    A.i_descr,
                    B.i_descr,
                    (void*)&coefC,
                    C.i_descr,
                    (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                    CUSPARSE_SPMM_CSR_ALG1,
                    &bufferSize
                ));
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                    
//                 void * buffer;
//                 if (bufferSize > 0) {
//                         CHECK_CUDART_ERROR(cudaMalloc(&buffer, bufferSize));
//                 }
                cudaMemory<unsigned char> buffer(bufferSize);
                
                CHECK_CUSPARSE_ERROR(cusparseSpMM(
                    cudaHandles,
                    (transposeA)?(CUSPARSE_OPERATION_TRANSPOSE):(CUSPARSE_OPERATION_NON_TRANSPOSE),
                    (transposeB)?(CUSPARSE_OPERATION_TRANSPOSE):(CUSPARSE_OPERATION_NON_TRANSPOSE),
                    (void*)&coefA,
                    A.i_descr,
                    B.i_descr,
                    (void*)&coefC,
                    C.i_descr,
                    (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F),
                    CUSPARSE_SPMM_CSR_ALG1,
                    (void *) buffer
                ));
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                
                return C;
            }

//BEGIN old sparse sparse mult
//             template<typename T>
//             cudaSparse<T> & multiply_cuda(const cudaSparse<T> & A,const cudaSparse<T> & B,cudaSparse<T> & C,T coefA, autoCudaHandles & cudaHandles){
//
//                 int baseC, nnzC;
//                 int *nnzTotalDevHostPtr = &nnzC;
//                 T alpha = coefA;
//                 T beta  =  0;
//                 size_t bufferSize;
//
//                 class optinfo{
//                 public:
//                     csrgemm2Info_t info = NULL;
//                     cusparseMatDescr_t dA,dB,dC,dD;
//
//                     optinfo(){
//                         CHECK_CUSPARSE_ERROR(cusparseCreateMatDescr(&dA));
//                         CHECK_CUSPARSE_ERROR(cusparseCreateMatDescr(&dB));
//                         CHECK_CUSPARSE_ERROR(cusparseCreateMatDescr(&dC));
//                         CHECK_CUSPARSE_ERROR(cusparseCreateMatDescr(&dD));
//                         CHECK_CUSPARSE_ERROR(cusparseCreateCsrgemm2Info(&info));
//                         CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                     }
//
//                     ~optinfo(){
//                         CHECK_CUSPARSE_ERROR(cusparseDestroyCsrgemm2Info(info));
//                         CHECK_CUSPARSE_ERROR(cusparseDestroyMatDescr(dA));
//                         CHECK_CUSPARSE_ERROR(cusparseDestroyMatDescr(dB));
//                         CHECK_CUSPARSE_ERROR(cusparseDestroyMatDescr(dC));
//                         CHECK_CUSPARSE_ERROR(cusparseDestroyMatDescr(dD));
//                         CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                     }
//                 } opt;
//
//                 cudaSparse<T> D(A.i_nrows,B.i_ncols,0);
//
//                 CHECK_CUSPARSE_ERROR(cusparseTcsrgemm2_bufferSizeExt(
//                     cudaHandles,
//                     A.i_nrows,
//                     B.i_ncols,
//                     A.i_ncols,
//                     &alpha,
//                     opt.dA,
//                     A.i_nnz,
//                     (const int *) (const cudaIndexType *) A.i_csrRowPtr,
//                                                                      (const int *)  (const cudaIndexType *) A.i_csrColInd,
//                                                                      opt.dB,
//                                                                      B.i_nnz,
//                                                                      (const int *)  (const cudaIndexType *) B.i_csrRowPtr,
//                                                                      (const int *)  (const cudaIndexType *) B.i_csrColInd,
//                                                                      &beta,
//                                                                      opt.dD,
//                                                                      0,
//                                                                      (const int *)  (const cudaIndexType *) D.i_csrRowPtr,
//                                                                      (const int *)  (const cudaIndexType *) D.i_csrColInd,
//                                                                      opt.info,
//                                                                      &bufferSize
//                 ));
//                 CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                 cudaMemory<unsigned char> buffer(bufferSize);
//
//                 CHECK_CUSPARSE_ERROR(cusparseXcsrgemm2Nnz(
//                     cudaHandles,
//                     A.i_nrows,
//                     B.i_ncols,
//                     A.i_ncols,
//                     opt.dA,
//                     A.i_nnz,
//                     (const int *) (const cudaIndexType *) A.i_csrRowPtr,
//                                                           (const int *) (const cudaIndexType *) A.i_csrColInd,
//                                                           opt.dB,
//                                                           B.i_nnz,
//                                                           (const int *) (const cudaIndexType *)B.i_csrRowPtr,
//                                                           (const int *) (const cudaIndexType *)B.i_csrColInd,
//                                                           opt.dD,
//                                                           0,
//                                                           (const int *) (const cudaIndexType *)D.i_csrRowPtr,
//                                                           (const int *) (const cudaIndexType *)D.i_csrColInd,
//                                                           opt.dC,
//                                                           (int *) (cudaIndexType *)C.i_csrRowPtr,
//                                                           nnzTotalDevHostPtr,
//                                                           opt.info,
//                                                           (void *) buffer
//                 ));
//                 CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//
//                 if (NULL != nnzTotalDevHostPtr){
//                     C.i_nnz = *nnzTotalDevHostPtr;
//                 }else{
//                     cudaMemcpy(&C.i_nnz, C.i_csrRowPtr+A.i_nrows, sizeof(cudaIndexType), cudaMemcpyDeviceToHost);
//                     cudaMemcpy(&baseC, C.i_csrRowPtr, sizeof(cudaIndexType), cudaMemcpyDeviceToHost);
//                     C.i_nnz -= baseC;
//                 }
//
//                 C.i_csrColInd=cudaMemory<cudaIndexType>(C.i_nnz);
//                 C.i_csrValues=cudaMemory<T>(C.i_nnz);
//
//                 CHECK_CUSPARSE_ERROR(cusparseTcsrgemm2(
//                     cudaHandles,
//                     A.i_nrows,
//                     B.i_ncols,
//                     A.i_ncols,
//                     &alpha,
//                     opt.dA,
//                     A.i_nnz,
//                     (const T *) A.i_csrValues,
//                                                        (const int *) (const cudaIndexType *)A.i_csrRowPtr,
//                                                        (const int *) (const cudaIndexType *)A.i_csrColInd,
//                                                        opt.dB,
//                                                        B.i_nnz,
//                                                        (const T *) B.i_csrValues,
//                                                        (const int *) (const cudaIndexType *)B.i_csrRowPtr,
//                                                        (const int *) (const cudaIndexType *)B.i_csrColInd,
//                                                        &beta,
//                                                        opt.dD,
//                                                        0,
//                                                        (const T *) D.i_csrValues,
//                                                        (const int *) (const cudaIndexType *)D.i_csrRowPtr,
//                                                        (const int *) (const cudaIndexType *)D.i_csrColInd,
//                                                        opt.dC,
//                                                        (T *) C.i_csrValues,
//                                                        (int *) (cudaIndexType *)C.i_csrRowPtr,
//                                                        (int *) (cudaIndexType *)C.i_csrColInd,
//                                                        opt.info,
//                                                        (void *) buffer
//                 ));
//                 CHECK_CUDART_ERROR(cudaDeviceSynchronize());
//                 return C;
//             }
//END old sparse sparse mult

            template<typename T>
            cudaSparse<T> & multiply_cuda_default(const cudaSparse<T> & A,const cudaSparse<T> & B,cudaSparse<T> & C,T coefA, autoCudaHandles & cudaHandles){
                // std::cout << "CALLED SPARSE SPARSE MULTIPLICATION!" << std::endl;

                T               alpha       = coefA;
                T               beta        = 0.0;
                cusparseOperation_t opA         = CUSPARSE_OPERATION_NON_TRANSPOSE;
                cusparseOperation_t opB         = CUSPARSE_OPERATION_NON_TRANSPOSE;
                cudaDataType        computeType = (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F);


                size_t bufferSize1 = 0,    bufferSize2 = 0;

                C=cudaSparse<T>(A.i_nrows,B.i_ncols,0);

                cusparseSpGEMMDescr_t spgemmDesc;
                CHECK_CUSPARSE_ERROR( cusparseSpGEMM_createDescr(&spgemmDesc) );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                // std::cout << "Created descr" << std::endl;

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_workEstimation(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType, CUSPARSE_SPGEMM_DEFAULT,
                    spgemmDesc,
                    &bufferSize1,
                    NULL)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                cudaMemory<unsigned char> dBuffer1(bufferSize1);
                // std::cout << "bufferSize1 ok" << std::endl;

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_workEstimation(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType, CUSPARSE_SPGEMM_DEFAULT,
                    spgemmDesc,
                    &bufferSize1,
                    (void *)dBuffer1)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                // std::cout << "buffer1 ok" << std::endl;

                // ask bufferSize2 bytes for external memory
                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_compute(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType, CUSPARSE_SPGEMM_DEFAULT,
                    spgemmDesc,
                    &bufferSize2,
                    NULL)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());

                // double Gb=(1024*1024*1024);
                // std::cout << bufferSize1 << "\t" << bufferSize2 << std::endl;
                // std::cout << (double) bufferSize1/Gb << "\t" << (double)bufferSize2/Gb << std::endl;

                cudaMemory<unsigned char> dBuffer2(bufferSize2);
                // std::cout << "bufferSize2 ok" << std::endl;

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_compute(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType, CUSPARSE_SPGEMM_DEFAULT,
                    spgemmDesc,
                    &bufferSize2,
                    dBuffer2)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                // std::cout << "compute ok" << std::endl;

                int64_t C_nrows, C_ncols, C_nnz;
                CHECK_CUSPARSE_ERROR( cusparseSpMatGetSize(C.i_descr, &C_nrows, &C_ncols,&C_nnz));
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());

                // std::cout << C_nrows << "\t" << C_ncols << "\t" << C_nnz << std::endl;
                C=cudaSparse<T> (C_nrows,C_ncols,C_nnz);

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_copy(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType,
                    CUSPARSE_SPGEMM_DEFAULT,
                    spgemmDesc)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                // std::cout << "results copied" << std::endl;

                CHECK_CUSPARSE_ERROR( cusparseSpGEMM_destroyDescr(spgemmDesc) );
                return C;
            }

            template<typename T>
            cudaSparse<T> & multiply_cuda(const cudaSparse<T> & A,const cudaSparse<T> & B,cudaSparse<T> & C,T coefA, autoCudaHandles & cudaHandles){
                std::cout << "CALLED SPARSE SPARSE MULTIPLICATION!" << std::endl;

                T alpha = coefA;
                T beta = 0.0;
                cusparseOperation_t opA = CUSPARSE_OPERATION_NON_TRANSPOSE;
                cusparseOperation_t opB = CUSPARSE_OPERATION_NON_TRANSPOSE;
                cusparseSpGEMMAlg_t  alg = CUSPARSE_SPGEMM_ALG3;
                cudaDataType        computeType = (typeid(T)==typeid(float))?(CUDA_R_32F):(CUDA_R_64F);

                int64_t              num_prods;
                T                    chunk_fraction = 0.2;
                size_t bufferSize1 = 0;
                size_t bufferSize2 = 0;
                size_t bufferSize3 = 0;
                double Gb = (1024*1024*1024);


                C=cudaSparse<T>(A.i_nrows,B.i_ncols,0);

                cusparseSpGEMMDescr_t spgemmDesc;
                CHECK_CUSPARSE_ERROR( cusparseSpGEMM_createDescr(&spgemmDesc) );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                std::cout << "Created descr" << std::endl;

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_workEstimation(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType,
                    alg,
                    spgemmDesc,
                    &bufferSize1,
                    NULL)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                std::cout << "bufferSize1 " << bufferSize1 << "\t(" << (double)bufferSize1/Gb << ")" << std::endl;
                cudaMemory<unsigned char> dBuffer1(bufferSize1);

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_workEstimation(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType,
                    alg,
                    spgemmDesc,
                    &bufferSize1,
                    (void *)dBuffer1)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                std::cout << "buffer1 ok" << std::endl;

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_getNumProducts(
                    spgemmDesc,
                    &num_prods)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_estimateMemory(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType,
                    alg,
                    spgemmDesc,
                    chunk_fraction,
                    &bufferSize3,
                    NULL,
                    NULL)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                std::cout << "bufferSize3 " << bufferSize3 << "\t(" << (double)bufferSize3/Gb << ")" << std::endl;

                // {
                    //cudaMemory<unsigned char> dBuffer3(bufferSize3);
                    // std::cout << "buffer3 ok" << std::endl;
                    //CHECK_CUDA( cudaFree(dBuffer3) ) // dBuffer3 can be safely freed to
                    // save more memory

                    CHECK_CUSPARSE_ERROR(cusparseSpGEMM_estimateMemory(
                        cudaHandles,
                        opA,
                        opB,
                        &alpha,
                        A.i_descr,
                        B.i_descr,
                        &beta,
                        C.i_descr,
                        computeType,
                        alg,
                        spgemmDesc,
                        chunk_fraction,
                        &bufferSize3,
                        (void *)cudaMemory<unsigned char>(bufferSize3),
                        &bufferSize2)
                    );
                    CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                // }


                std::cout << "bufferSize2 " << bufferSize2 << "\t(" << (double)bufferSize2/Gb << ")" << std::endl;
                //dBuffer2=cudaMemory<unsigned char>(bufferSize2);
                cudaMemory<unsigned char> dBuffer2(bufferSize2);
                std::cout << "buffer2 ok" << std::endl;



                // compute the intermediate product of A * B
                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_compute(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType,
                    alg,
                    spgemmDesc,
                    &bufferSize2,
                    (void *)dBuffer2)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());

                int64_t C_nrows, C_ncols, C_nnz;
                CHECK_CUSPARSE_ERROR( cusparseSpMatGetSize(C.i_descr, &C_nrows, &C_ncols,&C_nnz));
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());

                // std::cout << C_nrows << "\t" << C_ncols << "\t" << C_nnz << std::endl;
                C=cudaSparse<T> (C_nrows,C_ncols,C_nnz);

                CHECK_CUSPARSE_ERROR(cusparseSpGEMM_copy(
                    cudaHandles,
                    opA,
                    opB,
                    &alpha,
                    A.i_descr,
                    B.i_descr,
                    &beta,
                    C.i_descr,
                    computeType,
                    CUSPARSE_SPGEMM_DEFAULT,
                    spgemmDesc)
                );
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                // std::cout << "results copied" << std::endl;

                CHECK_CUSPARSE_ERROR( cusparseSpGEMM_destroyDescr(spgemmDesc) );
                return C;
            }

            
            template<typename T>
            cudaDense<T> & elementwise_multiply(const cudaDense<T> & A,const cudaDense<T> & B, cudaDense<T> & R);
            
            template<typename T>
            cudaDense<T> & elementwise_divide(const cudaDense<T> & A,const cudaDense<T> & B, cudaDense<T> & R);
            
            template<typename T>
            cudaDense<T> & elementwise_add(const cudaDense<T> & A,const cudaDense<T> & B, cudaDense<T> & R);
            
            template<typename T>
            cudaDense<T> & elementwise_add_constant(const cudaDense<T> & A,T const, cudaDense<T> & R);
            
            template<typename T>
            cudaDense<T> & elementwise_multiply_constant(const cudaDense<T> & A,T const, cudaDense<T> & R);
            
            template<typename T>
            T elementwise_max(const cudaDense<T> & A);
            
            template<typename T>
            T elementwise_min(const cudaDense<T> & A);
            
            template<typename T>
            std::vector<cudaIndexType> find(const cudaDense<T> & A, T threshold,bool greater);
            
            template<typename T>
            cudaDense<T> & select_columns_inplace(cudaDense<T> & A,const std::vector<cudaIndexType> & columns);
            
            template<template<typename> class M, typename T>
            M<T> select_rows(const M<T> & A, const std::vector<cudaIndexType> & rows);
            
            //ATTENTION: THIS FUNCTIOn RETURNS A lower (or upper) VERSION OF THE SYMMETRIC RESULT op(A)A';
            template<typename T>
            cudaDense<T> & syrk(const cudaDense<T> & A, cudaDense<T> & R,bool transpose,bool lower,T alpha, T beta, autoCudaHandles & cudaHandles){
                CHECK_CUBLAS_ERROR(cublasTsyrk(
                    cudaHandles,
                    lower?CUBLAS_FILL_MODE_LOWER:CUBLAS_FILL_MODE_UPPER, 
                    transpose?CUBLAS_OP_T:CUBLAS_OP_N,
                    transpose?A.i_ncols:A.i_nrows, 
                    transpose?A.i_nrows:A.i_ncols,
                    &alpha,
                    (const T *) A.i_values,
                    A.i_nrows,
                    &beta,
                    (T *) R.i_values, 
                    R.i_nrows
                ));
                cudaDeviceSynchronize();
                return R;
            }
            
            template<typename T> 
            cudaDense<T> & syevd(cudaDense<T> & A, cudaDense<T> & R,bool lower, autoCudaHandles & cudaHandles){
                cudaIndexType lwork=0;
                
                CHECK_CUSOLVER_ERROR(cusolverDnTsyevd_bufferSize(
                    cudaHandles,
                    CUSOLVER_EIG_MODE_NOVECTOR,
                    lower?CUBLAS_FILL_MODE_LOWER:CUBLAS_FILL_MODE_UPPER,
                    A.i_ncols,
                    (const T *) nullptr,
                    A.i_ncols,
                    (const T *) nullptr,
                    (int *) &lwork
                ),0);
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                
                cudaMemory<T> work(lwork);
                cudaMemory<int> devInfo(1);
                
                CHECK_CUDART_ERROR(cudaMalloc((void **)&devInfo, sizeof(int)));
                CHECK_CUSOLVER_ERROR(
                cusolverDnTsyevd(
                    cudaHandles,
                    CUSOLVER_EIG_MODE_NOVECTOR,
                    CUBLAS_FILL_MODE_LOWER,
                    A.i_ncols,
                    (T *) A.i_values,
                    A.i_ncols,
                    (T *) R.i_values,
                    (T *) work,
                    lwork,
                    (int *) devInfo
                ),0);
                CHECK_CUDART_ERROR(cudaDeviceSynchronize());
                
                return R;
            }
            
            template<typename T>
            cudaDense<T> & normalise_inplace(cudaDense<T> & A);
            
            template<typename T>
            cudaDense<T> & get_diag(const cudaDense<T> & M,cudaDense<T> & diag);
            
            template<typename T>
            cudaDense<T> & sqrt(const cudaDense<T> & M,cudaDense<T> & R);
            
            
            
// these seem not to be better than using multiply on the same matrix            
//             template<typename T>
//             cudaDense<T> & multiply_transpose_cuda(const cudaDense<T> & M,cudaDense<T> & C,bool transpose,autoCublasHandle & cublasHandle){
//                 T alpha=1;
//                 T beta=0;
//                 indexType n=(transpose)?(M.i_ncols):(M.i_nrows);
//                 indexType k=(transpose)?(M.i_nrows):(M.i_ncols);
//                 
//                 CHECK_CUBLAS_ERROR(cublasTsyrk(
//                     cublasHandle,
//                     CUBLAS_FILL_MODE_UPPER, 
//                     (transpose)?(CUBLAS_OP_T):(CUBLAS_OP_N),
//                     n,
//                     k,
//                     &alpha,
//                     (const T *) M.i_values, 
//                     (transpose)?(k):(n),
//                     &beta,
//                     (T *) C.i_values, 
//                     n
//                 ));
//                 return C;
//             }
//             
//             template<typename T>
//             cudaSparse<T> multiply_transpose_cuda(const cudaSparse<T> & M,cudaSparse<T> & C,bool transpose){
//             }

            
//END cuda methods
        }
	}
}
#pragma GCC diagnostic pop
