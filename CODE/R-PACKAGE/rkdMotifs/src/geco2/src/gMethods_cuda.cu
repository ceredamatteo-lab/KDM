#include <cuda_device_runtime_api.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/transform.h>
#include <thrust/find.h>
#include <thrust/copy.h>
#include <thrust/gather.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>



#include<list>

//#include <thrust/device_vector.h>
// #include <thrust/inner_product.h>
// #include <thrust/reduce.h>
// #include "thrust/for_each.h"
#include "gMethods_cuda.h"

//using namespace geco::methods;
//using namespace geco::methods::cuda;

template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::elementwise_multiply(const cudaDense<T> & A,const cudaDense<T> & B, cudaDense<T> & R){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    const thrust::device_ptr<const T> Bp=thrust::device_pointer_cast<const T>(B.i_values);
    thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);

    thrust::transform(Ap,Ap+A.i_nrows*A.i_ncols,Bp,Rp,thrust::multiplies<T>());
    return R;
}

template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::elementwise_divide(const cudaDense<T> & A,const cudaDense<T> & B, cudaDense<T> & R){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    const thrust::device_ptr<const T> Bp=thrust::device_pointer_cast<const T>(B.i_values);
    thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);
    
    thrust::transform(Ap,Ap+A.i_nrows*A.i_ncols,Bp,Rp,thrust::divides<T>());
    return R;
}

template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::elementwise_add(const cudaDense<T> & A,const cudaDense<T> & B, cudaDense<T> & R){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    const thrust::device_ptr<const T> Bp=thrust::device_pointer_cast<const T>(B.i_values);
    thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);

    thrust::transform(Ap,Ap+A.i_nrows*A.i_ncols,Bp,Rp,thrust::plus<T>());
    return R;
}




template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::elementwise_add_constant(const cudaDense<T> & A,T constant, cudaDense<T> & R){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);
    thrust::transform(Ap,Ap+A.i_nrows*A.i_ncols,thrust::make_constant_iterator(constant),Rp,thrust::plus<T>());
    return R;
}

template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::elementwise_multiply_constant(const cudaDense<T> & A,T const constant, cudaDense<T> & R){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);
    thrust::transform(Ap,Ap+A.i_nrows*A.i_ncols,thrust::make_constant_iterator(constant),Rp,thrust::multiplies<T>());
    return R;
}

template<typename T>
T geco::methods::cuda::elementwise_max(const cudaDense<T> & A){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    return *thrust::max_element(Ap,Ap+A.i_nrows*A.i_ncols);
}
            
template<typename T>
T geco::methods::cuda::elementwise_min(const cudaDense<T> & A){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    return *(thrust::min_element(Ap,Ap+A.i_nrows*A.i_ncols));
}

template<typename T>
struct gt{
public:
    T i_threshold;
    gt(T threshold):i_threshold(threshold){}
    __host__ __device__
    bool operator()(T x){
        return x > i_threshold;
    }
};

template<typename T>
struct lt{
public:
    T i_threshold;
    lt(T threshold):i_threshold(threshold){}
    __host__ __device__
    bool operator()(T x){
        return x < i_threshold;
    }
};

template<typename T>
std::vector<geco::methods::cuda::cudaIndexType> geco::methods::cuda::find(const cudaDense<T> & A, T threshold,bool greater){
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    
    std::list<cudaIndexType> rr;
    if(greater){
        auto i=thrust::find_if(Ap,Ap+A.i_nrows*A.i_ncols,gt<T>(threshold));
        while(i!=Ap+A.i_nrows*A.i_ncols){
            rr.push_back(thrust::distance(Ap,i));
            i=thrust::find_if(i+1,Ap+A.i_nrows*A.i_ncols,gt<T>(threshold));
        }
    }else{
        auto i=thrust::find_if(Ap,Ap+A.i_nrows*A.i_ncols,lt<T>(threshold));
        while(i!=Ap+A.i_nrows*A.i_ncols){
            
            rr.push_back(thrust::distance(Ap,i));
            i=thrust::find_if(i+1,Ap+A.i_nrows*A.i_ncols,lt<T>(threshold));
        }
    }
    std::vector<cudaIndexType> ret(rr.begin(),rr.end());
    return ret;
}

//ALTERNTIVE IMPLIES COPYING TOO MUCH
// template<typename T>
// std::vector<cudaIndexType> geco::methods::cuda::find(const cudaDense<T> & A, T threshold,bool greater){
//     const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
//     thrust::device_vector<T> TT(10);
//     thrust::device_vector<cudaIndexType> res;
//     if(greater){
//         cudaIndexType num=thrust::count_if(Ap,Ap+A.i_nrows*A.i_ncols,thrust::placeholders::_1 > threshold);
//         res=thrust::device_vector<cudaIndexType>(num);
//         //thrust::copy_if(thrust::make_counting_iterator<int>(0),thrust::make_counting_iterator<int>(A.i_nrows*A.i_ncols),Ap,res.begin(),thrust::placeholders::_1 > threshold);
//         thrust::copy_if(thrust::make_counting_iterator<int>(0),thrust::make_counting_iterator<int>(10),TT.begin(),res.begin(),thrust::placeholders::_1 > threshold);
//     }else{
//         cudaIndexType num=thrust::count_if(Ap,Ap+A.i_nrows*A.i_ncols,thrust::placeholders::_1 < threshold);
//         res=thrust::device_vector<cudaIndexType>(num);
//         //thrust::copy_if(thrust::make_counting_iterator<int>(0),thrust::make_counting_iterator<int>(A.i_nrows*A.i_ncols),Ap,res.begin(),thrust::placeholders::_1 < threshold);
//         thrust::copy_if(thrust::make_counting_iterator<int>(0),thrust::make_counting_iterator<int>(10),TT.begin(),res.begin(),thrust::placeholders::_1 < threshold);        
//     }
//     return std::vector<cudaIndexType>(res.begin(),res.end());

// }


template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::select_columns_inplace(cudaDense<T> & A,const std::vector<cudaIndexType> & columns){
    thrust::device_ptr<T> Ap=thrust::device_pointer_cast<T>(A.i_values);
    cudaIndexType newpos=0;
//     cudaDense<T> R(A.i_nrows,columns.size());
//     thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);
    for(size_t i=0;i<columns.size();i++){
        cudaIndexType offset=columns[i]*A.i_nrows;
        thrust::copy(Ap + offset,Ap + offset + A.i_nrows,Ap+newpos);
//         thrust::copy(Ap + offset,Ap + offset + A.i_nrows,Rp+newpos);
        newpos+=A.i_nrows;
    }
    // WE CHANGE THE NUMBER OF COLUMNS BUT WE DO NOT RESIZE THE ARRAY SINCE IT WILL BE FULLY RELEASED REGARDLESS
    A.updateNcols(columns.size());
//     A=R;
    return A;
}


struct inIndex{
    geco::methods::cuda::cudaIndexType * rows;
    geco::methods::cuda::cudaIndexType nrows;
    geco::methods::cuda::cudaIndexType ncols;
    
    inIndex(geco::methods::cuda::cudaIndexType * r,geco::methods::cuda::cudaIndexType nr,geco::methods::cuda::cudaIndexType nc):rows(r),nrows(nr),ncols(nc){
    }
    
  __host__ __device__
  geco::methods::cuda::cudaIndexType operator()(geco::methods::cuda::cudaIndexType x) const{
      return  (x % ncols) * nrows + rows[x/ncols];
  }
};

struct outIndex{
    geco::methods::cuda::cudaIndexType ncols;
    geco::methods::cuda::cudaIndexType nrows;
    
    outIndex(geco::methods::cuda::cudaIndexType nr,geco::methods::cuda::cudaIndexType nc):nrows(nr),ncols(nc){
    }
    
  __host__ __device__
  geco::methods::cuda::cudaIndexType operator()(geco::methods::cuda::cudaIndexType x) const{
      return  (x % ncols) * nrows + x/ncols;
  }
};


template<typename T>
geco::methods::cuda::cudaDense<T> select_rows_dense(const geco::methods::cuda::cudaDense<T> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows){
    const thrust::device_ptr<const T> Av=thrust::device_pointer_cast<const T>(A.i_values);
    
    geco::methods::cuda::cudaDense<T> R(rows.size(),A.i_ncols);
    thrust::device_ptr<T> Rv=thrust::device_pointer_cast<T>(R.i_values);
    
    thrust::device_vector<geco::methods::cuda::cudaIndexType> select(rows.size());
    thrust::copy(rows.begin(),rows.end(),select.begin());
    geco::methods::cuda::cudaIndexType * selp = thrust::raw_pointer_cast(&select[0]);
    
    
    auto inputIterator=thrust::make_permutation_iterator(
        Av,
        thrust::make_transform_iterator(thrust::make_counting_iterator(0),inIndex(selp,A.i_nrows,A.i_ncols))
    );
    
    auto outputIterator=thrust::make_permutation_iterator(
        Rv,
        thrust::make_transform_iterator(thrust::make_counting_iterator(0),outIndex(R.i_nrows,R.i_ncols))
    );
    
    thrust::copy(inputIterator,inputIterator + A.i_ncols * rows.size(), outputIterator);
    
//     thrust::device_vector<cudaIndexType> select(rows.size());
//     thrust::copy(rows.begin(),rows.end(),select.begin());
//     cudaIndexType * selp = thrust::raw_pointer_cast(&select[0]);
//     dim3 blockdim = dim3(128);
//     dim3 griddim = dim3(rows.size());
//     rowcopykernel<<<griddim,blockdim>>>((const T *) A.i_values, (T *) R.i_values, selp, A.i_nrows, A.i_ncols, rows.size());    
    
//     const thrust::device_ptr<const T> Av=thrust::device_pointer_cast<const T>(A.i_values);
//     cudaDense<T> R(rows.size(),A.i_ncols);
//     thrust::device_ptr<T> Rv=thrust::device_pointer_cast<T>(R.i_values);
//     cudaIndexType newpos=0;
//     for(cudaIndexType c=0;c<A.i_ncols;c++){
//         for(cudaIndexType i=0;i<rows.size();i++){
//             Rv[newpos++]=Av[c*A.i_nrows+rows[i]];
//         }
//     }
    
    geco::methods::cuda::cudaIndexType nob2copy=A.i_ncols * rows.size();
    
    
    
    return std::move(R);
}

template<typename T>
geco::methods::cuda::cudaSparse<T> select_rows_sparse(const geco::methods::cuda::cudaSparse<T> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows){
    //uint32_t cuda_select(const gKIndex * columns,size_t ncols,gKType * csrValuesV,uint32_t * csrRowOffsetsV, uint32_t * csrColIndV,	gKType ** csrValuespV, uint32_t ** csrRowOffsetspV, uint32_t ** csrColIndpV){
    const thrust::device_ptr<const T> vV=thrust::device_pointer_cast<const T>(A.i_csrValues);
    const thrust::device_ptr<const geco::methods::cuda::cudaIndexType> roV=thrust::device_pointer_cast<const geco::methods::cuda::cudaIndexType>(A.i_csrRowPtr);
    const thrust::device_ptr<const geco::methods::cuda::cudaIndexType> ciV=thrust::device_pointer_cast<const geco::methods::cuda::cudaIndexType>(A.i_csrColInd);
    //cudaMalloc((void **)csrRowOffsetspV,(ncols+1)*sizeof(uint32_t));            
    
    geco::methods::cuda::cudaMemory<geco::methods::cuda::cudaIndexType> rop(rows.size()+1);
    thrust::device_ptr<geco::methods::cuda::cudaIndexType> ropV=thrust::device_pointer_cast<geco::methods::cuda::cudaIndexType>(rop);
    ropV[0]=0;
    for(size_t i=0;i<rows.size();i++){
        ropV[i+1] = ropV[i] + roV[rows[i]+1]-roV[rows[i]];
    }

    geco::methods::cuda::cudaIndexType nelm=ropV[rows.size()];
    geco::methods::cuda::cudaSparse<T> R(rows.size(),A.i_ncols,nelm);
    //cudaMalloc((void **)csrValuespV,nelm*sizeof(gKType));
    //cudaMalloc((void **)csrColIndpV,nelm*sizeof(uint32_t));
    thrust::device_ptr<T> vpV=thrust::device_pointer_cast<T>(R.i_csrValues);

    thrust::device_ptr<geco::methods::cuda::cudaIndexType> cipV=thrust::device_pointer_cast<geco::methods::cuda::cudaIndexType>(R.i_csrColInd);
    for(size_t i=0;i<rows.size();i++){
        thrust::copy(vV+roV[rows[i]],vV+roV[rows[i]+1],vpV+ropV[i]);
        thrust::copy(ciV+roV[rows[i]],ciV+roV[rows[i]+1],cipV+ropV[i]);
    }
    R.i_csrRowPtr=std::move(rop);
    R.updateDescr();

    return std::move(R);
}

template<typename F1, typename F2, template<typename> class M, typename T>
auto invoke_select_rows(F1 f1, F2 f2, const M<T> & A,const std::vector<geco::methods::cuda::cudaIndexType> & rows) -> decltype(f1(A,rows))
{
    return f1(A,rows);
}

template<typename F1, typename F2, template<typename> class M, typename T>
auto invoke_select_rows(F1 f1, F2 f2, const M<T> & A,const std::vector<geco::methods::cuda::cudaIndexType> & rows) -> decltype(f2(A,rows))
{
    return f2(A,rows);
}

template<template<typename> class M, typename T>
M<T> geco::methods::cuda::select_rows(const M<T> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows){
    return invoke_select_rows(select_rows_dense<T>,select_rows_sparse<T>,A,rows);
}


template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::normalise_inplace(geco::methods::cuda::cudaDense<T> & A){
    thrust::device_ptr<T> Ap=thrust::device_pointer_cast<T>(A.i_values);
    for(geco::methods::cuda::cudaIndexType i=0;i<A.i_ncols;i++){
        T norm=::sqrt(thrust::transform_reduce(Ap+i*A.i_nrows,Ap+(i+1)*A.i_nrows,thrust::square<T>(),(T)0,thrust::plus<T>()));
        thrust::transform(Ap+i*A.i_nrows,Ap+(i+1)*A.i_nrows,thrust::make_constant_iterator(norm),Ap+i*A.i_nrows,thrust::divides<T>());
    }
    return A;
}

struct calcDiagIndex//: public thrust::unary_function<geco::methods::cuda::cudaIndexType,geco::methods::cuda::cudaIndexType>
{
    geco::methods::cuda::cudaIndexType i_nrows_plus_one;
    calcDiagIndex(geco::methods::cuda::cudaIndexType nrows):i_nrows_plus_one(nrows+1){
    }
    
  __host__ __device__
  geco::methods::cuda::cudaIndexType operator()(geco::methods::cuda::cudaIndexType x) const
  {
    return x * ( i_nrows_plus_one);  //this is x * i_nrows + x = x * (i_nrows + 1)
  }
};

template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::get_diag(const geco::methods::cuda::cudaDense<T> & M,geco::methods::cuda::cudaDense<T> & diag){
    const thrust::device_ptr<const T> Mp=thrust::device_pointer_cast<const T>(M.i_values);
    thrust::device_ptr<T> diagp=thrust::device_pointer_cast<T>(diag.i_values);

    auto Mbegin=make_permutation_iterator(Mp,thrust::make_transform_iterator(thrust::make_counting_iterator(0),calcDiagIndex(M.i_nrows)));
    thrust::copy(Mbegin,Mbegin+M.i_nrows,diagp);

    return diag;
    
}

template <typename T>
struct sqrt_functor//: public thrust::unary_function<T,T>
{
  __host__ __device__
  T operator()(const T & x) const
  {
    return sqrt(x);
  }
};

template<typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::sqrt(const geco::methods::cuda::cudaDense<T> & M,geco::methods::cuda::cudaDense<T> & R){
    const thrust::device_ptr<const T> Mp=thrust::device_pointer_cast<const T>(M.i_values);
    thrust::device_ptr<T> Rp=thrust::device_pointer_cast<T>(R.i_values);
    thrust::transform(Mp,Mp + M.i_nrows * M. i_ncols,Rp,sqrt_functor<T>());
    return R;
}





template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::elementwise_multiply(const geco::methods::cuda::cudaDense<float> & A,const geco::methods::cuda::cudaDense<float> & B, cudaDense<float> & R);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::elementwise_multiply(const geco::methods::cuda::cudaDense<double> & A,const geco::methods::cuda::cudaDense<double> & B, geco::methods::cuda::cudaDense<double> & R);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::elementwise_divide(const geco::methods::cuda::cudaDense<float> & A,const geco::methods::cuda::cudaDense<float> & B, geco::methods::cuda::cudaDense<float> & R);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::elementwise_divide(const geco::methods::cuda::cudaDense<double> & A,const geco::methods::cuda::cudaDense<double> & B, geco::methods::cuda::cudaDense<double> & R);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::elementwise_add(const geco::methods::cuda::cudaDense<float> & A,const geco::methods::cuda::cudaDense<float> & B, geco::methods::cuda::cudaDense<float> & R);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::elementwise_add(const geco::methods::cuda::cudaDense<double> & A,const geco::methods::cuda::cudaDense<double> & B, geco::methods::cuda::cudaDense<double> & R);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::elementwise_add_constant(const geco::methods::cuda::cudaDense<float> & A,float constant, geco::methods::cuda::cudaDense<float> & R);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::elementwise_add_constant(const geco::methods::cuda::cudaDense<double> & A,double constant, geco::methods::cuda::cudaDense<double> & R);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::elementwise_multiply_constant(const geco::methods::cuda::cudaDense<float> & A,float constant, geco::methods::cuda::cudaDense<float> & R);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::elementwise_multiply_constant(const geco::methods::cuda::cudaDense<double> & A,double constant, geco::methods::cuda::cudaDense<double> & R);

template float geco::methods::cuda::elementwise_min(const geco::methods::cuda::cudaDense<float> & A);
template double geco::methods::cuda::elementwise_min(const geco::methods::cuda::cudaDense<double> & A);

template float geco::methods::cuda::elementwise_max(const geco::methods::cuda::cudaDense<float> & A);
template double geco::methods::cuda::elementwise_max(const geco::methods::cuda::cudaDense<double> & A);

template std::vector<geco::methods::cuda::cudaIndexType> geco::methods::cuda::find(const geco::methods::cuda::cudaDense<float> & A, float threshold,bool greater);
template std::vector<geco::methods::cuda::cudaIndexType> geco::methods::cuda::find(const geco::methods::cuda::cudaDense<double> & A, double threshold,bool greater);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::select_columns_inplace(geco::methods::cuda::cudaDense<float> & A,const std::vector<geco::methods::cuda::cudaIndexType> & columns);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::select_columns_inplace(geco::methods::cuda::cudaDense<double> & A,const std::vector<geco::methods::cuda::cudaIndexType> & columns);

template geco::methods::cuda::cudaDense<float> geco::methods::cuda::select_rows(const geco::methods::cuda::cudaDense<float> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows);
template geco::methods::cuda::cudaDense<double> geco::methods::cuda::select_rows(const geco::methods::cuda::cudaDense<double> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows);
template geco::methods::cuda::cudaSparse<float> geco::methods::cuda::select_rows(const geco::methods::cuda::cudaSparse<float> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows);
template geco::methods::cuda::cudaSparse<double> geco::methods::cuda::select_rows(const geco::methods::cuda::cudaSparse<double> & A, const std::vector<geco::methods::cuda::cudaIndexType> & rows);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::normalise_inplace(geco::methods::cuda::cudaDense<float> & A);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::normalise_inplace(geco::methods::cuda::cudaDense<double> & A);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::get_diag(const geco::methods::cuda::cudaDense<float> & M,geco::methods::cuda::cudaDense<float> & diag);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::get_diag(const geco::methods::cuda::cudaDense<double> & M,geco::methods::cuda::cudaDense<double> & diag);

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::sqrt(const geco::methods::cuda::cudaDense<float> & M,geco::methods::cuda::cudaDense<float> & R);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::sqrt(const geco::methods::cuda::cudaDense<double> & M,geco::methods::cuda::cudaDense<double> & R);
