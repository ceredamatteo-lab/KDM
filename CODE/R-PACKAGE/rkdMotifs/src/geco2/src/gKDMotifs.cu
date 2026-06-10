#include <cuda_device_runtime_api.h>
#include <thrust/device_ptr.h>
#include <thrust/transform.h>

#include "gMethods_cuda.h"
using namespace geco::methods;

namespace geco{
    namespace methods{
        namespace cuda{
            template <typename T>
            cudaDense<T> & isGreater(const cudaDense<T> & A, cudaDense<T> & B, cudaDense<T> & C,indexType nrep);
        }
    }
};


template <typename T>
struct igfunctor//: public thrust::binary_function<T,T,T>
{
    T i_nrep;
    igfunctor(T nrep):i_nrep(nrep){
    }
    
  __host__ __device__
  T operator()(const T & x,const T y) const
  {
    return ((x>0)?(y>x):(y<x));
  }
};




template <typename T>
geco::methods::cuda::cudaDense<T> & geco::methods::cuda::isGreater(const cuda::cudaDense<T> & A, cuda::cudaDense<T> & B, cuda::cudaDense<T> & C,geco::methods::indexType nrep){
    cuda::cudaIndexType nelm=A.i_nrows*A.i_ncols;
    const thrust::device_ptr<const T> Ap=thrust::device_pointer_cast<const T>(A.i_values);
    thrust::device_ptr<T> Bp=thrust::device_pointer_cast<T>(B.i_values);
    thrust::device_ptr<T> Cp=thrust::device_pointer_cast<T>(C.i_values);
    thrust::transform(Ap,Ap+nelm,Bp,Bp,igfunctor<T>(nrep));
    thrust::transform(Cp,Cp+nelm,Bp,Cp,thrust::plus<T>());
    return C;
}

template geco::methods::cuda::cudaDense<float> & geco::methods::cuda::isGreater(const cuda::cudaDense<float> & A, cuda::cudaDense<float> & B, cuda::cudaDense<float> & C,geco::methods::indexType nrep);
template geco::methods::cuda::cudaDense<double> & geco::methods::cuda::isGreater(const cuda::cudaDense<double> & A, cuda::cudaDense<double> & B, cuda::cudaDense<double> & C,geco::methods::indexType nrep);
