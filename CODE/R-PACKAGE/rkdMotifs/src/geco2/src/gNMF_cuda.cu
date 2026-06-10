#include "gNMF_cuda.h"
#include <cuda_device_runtime_api.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/permutation_iterator.h>



#include <thrust/inner_product.h>
#include <thrust/generate.h>
// #include <thrust/reduce.h>
// #include "thrust/for_each.h"

using namespace geco::methods::cuda;

template<typename T>
struct binf //: thrust::binary_function<T,const T &,const T &>
{
public:
__host__ __device__
    T operator()(const T & w, const T & oldw){
        return (oldw-w)*(oldw-w);
    }
};

template<typename T>
T geco::methods::cuda::nmf_change(const cudaDense<T> & W,const cudaDense<T> & oldW){
    
    const thrust::device_ptr<const T> Wp=thrust::device_pointer_cast<const T>(W.i_values);
    const thrust::device_ptr<const T> oldWp=thrust::device_pointer_cast<const T>(oldW.i_values);

    cudaIndexType size=W.i_nrows*W.i_ncols;
	T res1=::sqrt(thrust::inner_product(Wp,Wp+size,oldWp,(T) 0.0,thrust::plus<T>(),binf<T>()));
	T res2=::sqrt(thrust::inner_product(oldWp,oldWp+size,oldWp,(T)0.0,thrust::plus<T>(),thrust::multiplies<T>()));
	return pow((res1)/(res2),2);
}


struct calcDiagIndex//: public thrust::unary_function<cudaIndexType,cudaIndexType>
{
    cudaIndexType i_nrows_plus_one;
    calcDiagIndex(cudaIndexType nrows):i_nrows_plus_one(nrows+1){
    }
    
  __host__ __device__
  cudaIndexType operator()(cudaIndexType x) const
  {
    return x *( i_nrows_plus_one);  //this is x * i_nrows +x
  }
};

template<typename T>
cudaDense<T> & geco::methods::cuda::get_sigma(const cudaDense<T> & M,cudaDense<T> & sigma){
    const thrust::device_ptr<const T> Mp=thrust::device_pointer_cast<const T>(M.i_values);
    thrust::device_ptr<T> sigmap=thrust::device_pointer_cast<T>(sigma.i_values);

    auto Mbegin=make_permutation_iterator(Mp,thrust::make_transform_iterator(thrust::make_counting_iterator(0),calcDiagIndex(M.i_nrows)));
    thrust::copy(Mbegin,Mbegin+sigma.i_ncols,sigmap);

    return sigma;
}



template<typename T>
struct invsigma//: public thrust::unary_function<T,T>
{
  __host__ __device__
  T operator()(T x) const
  {
    return 1/x;
  }
};

template <typename T>
cudaDense<T> & geco::methods::cuda::set_dgma(const cudaDense<T> & sigma, cudaDense<T> & dsigma){
    CHECK_CUDART_ERROR(cudaMemset((T *) dsigma.i_values,0,dsigma.i_nrows * dsigma.i_ncols * sizeof(T)));
    
    thrust::device_ptr<T> dsigmap=thrust::device_pointer_cast<T>(dsigma.i_values);
    const thrust::device_ptr<const T> sigmap=thrust::device_pointer_cast<const T>(sigma.i_values);
    
    auto dsigmabegin=thrust::make_permutation_iterator(
        dsigmap,
        thrust::make_transform_iterator(thrust::make_counting_iterator(0),calcDiagIndex(dsigma.i_nrows))
    );
    thrust::transform(sigmap,sigmap+sigma.i_ncols,dsigmabegin,invsigma<T>());
    return dsigma;
}


template<typename T>
struct norsigma//: public thrust::unary_function<T,T>
{
  __host__ __device__
  T operator()(T x) const
  {
    return x;
  }
};

template <typename T>
cudaDense<T> & geco::methods::cuda::set_dgma2(const cudaDense<T> & sigma, cudaDense<T> & dsigma){
    CHECK_CUDART_ERROR(cudaMemset((T *) dsigma.i_values,0,dsigma.i_nrows * dsigma.i_ncols * sizeof(T)));
    
    thrust::device_ptr<T> dsigmap=thrust::device_pointer_cast<T>(dsigma.i_values);
    const thrust::device_ptr<const T> sigmap=thrust::device_pointer_cast<const T>(sigma.i_values);
    
    auto dsigmabegin=thrust::make_permutation_iterator(
        dsigmap,
        thrust::make_transform_iterator(thrust::make_counting_iterator(0),calcDiagIndex(dsigma.i_nrows))
    );
    thrust::copy(sigmap,sigmap+sigma.i_ncols,dsigmabegin);
    return dsigma;
}


template <typename T>
void geco::methods::cuda::printout(const cudaDense<T> & M,const autoCudaHandles & cudaHandles){
    const thrust::device_ptr<const T> Mp=thrust::device_pointer_cast<const T>(M.i_values);
    
    for(cudaIndexType r=0;r<M.i_nrows;r++){
        for(cudaIndexType c=0;c<M.i_ncols;c++){
            std::cout << Mp[c *M.i_nrows+r] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


template float geco::methods::cuda::nmf_change<float>(const cudaDense<float> & W,const cudaDense<float> & oldW);
template double geco::methods::cuda::nmf_change<double>(const cudaDense<double> & W,const cudaDense<double> & oldW);

template cudaDense<float> & geco::methods::cuda::get_sigma(const cudaDense<float> & M,cudaDense<float> & sigma);
template cudaDense<double> & geco::methods::cuda::get_sigma(const cudaDense<double> & M,cudaDense<double> & sigma);

template cudaDense<float> & geco::methods::cuda::set_dgma(const cudaDense<float> & sigma, cudaDense<float> & dsigma);
template cudaDense<double> & geco::methods::cuda::set_dgma(const cudaDense<double> & sigma, cudaDense<double> & dsigma);

template cudaDense<float> & geco::methods::cuda::set_dgma2(const cudaDense<float> & sigma, cudaDense<float> & dsigma);
template cudaDense<double> & geco::methods::cuda::set_dgma2(const cudaDense<double> & sigma, cudaDense<double> & dsigma);


template void geco::methods::cuda::printout(const cudaDense<float> & M,const autoCudaHandles & cudaHandles);
template void geco::methods::cuda::printout(const cudaDense<double> & M,const autoCudaHandles & cudaHandles);


