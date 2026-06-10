#include "gKmers_cuda.h"
#include <cuda_device_runtime_api.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>


//#include <fstream>

using namespace std;
using namespace geco;
using namespace geco::methods;
using namespace geco::methods::kmers;









struct getLmerGkmer{
    indexType i_nmasks;
    indexType i_gkmer_mask;
    
    getLmerGkmer(indexType nmasks,indexType gkmer_mask):i_nmasks(nmasks),i_gkmer_mask(gkmer_mask){
    }
    
    __host__ __device__
    indexType operator() (const indexType & lmer) const {
        return lmer * i_nmasks + i_gkmer_mask;
    }
};

struct getMismatches{
    indexType i_k;
    indexType i_gkmer_kmer;
    indexType i_nmasks;
    getMismatches(indexType k,indexType gkmer_kmer,indexType nmasks):i_k(k),i_gkmer_kmer(gkmer_kmer),i_nmasks(nmasks){
    }
    
    __host__ __device__
    indexType operator() (const indexType & lmer_kmer) const{
        return i_gkmer_kmer ^ (lmer_kmer/i_nmasks);
    }
};

template <typename T>
struct summer{
    T i_val;
        summer(const T & val):i_val(val){
        }
    
    __host__ __device__
    T operator() (const T val1, const T val2) const{
        return val1 * i_val + val2;
    }
};



template<typename T>
struct cutter{
    T eps=std::numeric_limits<T>::epsilon();
    __host__ __device__
    T operator ()(const T & val){
        return (val>eps)?val:0;
    }
};

template<typename  T>
gCudaCounter<T>::gCudaCounter(const gLmerEstimator<T> & counter):i_counter(counter),i_masks(counter.i_masks),i_mmcoeffs(counter.i_mmcoeffs),i_mlmerGkmer(counter.i_mlmerGkmer),i_lmerDSLmer(counter.i_lmerDSLmer){
}


template<typename T>
bool gCudaCounter<T>::estimatelmers(const std::map<size_t,std::list<unsigned char> > & nseq, std::vector<T> & elmers,indexType base,bool doubleStrand,bool strict) const{
    //vector<list<unsigned char> > nseq = i_counter.seqCode(seq,0,strict);
    bool nnzero=false;
    thrust::device_vector<T> dev_elmers(i_counter.i_nlmers,0.0);
    indexType i_mask=pow(4,i_counter.i_l)-1;
    
    for(auto u=nseq.begin();u!=nseq.end();u++){
        if(u->second.size()>=i_counter.i_l){
            nnzero=true;
            auto bb=u->second.begin();
            indexType lmer=*bb++;
            for(indexType i=1;i<i_counter.i_l-1;i++){
                lmer <<= 2;
                lmer|=*bb++;
            }
            while(bb!=u->second.end()){
                lmer <<= 2;
                lmer|=*bb++;
                indexType mlmer=(lmer & i_mask) * i_counter.i_nmasks;
                for(indexType mask=0;mask<i_counter.i_nmasks;mask++){
                    indexType gkmer_kmer = i_counter.i_mlmerGkmer[mlmer++] / i_counter.i_nmasks;
                    auto Cbegin = thrust::make_permutation_iterator(
                        i_mmcoeffs.begin(),
                        thrust::make_transform_iterator(
                            thrust::make_permutation_iterator(
                                i_mlmerGkmer.begin(),
                                thrust::make_transform_iterator(
                                    thrust::make_counting_iterator(0),
                                    getLmerGkmer(i_counter.i_nmasks,mask)
                                )
                            ),
                            getMismatches(i_counter.i_k,gkmer_kmer,i_counter.i_nmasks)
                        )
                    );
                    thrust::transform(Cbegin,Cbegin+i_counter.i_nlmers,dev_elmers.begin(),dev_elmers.begin(),thrust::plus<T>());
                }
            }
        }
    }
        
    if(doubleStrand){
        thrust::host_vector<T> helmers=dev_elmers;
        thrust::host_vector<T> hdselmers(i_counter.i_nDSlmers,0);
        for(indexType lmer=0;lmer<i_counter.i_nlmers;lmer++){
            hdselmers[i_counter.i_lmerDSLmer[lmer]] += helmers[lmer];
        }
        thrust::transform(hdselmers.begin(),hdselmers.end(),hdselmers.begin(),cutter<T>());
        std::copy(hdselmers.begin(),hdselmers.end(),elmers.begin()+base);
    }else{
        thrust::transform(dev_elmers.begin(),dev_elmers.end(),dev_elmers.begin(),cutter<T>());
        thrust::host_vector<T> helmers=dev_elmers;
        std::copy(helmers.begin(),helmers.end(),elmers.begin()+base);
    }

    return nnzero;
}


template<typename T>
bool gCudaCounter<T>::estimatelmers(const std::map<indexType,T> gkmers, std::vector<T> & elmers,indexType base,bool doubleStrand) const{
    thrust::device_vector<T> dev_elmers(i_counter.i_nlmers,0.0);
    

        //indexType mlmer=lmer * i_counter.i_nmasks;
        for(auto gkmer=gkmers.begin();gkmer!=gkmers.end();gkmer++){
            indexType gkmer_kmer = gkmer->first  / i_counter.i_nmasks;
            indexType gkmer_mask = gkmer->first  % i_counter.i_nmasks;
            T val = gkmer->second;
            //indexType lmer_kmer=(gGappedKmerCounter<T>::i_mlmerGkmer[lmer * gGappedKmerCounter<T>::i_nmasks + gkmer_mask] / gGappedKmerCounter<T>::i_nmasks);
            auto Cbegin = thrust::make_permutation_iterator(
                i_mmcoeffs.begin(),
                thrust::make_transform_iterator(
                    thrust::make_permutation_iterator(
                        i_mlmerGkmer.begin(),
                        thrust::make_transform_iterator(
                            thrust::make_counting_iterator(0),
                            getLmerGkmer(i_counter.i_nmasks,gkmer_mask)
                        )
                    ),
                    getMismatches(i_counter.i_k,gkmer_kmer,i_counter.i_nmasks)
                )
            );
            thrust::transform(Cbegin,Cbegin+i_counter.i_nlmers,dev_elmers.begin(),dev_elmers.begin(),summer<T>(val));
        }

        
    if(doubleStrand){
        thrust::host_vector<T> helmers=dev_elmers;
        thrust::host_vector<T> hdselmers(i_counter.i_nDSlmers,0);
        for(indexType lmer=0;lmer<i_counter.i_nlmers;lmer++){
            hdselmers[i_counter.i_lmerDSLmer[lmer]] += helmers[lmer];
        }
        thrust::transform(hdselmers.begin(),hdselmers.end(),hdselmers.begin(),cutter<T>());
        std::copy(hdselmers.begin(),hdselmers.end(),elmers.begin()+base);
    }else{
        thrust::transform(dev_elmers.begin(),dev_elmers.end(),dev_elmers.begin(),cutter<T>());
        thrust::host_vector<T> helmers=dev_elmers;
        std::copy(helmers.begin(),helmers.end(),elmers.begin()+base);
    }

    return true;
}



template class geco::methods::kmers::gCudaCounter<double>;
template class geco::methods::kmers::gCudaCounter<float>;


