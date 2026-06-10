#pragma once
#include "gKmers.h"

#define THRUST_IGNORE_CUB_VERSION_CHECK
#include "thrust/device_vector.h"




namespace geco{
	namespace methods{
		namespace kmers{

            template<typename T>
            class gCudaCounter{
                const gLmerEstimator<T> & i_counter;
                thrust::device_vector<indexType> i_masks; 
                thrust::device_vector<T> i_mmcoeffs; 
                thrust::device_vector<indexType> i_mlmerGkmer;
                thrust::device_vector<indexType> i_lmerDSLmer;
            public:
                
                gCudaCounter(const gLmerEstimator<T> & counter);                
                bool estimatelmers(const std::map<size_t,std::list<unsigned char> > & nseq, std::vector<T> & elmers,indexType base,bool doubleStrand=false,bool strict=false) const;
                bool estimatelmers(const std::map<indexType,T> gkmers, std::vector<T> & elmers,indexType base,bool doubleStrand) const;
            };
        }
    }
}

