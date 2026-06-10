#pragma once
#include "gMethods.h"

namespace geco{
	namespace methods{
        namespace cluster{
        
            typedef enum{
                //single
                METHOD_METR_SINGLE           = 0,
                //complete
                METHOD_METR_COMPLETE         = 1,
                //average
                METHOD_METR_AVERAGE          = 2,
                //mcquitty
                METHOD_METR_WEIGHTED         = 3,
                //ward.D
                METHOD_METR_WARD             = 4,
                METHOD_METR_WARD_D           = METHOD_METR_WARD,
                //centroid
                METHOD_METR_CENTROID         = 5,
                //median
                METHOD_METR_MEDIAN           = 6,
                //ward.D2
                METHOD_METR_WARD_D2          = 7,
                
                MIN_METHOD_CODE              = 0,
                MAX_METHOD_CODE              = 7
            } gHClusterMethod;
        
            class gHCluster{
            public:
                template<template<typename> class D,typename T>
                gHCluster(const D<T> & dist,gHClusterMethod method=METHOD_METR_COMPLETE);

                gHCluster(double * dist,indexType npoints,gHClusterMethod method=METHOD_METR_COMPLETE);
                ~gHCluster();
                
                gMDense<indexType> cutree_k(int nclust) const;
                gMDense<indexType> cutree_h(double height) const;
                
                const int * getMerge() const;
                const double * getHeight() const;
                const int * getOrder() const;
                indexType getNPoints() const;
                
//                 template<template<typename> class D,typename T>
//                 D<T> getCentroids(const D<T> & data,indexType rank);
                
            private:
                gHCluster();
                void do_clustering(gHClusterMethod method,double * D);
                
                indexType i_npoints;
                int * i_merge;
                double * i_height;
                int * i_order;
                
            };
            
            template<typename T>
            std::pair<gMDense<T>,std::pair<gMDense<indexType>,gMDense<T>>> getSilhouette(const gMDense<T> & dist,const gMDense<indexType> & clusterLabels);

            
            std::string getClusteringMethodName(gHClusterMethod method);
        }
	}
}

