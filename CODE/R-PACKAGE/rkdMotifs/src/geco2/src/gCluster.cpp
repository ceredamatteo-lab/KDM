#include "gCluster.h"
#include <math.h>
bool fc_isnan(double x) { return isnan(x); }

#include "fastcluster/fastcluster_dm.cpp"
#include "fastcluster/fastcluster_R_dm.cpp"
#include "fastcluster/fastcluster_cd.h"

#include "gMethods_armadillo.h"
#include <boost/math/special_functions/binomial.hpp>

using namespace geco;
using namespace geco::methods;
using namespace geco::methods::cluster;
using namespace arma;
using namespace std;


indexType choose(indexType  n, indexType  k){
    indexType result;    
    if(n<k){
        result=0;
    }else{
        result=boost::math::binomial_coefficient<double>(n,k);
    }
    return result;
}

//BEGIN gHCluster definition

gHCluster::gHCluster(){
	i_npoints=0;
	i_merge=nullptr;
	i_height=nullptr;
    i_order=nullptr;
}

template<typename T>
void doDist(t_float * distmat,const gMSparse<T> & D,indexType npoints){
    SpMat<T> Dist=armadillo::convert<gMSparse,SpMat,T>(D);
#pragma omp parallel for
    for (size_t i=0; i<(npoints-1); i++) {
        size_t base=npoints*i-choose(i+1,2);
        for (size_t j=i+1; j<npoints; j++) {
            distmat[base++] = (t_float) Dist(i,j);
        }
    }
}

template<typename T>
void doDist(t_float * distmat,const gMDense<T> & D,indexType npoints){
    Mat<T> Dist=armadillo::convert<gMDense,Mat,T>(D);
#pragma omp parallel for
    for (size_t i=0; i<(npoints-1); i++) {
        size_t base=npoints*i-choose(i+1,2);
        for (size_t j=i+1; j<npoints; j++) {
            distmat[base++] = (t_float) Dist(i,j);
        }
    }
}

template<template<typename> class D,typename T>
gHCluster::gHCluster(const D<T> & dist,gHClusterMethod method):gHCluster(){
    t_float * distmat=nullptr;
    try{
        i_npoints=dist.i_nrows;
        distmat = new t_float[(i_npoints*(i_npoints-1))/2];
        doDist(distmat,dist,i_npoints);
        do_clustering(method,distmat);
        delete [] distmat;    
    }catch(exception & e){
        if(distmat!=nullptr) delete [] distmat;
        throw gMethodException(e.what());
    }
}
template geco::methods::cluster::gHCluster::gHCluster(const gMSparse<float> & dist,gHClusterMethod method);
template geco::methods::cluster::gHCluster::gHCluster(const gMSparse<double> & dist,gHClusterMethod method);
template geco::methods::cluster::gHCluster::gHCluster(const gMDense<float> & dist,gHClusterMethod method);
template geco::methods::cluster::gHCluster::gHCluster(const gMDense<double> & dist,gHClusterMethod method);


gHCluster::gHCluster(double * dist,indexType npoints,gHClusterMethod method):gHCluster(){
    i_npoints=npoints;
	t_float * distmat=dist;
    try{
        do_clustering(method,distmat);
    }catch(exception & e){
        throw gMethodException(e.what());
    }
}



gHCluster::~gHCluster(){
	if(i_merge!=nullptr){
		delete [] i_merge;
	}
	
	if(i_height!=nullptr){
		delete [] i_height;
	}
	
	if(i_order!=nullptr){
        delete [] i_order;
    }
}

void gHCluster::do_clustering(gHClusterMethod method,double * D){
    double * members = nullptr;
    if (method==METHOD_METR_AVERAGE || method==METHOD_METR_WARD_D || method==METHOD_METR_WARD_D2 || method==METHOD_METR_CENTROID){
        members = new double[i_npoints];
        for (int i=0; i<i_npoints; i++) members[i] = 1;
    }
    if (method==METHOD_METR_WARD_D2) {
        for (double * DD = D; DD!=D+static_cast<std::ptrdiff_t>(i_npoints)*(i_npoints-1)/2;++DD){
            *DD *= *DD;
        }
    }
    
    cluster_result Z2(i_npoints-1);
    switch (method) {
    case METHOD_METR_SINGLE:
    MST_linkage_core(i_npoints, D, Z2);
    break;
    case METHOD_METR_COMPLETE:
    NN_chain_core<method_codes::METHOD_METR_COMPLETE, t_float>(i_npoints, D, NULL, Z2);
    break;
    case METHOD_METR_AVERAGE:
    NN_chain_core<method_codes::METHOD_METR_AVERAGE, t_float>(i_npoints, D, members, Z2);
    break;
    case METHOD_METR_WEIGHTED:
    NN_chain_core<method_codes::METHOD_METR_WEIGHTED,t_float>(i_npoints, D, NULL, Z2);
    break;
    case METHOD_METR_WARD_D:
    case METHOD_METR_WARD_D2:
    NN_chain_core<method_codes::METHOD_METR_WARD, t_float>(i_npoints, D, members, Z2);
    break;
    case METHOD_METR_CENTROID:
    generic_linkage<method_codes::METHOD_METR_CENTROID, t_float>(i_npoints, D, members, Z2);
    break;
    case METHOD_METR_MEDIAN:
    generic_linkage<method_codes::METHOD_METR_MEDIAN, t_float>(i_npoints, D, NULL, Z2);
    break;
    default:
    throw gMethodException("gHCluster: Invalid method");
    }

    if(members!=nullptr) delete [] members;
    
	i_merge = new int[2*(i_npoints-1)];
    i_height = new double[i_npoints-1];
    i_order = new int[i_npoints];

    if (method==METHOD_METR_WARD_D2) {
        Z2.sqrt();
    }
    if (method==METHOD_METR_CENTROID || method==METHOD_METR_MEDIAN){
        generate_R_dendrogram<true>(i_merge, i_height, i_order, Z2, i_npoints);
    }else{
        generate_R_dendrogram<false>(i_merge, i_height, i_order, Z2, i_npoints);
    }
}


gMDense<indexType> gHCluster::cutree_k(int nclust) const{
	int * labels = new int[i_npoints];
	::cutree_k((int) i_npoints, i_merge, nclust, labels);
	gMDense<indexType> ret(i_npoints,1,labels);
	return ret;
}

gMDense<indexType> gHCluster::cutree_h(double height) const{
	int * labels = new int[i_npoints];
	::cutree_cdist(i_npoints,i_merge, i_height, height,labels);
    gMDense<indexType> ret(i_npoints,1,labels);
	return ret;
}

const int * gHCluster::getMerge() const{
    return i_merge;
}

const double * gHCluster::getHeight() const{
    return i_height;
}

const int * gHCluster::getOrder() const{
    return i_order;
}

indexType gHCluster::getNPoints() const{
    return i_npoints;
}


template<typename T>
pair<gMDense<T>,pair<gMDense<indexType>,gMDense<T>>> geco::methods::cluster::getSilhouette(const gMDense<T> & dist,const gMDense<indexType> & clusterLabels){
    arma::Mat<T> D=armadillo::convert<gMDense,arma::Mat,T>(dist);
    arma::Mat<indexType> tmpcl=armadillo::convert<gMDense,arma::Mat,indexType>(clusterLabels);
    arma::Col<indexType> labs=arma::unique(tmpcl);
    indexType nclusters=labs.n_elem;
    indexType npoints=tmpcl.n_elem;
    
    vector<arma::uvec> clids(nclusters);
    arma::uvec cl(tmpcl.n_elem);    
#pragma omp parallel for
    for(indexType i=0;i<nclusters;i++){
        clids[i]=find(tmpcl==labs[i]);
        uvec tmpf(clids[i].n_elem);
        tmpf.fill(i);
        cl(clids[i])=tmpf;
    }
    
    arma::Mat<T> s(npoints,1,arma::fill::zeros);

#pragma omp parallel for    
    for(indexType i=0;i<npoints;i++){
        T a=0;
        T b=arma::datum::inf;
        indexType cli=cl(i);
        Col<T> Di=D.row(i); 
        if(clids[cli].n_elem>1){
            arma::uvec ss=clids[cli];
            a=(arma::sum(Di(ss)))/(ss.n_elem-1);
        }
        
        for(indexType k=0;k<nclusters;k++){
            if(k!=cli){
                T mm=arma::mean(Di(clids[k]));
                if(mm<b){
                    b=mm;
                }
            }
        }
        if(clids[cli].n_elem>1){
            T M=max(a,b);
            if(M!=0){
                s(i) = (b-a) / M;
            }
        }
    }
    arma::Mat<T> sc(nclusters,1,arma::fill::zeros);
#pragma omp parallel for
    for(indexType i=0;i<nclusters;i++){
        sc(i) = mean(s(clids[i]));
    }
    
    return make_pair(move(armadillo::convert<arma::Mat,gMDense,T>(s)),make_pair(armadillo::convert<arma::Mat,gMDense,indexType>(labs),armadillo::convert<arma::Mat,gMDense,T>(sc)));
}

template pair<gMDense<float>,pair<gMDense<indexType>,gMDense<float>>>  geco::methods::cluster::getSilhouette(const gMDense<float> & dist,const gMDense<indexType> & clusterLabels);
template pair<gMDense<double>,pair<gMDense<indexType>,gMDense<double>>>  geco::methods::cluster::getSilhouette(const gMDense<double> & dist,const gMDense<indexType> & clusterLabels);



std::string geco::methods::cluster::getClusteringMethodName(gHClusterMethod method){
    std::string ret;
    switch (method) {
        case METHOD_METR_SINGLE:    ret="single";break;
        case METHOD_METR_COMPLETE:  ret="complete";break;
        case METHOD_METR_AVERAGE:   ret="average";break;
        case METHOD_METR_WEIGHTED:  ret="mcquitty";break;
        case METHOD_METR_WARD_D:    
        case METHOD_METR_WARD_D2:   ret="ward.D2";break;
        case METHOD_METR_CENTROID:  ret="centroid";break;
        case METHOD_METR_MEDIAN:    ret="median";break;
        default:                    ret="unknown method";break;
    }
    return ret;
}

//END gHClust definition

