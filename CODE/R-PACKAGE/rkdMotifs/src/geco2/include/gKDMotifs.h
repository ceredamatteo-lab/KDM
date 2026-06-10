#pragma once
#include "gPWM.h"
#include "gNMF.h"
#include "gKmers.h"
#include "gCluster.h"
#include "gBed.h"


namespace geco{
    namespace methods{
        namespace kdmotifs{

            
            template <typename entryTypeSelected,typename entryTypeExcluded>
            std::pair<geco::bed::gBedData<entryTypeSelected>,geco::bed::gBedData<geco::bed::gBedBaseEntry>> getBackground(const geco::bed::gBedData<entryTypeSelected> & selected,const geco::bed::gBedData<entryTypeExcluded> & excluded,unsigned int nregions,unsigned int ngcbins,unsigned int nrpbins,const std::string & fileName2bits);
            
            
            template<template<typename> class D, typename T>
            class kcountsNMFInit:public nmf::baseNMFInit<D,T>{
            public:
                cluster::gHClusterMethod i_method;
                bool i_bhattacharyya;
                
                kcountsNMFInit(cluster::gHClusterMethod method,bool bhattacharyya):nmf::baseNMFInit<D,T>(true,false),i_method(method),i_bhattacharyya(bhattacharyya){
                }
                
                bool operator ()(const D<T> & V, nmf::gNMF<D,T> & M, indexType rank) const;
            };
            
            class gKDMInfo{
            public:
                indexType i_l;
                indexType i_k;
                bool i_doubleStrand;
                bool i_estimated;
                gKDMInfo(indexType l, indexType k, bool doubleStrand, bool estimated):i_l(l),i_k(k),i_doubleStrand(doubleStrand),i_estimated(estimated){
                }
            };
            
            template<template<typename> class D, typename T>
            class kcountsNMFInit2:public nmf::baseNMFInit2<D,T>{
            public:
                cluster::gHClusterMethod i_method;
                bool i_bhattacharyya;
                
                kcountsNMFInit2(cluster::gHClusterMethod method,bool bhattacharyya):nmf::baseNMFInit2<D,T>(true,false),i_method(method),i_bhattacharyya(bhattacharyya){
                }
                
                bool operator ()(const nmf::dataFeeder<D,T> & feeder, nmf::gNMF2<D,T> & M, indexType rank) const;
            };
            
            template<typename T>
            class gKDistr:public gMSparse<T>,public gKDMInfo{
            public:
              gKDistr(const gMSparse<T> & M,indexType l, indexType k, bool doubleStrand, bool estimated):gMSparse<T>(M),gKDMInfo(l,k,doubleStrand,estimated){
              }
              
              gKDistr(const gMSparse<T> & M,const gKDMInfo &  info):gMSparse<T>(M),gKDMInfo(info){
              }
              
              gKDistr(gMSparse<T> && M,indexType l, indexType k, bool doubleStrand, bool estimated):gMSparse<T>(M),gKDMInfo(l,k,doubleStrand,estimated){
              }
              
              gKDistr(gMSparse<T> && M,const gKDMInfo &  info):gMSparse<T>(M),gKDMInfo(info){
              }
              
            };
            
            template<typename T>
            class gKDMotifs:public gMDense<T>,public gKDMInfo{
            public:
              gKDMotifs(const gMDense<T> & M,indexType l, indexType k, bool doubleStrand, bool estimated):gMDense<T>(M),gKDMInfo(l,k,doubleStrand,estimated){
              }
              
              gKDMotifs(const gMDense<T> & M, const gKDMInfo & info):gMDense<T>(M),gKDMInfo(info){
              }
              
              gKDMotifs(gMDense<T> && M,indexType l, indexType k, bool doubleStrand, bool estimated):gMDense<T>(M),gKDMInfo(l,k,doubleStrand,estimated){
              }
              
              gKDMotifs(gMDense<T> && M,const gKDMInfo & info):gMDense<T>(M),gKDMInfo(info){
              }
              
            };
            
            template<typename T>
            class gKDMFeatures:public gMDense<T>,public gKDMInfo{
            public:
              gKDMFeatures(const gMDense<T> & M,indexType l, indexType k, bool doubleStrand, bool estimated):gMDense<T>(M),gKDMInfo(l,k,doubleStrand,estimated){
              }
              
              gKDMFeatures(const gMDense<T> & M,const gKDMInfo & info):gMDense<T>(M),gKDMInfo(info){
              }
              
              
              gKDMFeatures(gMDense<T> && M,indexType l, indexType k, bool doubleStrand, bool estimated):gMDense<T>(M),gKDMInfo(l,k,doubleStrand,estimated){
              }
              
              gKDMFeatures(gMDense<T> && M,const gKDMInfo & info):gMDense<T>(M),gKDMInfo(info){
              }

            };
            
            template<typename T>
            class gKDModel:public gKDMotifs<T>{
            public:
                gMDense<T> i_coefficients;
                std::vector<T> i_intercept;
                gKDModel(const gKDMotifs<T> & W, const gMDense<T> & coefficients,std::vector<T> intercept):gKDMotifs<T>(W),i_coefficients(coefficients),i_intercept(intercept){
                }
            };
            
            
            template<typename T>
            gMDense<T> predict(const gKDModel<T> & model,const gKDMFeatures<T> & features);

            template<typename T>
            gMDense<T> predict(const gKDModel<T> & model,const gKDistr<T> & kdistr);

            template<typename T>
            gMDense<T> predict(const gKDModel<T> & model,const std::vector<std::string> & sequences,bool strict=false);
            
            template<typename T,typename entryType>
            gMDense<T> predict(const gKDModel<T> & model,const geco::bed::gBedData<entryType> & bedData, const std::string & genomeFile,bool strict);

            
            
            template<typename T>
            gKDistr<T> toKDistr(const gKDMotifs<T> & kmotifs);
            
            template<typename T> 
            gKDMotifs<T> toKDMotifs(const gKDistr<T> & kdistr);
            
            typedef enum{random=0,nndsvd=1,cluster=2,custom=3} initMethod;


            
            template<typename T>
            gKDistr<T> kdmDistr(const std::vector<std::string> & sequences,indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers=false,bool strict=false);
            
            template<typename T>
            gKDistr<T> kdmDistrProfile(const std::string & sequence,indexType win, indexType kmerlength,indexType ngaps,bool doubleStrand,bool estimatekmers=false,bool strict=false);
            
            template<typename T>
            gKDMotifs<T> kdmMotifs(const gKDistr<T> & kdistr,indexType rank, indexType np, initMethod imethod=cluster,const gMDense<T> & W0=gMDense<T>(),cluster::gHClusterMethod clustering_method=cluster::gHClusterMethod::METHOD_METR_COMPLETE, bool bhattacharyya=false,bool initOnly = false, T tolerance=1e-6,indexType maxIterations=20000);
            
            template<typename T>
            gKDMotifs<T> kdmMotifs2(const std::vector<std::string> & sequences,const gKDMInfo & info,indexType rank, indexType np, initMethod imethod=cluster,const gMDense<T> & W0=gMDense<T>(),cluster::gHClusterMethod clustering_method=cluster::gHClusterMethod::METHOD_METR_COMPLETE, bool bhattacharyya=false,bool initOnly = false, T tolerance=1e-6,indexType maxIterations=20000);
            
            template<typename T>
            gKDMotifs<T> kdmMotifs(const gKDMotifs<T> & motifs,indexType rank, indexType np, initMethod imethod=cluster,const gMDense<T> & W0=gMDense<T>(),cluster::gHClusterMethod clustering_method=cluster::gHClusterMethod::METHOD_METR_COMPLETE, bool bhattacharyya=false, bool initOnly = false, T tolerance=1e-6,indexType maxIterations=20000);
            
            template<typename T>
            gKDMotifs<T> kdmEstimateMotifs(const gKDMotifs<T> & motifs);
            
            template<typename T>
            gKDMotifs<T> kdmMotifsFromPWM(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand,indexType padSize);
            
            template<typename T>
            gKDMotifs<T> kdmMotifsFromPWM2(const geco::kmers::gPWMSet & PWMs,indexType kmerlength,indexType ngaps,bool doubleStrand);
            
            template<typename T>
            std::vector<T> kdmGetPWMThresholds(const geco::kmers::gPWMSet & PWMs,T alpha,T initialGranularity=1e-4,bool forcedGranularity=true,T maxGranularity=1e-10);
            
            
            template<typename T>
            gKDMotifs<T> kdmMotifsReduce(const gKDMotifs<T>& kdmotifs,indexType rank=0,T tolerance=1e-6,indexType maxIterations=20000);
            

            template<typename T>
            gKDMFeatures<T> kdmFeatures(const gKDMotifs<T> & kmotifs,const gKDistr<T> & kdistr);
            
            template<typename T>
            gKDMFeatures<T> kdmFeatures(const gKDMotifs<T> & kmotifs,const std::vector<std::string> & sequences, bool strict=false);
            
            template<typename T>
            gKDMFeatures<T> kdmFeaturesProfile(const gKDMotifs<T> & kmotifs,const std::string & sequence, indexType win, bool strict=false);
            
            
            template< template<typename> class D, typename T>
            gMDense<T> kdmSimilarity(const D<T> & D1,const D<T> & D2);
            
            template< template<typename> class D, typename T>
            gMDense<T> kdmDistance(const D<T> & D1,const D<T> & D2,bool bhattacharyya=false);
            
            template<typename T>
            std::pair<gMDense<T>,gMDense<T>> kdmModelsCorrelation(const gKDModel<T> & models1,const gKDModel<T> & models2,indexType nrep=0);
            
            template<typename T>
            std::pair<gMDense<T>,gMDense<T>> kdmModelsCorrelation2(const gKDModel<T> & models1,const gKDModel<T> & models2,indexType nrep=0);
            
            
            template<typename T>
            std::pair<gMDense<T>,gMDense<T>> kdmModelsCorrelation(const gKDModel<T> & models,indexType nrep=0);

            template<typename T>            
            std::pair<gMDense<T>,gMDense<T>> kdmModelsCorrelation(const gKDModel<T> & models,const gKDMotifs<T> & motifs,indexType nrep=0);
            
            template<typename T>
            gKDistr<T> kdmSelectKDistr(const gKDistr<T> & d,const std::vector<indexType> & selection);
            
            template<typename T>
            gKDMotifs<T> kdmSelectKDMotifs(const gKDMotifs<T> & m,const std::vector<indexType> & selection);
            
            template<typename T>
            gKDistr<T> kdmMergeKDistr(const gKDistr<T> & d1,const gKDistr<T> & d2);
            
            template<typename T>
            gKDMotifs<T> kdmMergeKDMotifs(const gKDMotifs<T> & d1,const gKDMotifs<T> & d2);
            
            template<typename T>
            gKDMotifs<T> kdmMergeKDMotifs(const gKDMotifs<T> & d1,const gMDense<indexType> & labels,indexType rank);
            
            template<typename T>
            std::pair<gMDense<T>,std::pair<gMDense<indexType>,gMDense<T>>> kdmSilhouette(const gMDense<T> & dist,const gMDense<indexType> & clusterLabels);
            
            template<typename T>
            std::pair<gMDense<T>,gMDense<indexType>> kdmRankInfo(const gMDense<T> & Distance,cluster::gHClusterMethod clustering_method);
            
            template<typename T>
            std::pair<gMDense<T>,std::pair<gMDense<T>,gMDense<indexType>>> kdmGetW0(const gMSparse<T> & counts,cluster::gHClusterMethod clustering_method,bool goodclusters, indexType minRank, indexType maxRank,bool bhattacharyya);            
            
            template<typename T>
            std::pair<gMDense<T>,std::vector<T>> getMotifsProfile2(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin,bool strict);
            
            template<typename T>
            std::pair<std::pair<std::vector<gMDense<T>>,std::vector<gMDense<T>>>,std::vector<T>> getMotifsProfile3(const gKDMotifs<T> & motifs,const std::vector<std::string> & pseq, const std::vector<std::string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict);

            template<typename T>
            std::pair<std::pair<std::vector<gMDense<T>>,std::vector<gMDense<T>>>,std::vector<T>> getMotifsProfile3(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin,bool strict);
            
            template<typename T>
            std::pair<std::pair<std::vector<gMDense<T>>,std::vector<gMDense<T>>>,std::vector<T>> getMotifsProfile4(const gKDMotifs<T> & motifs,const std::vector<std::string> & pseq, const std::vector<std::string> & nseq, unsigned int halfInterval,unsigned int halfWin,bool strict);

            template<typename T>
            std::pair<std::pair<std::vector<gMDense<T>>,std::vector<gMDense<T>>>,std::vector<T>> getMotifsProfile4(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin,bool strict);
            
            template <typename T> 
            gMDense<T> getMotifsProfileParameters(const std::pair<gMDense<T>,std::vector<T>> & profile,unsigned int halfInterval, unsigned int halfWin,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval);
            
            template<typename T>
            std::pair<gMDense<T>,std::pair<gMDense<T>,std::vector<T>>>  getMotifsProfile(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval,bool strict);
            
//this is a first attempt
            template<typename T>
            class motifProfileInfo{
            public:
                    std::vector<T> x;    
                    std::vector< gMDense<T> > pfeat;
                    std::vector< gMDense<T> > nfeat;
                    std::vector<T> thresholds;
                    std::vector<T> aurocs;
                    geco::methods::gMDense<int> phits;
                    geco::methods::gMDense<int> nhits;
                    geco::methods::gMDense<int> pcounts;
                    geco::methods::gMDense<int> ncounts;
            };
            
            template<typename T>
            motifProfileInfo<T> kdmProfileInfo(const gKDMotifs<T> & motifs,const std::vector<std::string> & pseq, const std::vector<std::string> & nseq, unsigned int halfInterval,unsigned int halfWin, T tolerance=1e-10,bool strict=false);

            template<typename T>
            motifProfileInfo<T> kdmProfileInfo(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin, T tolerance=1e-10,bool strict=false);

/*  QUESTE ERANO DI PROVA PER L'APPROSSIMAZIONE
            template<typename T>
            motifProfileInfo<T> kdmProfileInfo2(const gKDMotifs<T> & motifs,const std::vector<std::string> & pseq, const std::vector<std::string> & nseq, unsigned int halfInterval,unsigned int halfWin, T tolerance=1e-10,bool strict=false);

            template<typename T>
            motifProfileInfo<T> kdmProfileInfo2(const gKDMotifs<T> & motifs,const geco::bed::gBedData<geco::bed::gBedBaseEntry> & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin, T tolerance=1e-10,bool strict=false);
*/                        
            enum class tailType {
                LOWER,
                UPPER,
                TWO_TAILED
            };

            template<typename T>
            std::vector< std::vector<std::vector<T>> > motifsEnrichment(const gMDense<int> & pcounts,const gMDense<int> & ncounts,unsigned int npseqs, unsigned int nnseqs,tailType tail,bool symmetric);
            
        };
    };
};
