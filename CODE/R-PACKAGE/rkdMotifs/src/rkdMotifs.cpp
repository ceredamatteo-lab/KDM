// #include <geco2/include/gKDMotifs.h>
// #include <geco2/include/gBed.h>
// #include <geco2/include/gPWM.h>
// #include <geco2/include/gKmers.h>

#include <gKDMotifs.h>
#include <gBed.h>
#include <gPWM.h>
#include <gKmers.h>

#include "rGeco.h"

using namespace Rcpp;

using namespace geco;
using namespace geco::methods;
using namespace geco::methods::kdmotifs;
using namespace geco::methods::kmers;
using namespace geco::kmers;
using namespace geco::bed;

template <typename T>
gKDistr<T> convertKDistrFromR(const Rcpp::S4 & mat){
  if(!mat.is("kdistr")){
    Rcpp::stop("An object of class \"kdistr\" should be provided" );
  }
  return gKDistr<T>(convertFromR<T>(mat),mat.slot("l"),mat.slot("k"),mat.slot("doubleStrand"),mat.slot("estimated"));
}

template <typename T>
Rcpp::S4 convertFromG(const gKDistr<T> & mat){
  S4 s=convertFromG<T>((const gMSparse<T> &) mat,std::string("kdistr"));
  s.slot("l")   = IntegerVector(1,mat.i_l);
  s.slot("k")   = IntegerVector(1,mat.i_k);
  s.slot("doubleStrand")   = mat.i_doubleStrand;
  s.slot("estimated") = mat.i_estimated;
  return s;
}

template <typename T>
gKDMotifs<T> convertKDMotifsFromR(const Rcpp::S4 & mat){
  if(!mat.is("kdmotifs")){
    Rcpp::stop("An object of class \"kmotifs\" should be provided" );
  }
  NumericMatrix m=mat.slot(".Data");
  return gKDMotifs<T>(convertFromR<T>(m),mat.slot("l"),mat.slot("k"),mat.slot("doubleStrand"),mat.slot("estimated"));
}

template <typename T>
Rcpp::S4 convertFromG(const gKDMotifs<T> & mat,const std::string & name=std::string("")){
  std::string klass="kdmotifs";
  if(name.length()>0){
    klass=name;
  }
  S4 s(klass);
  s.slot(".Data")   = convertFromG<T>((const gMDense<T> &) mat);
  s.slot("l")   = IntegerVector(1,mat.i_l);
  s.slot("k")   = IntegerVector(1,mat.i_k);
  s.slot("doubleStrand")   = mat.i_doubleStrand;
  s.slot("estimated") = mat.i_estimated;
  return s;
}

template <typename T>
gKDMFeatures<T> convertKDMFeaturesFromR(const Rcpp::S4 & mat){
  if(!mat.is("kdfeatures")){
    Rcpp::stop("An object of class \"kdfeatures\" should be provided" );
  }
  NumericMatrix m=mat.slot(".Data");
  return gKDMFeatures<T>(convertFromR<T>(m),mat.slot("l"),mat.slot("k"),mat.slot("doubleStrand"),mat.slot("estimated"));
}

template <typename T>
Rcpp::S4 convertFromG(const gKDMFeatures<T> & mat){
  S4 s("kdfeatures");
  s.slot(".Data")   = convertFromG<T>((const gMDense<T> &) mat);
  s.slot("l")   = IntegerVector(1,mat.i_l);
  s.slot("k")   = IntegerVector(1,mat.i_k);
  s.slot("doubleStrand")   = mat.i_doubleStrand;
  s.slot("estimated") = mat.i_estimated;
  return s;
}


template <typename T>
gKDModel<T> convertKDModelFromR(const Rcpp::S4 & model){
  if(!model.is("kdmodel")){
    Rcpp::stop("An object of class \"kdmodel\" must be provided" );
  }
  NumericMatrix coefs=model.slot("coefficients");
  NumericVector interc=model.slot("intercept");
  std::vector<T> ii(interc.begin(),interc.end());
  return gKDModel<T>(convertKDMotifsFromR<T>(model),convertFromR<T>(coefs),ii);
}

template <typename T>
Rcpp::S4 convertFromG(const gKDModel<T> & model){
  S4 s=convertFromG<T>(model,"kdmodel");
  s.slot("coefficients") = convertFromG(model.i_coefficients);
  s.slot("intercept") = NumericVector(model.i_intercept.begin(),model.i_intercept.end());
  return s;
}



S4 convert(const gPWM & pwm){
  S4 s("kpwm");
  s.slot(".Data") = convert(pwm.getPWM());
  s.slot("background")   = convert<gScore>(pwm.getBackground());
  return s;  
}

gPWM convertToPWM(const S4 & mat){
  NumericMatrix pwm=mat.slot(".Data");
  NumericVector background=mat.slot("background");
  
  return gPWM(convert<gScore>(pwm),convert<gScore>(background));
}

S4 convert(const gPWMSet & pwms){
  S4 ret(std::string("kpwmset"));
  List pwmlist;
  for(auto pwm=pwms.begin();pwm!=pwms.end();pwm++){
    pwmlist.push_back(convert(pwm->second),pwm->first);
  }
  ret.slot(".Data")=pwmlist;
  return ret;
}

gPWMSet convertToPWMSet(const S4 & pwmset){
  gPWMSet ret;
  List pwms=pwmset.slot(".Data");
  CharacterVector names=pwmset.slot("names");
  for(int i=0;i<pwms.size();i++){
    S4 pwm=pwms(i);
    gPWM gpwm=convertToPWM(pwm);
    ret.addMatrix(gpwm,std::string(names[i]));
  }
  return ret;
}


DataFrame convert(const gBedData<gBedBaseEntry> & data){
  S4 ret(std::string("kbed"));
  CharacterVector chrom(data.size()),name(data.size()),strand(data.size());
  NumericVector start(data.size()),end(data.size()),score(data.size());
  unsigned int dd=0;
  for(auto r:data){
    chrom[dd]=r.i_chrom;
    start[dd]=r.i_start;
    end[dd]=r.i_end;
    name[dd]=r.i_name;
    score[dd]=r.i_score;
    strand[dd]=(r.i_isForward)?("+"):("-");
    dd++;
  }
  
  // ret.slot(".Data") = Rcpp::DataFrame::create(
  //   _["chrom"] = chrom,
  //   _["start"] = start,
  //   _["end"] = end,
  //   _["name"] = name,
  //   _["score"] = score,
  //   _["strand"] = strand
  // );
  // ret.slot("type")="bed6";
  // return ret;
  
  return  Rcpp::DataFrame::create(
    _["chrom"] = chrom,
    _["start"] = start,
    _["end"] = end,
    _["name"] = name,
    _["score"] = score,
    _["strand"] = strand
  );
  
}

gBedData<gBedBaseEntry> convertToGBedDataBase(const DataFrame & bed){
    gBedData<gBedBaseEntry> ret;
    CharacterVector chrom=bed["chrom"],name=bed["name"],strand=bed["strand"];
    NumericVector start=bed["start"],end=bed["end"],score=bed["score"];
    
    for(int i=0;i<bed.nrow();i++){
      gBedBaseEntry entry(std::string(chrom[i]),(gSize)start[i],(gSize)end[i],std::string(name[i]),(unsigned int)score[i],strand[i][0]);
      ret.push_back(entry);
    }
    return ret;
}

  
DataFrame convert(const gBedData<gNarrowPeak> & data){
  S4 ret(std::string("kbed"));    
  CharacterVector chrom(data.size()),name(data.size()),strand(data.size());
  NumericVector start(data.size()),end(data.size()),score(data.size());
  NumericVector signal(data.size()),pvalue(data.size()),qvalue(data.size()),peakpos(data.size());
  unsigned int dd=0;
  for(auto r:data){
    chrom[dd]=r.i_chrom;
    start[dd]=r.i_start;
    end[dd]=r.i_end;
    name[dd]=r.i_name;
    score[dd]=r.i_score;
    strand[dd]=(r.i_isForward)?("+"):("-");
    signal[dd]=r.i_signal;
    pvalue[dd]=r.i_pvalue;
    qvalue[dd]=r.i_qvalue;
    peakpos[dd]=r.i_peak;
    dd++;
  }
  Rcpp::DataFrame df=Rcpp::DataFrame::create(
    _["chrom"] = chrom,
    _["start"] = start,
    _["end"] = end,
    _["name"] = name,
    _["score"] = score,
    _["strand"] = strand,
    _["signal"] = signal,
    _["pvalue"] = pvalue,
    _["qvalue"] = qvalue,
    _["peak"] = peakpos
  );
  return df;
  // ret.slot(".Data")=df;
  // ret.slot("type")="narrowPeak";
  
  // return ret;
}

gBedData<gNarrowPeak> convertToGBedDataNarrow(const DataFrame & bed){
    gBedData<gNarrowPeak> ret;
    CharacterVector chrom=bed["chrom"],name=bed["name"],strand=bed["strand"];
    NumericVector start=bed["start"],end=bed["end"],score=bed["score"];
    NumericVector signal=bed["signal"],pvalue=bed["pvalue"],qvalue=bed["qvalue"],peakpos=bed["peak"];
    for(int i=0;i<bed.nrow();i++){
      gNarrowPeak entry(std::string(chrom[i]),(gSize)start[i],(gSize)end[i],std::string(name[i]),(unsigned int)score[i],strand[i][0],(gScore)signal[i],(gScore)pvalue[i],(gScore)qvalue[i],(gPos)peakpos[i]);
      ret.push_back(entry);
    }
    return ret;
}

//[[Rcpp::export]]
std::vector<std::string> test_convert_internal(const DataFrame & bed){
  gBedData<gBedBaseEntry> u=convertToGBedDataBase(bed);
  std::vector<std::string> seqs=u.getSequences("/adat/database/2bit/hg38.2bit");
  // std::vector<double> res(2);
  // res[0]=0;
  // res[1]=0;
  // for(auto e:u){
  //   if(e.i_isForward){
  //     res[0]++;
  //   }else{
  //     res[1]++;
  //   }
  // }
  return seqs;  
}


//[[Rcpp::export]]
DataFrame readBedFile_internal(const std::string & filename,const std::string & type){
  if(type=="bed6"){
    gBedData<gBedBaseEntry> data(filename);
    return convert(data);
  }else if(type=="narrowPeak"){
    gBedData<gNarrowPeak> data(filename);
    return convert(data);
  }else{
    Rcpp::stop("unsupportef bed format");
  }
}

//[[Rcpp::export]]
DataFrame remove_overlapping_internal(const DataFrame & bedData, const std::string & type){
  // CharacterVector names=bedData.names();
  // for(int i=0;i<names.size();i++){
  //   std::cout << i << "\t" << names[i] << std::endl;
  // }
  if(type=="bed6"){
    std::vector<gBedData<gBedBaseEntry>> ret=removeOverlappingEntries(convertToGBedDataBase(bedData));
    return convert(ret[0]);
  }else if(type=="narrowPeak"){
    std::vector<gBedData<gNarrowPeak>> ret=removeOverlappingEntries(convertToGBedDataNarrow(bedData));
    return convert(ret[0]);
  }else{
    Rcpp::stop("unsupported bed format");
  }
}

//[[Rcpp::export]]
List background_internal(const DataFrame & regions,const std::string & regionsType,const DataFrame & excludedRegions,const std::string & excludedRegionsType, int ngcbins,int nrpbins,int nregions,const std::string & referenceFileName){
  List ret;
  if(regionsType=="bed6"){
    if(excludedRegionsType=="bed6"){
      std::pair<gBedData<gBedBaseEntry>,gBedData<gBedBaseEntry>> results=getBackground<gBedBaseEntry,gBedBaseEntry>(convertToGBedDataBase(regions),convertToGBedDataBase(excludedRegions),nregions,ngcbins,nrpbins,referenceFileName);
      ret=List::create(_["sel"]=convert(results.first),_["bkg"]=convert(results.second));
    }else if(excludedRegionsType=="narrowPeak"){
      std::pair<gBedData<gBedBaseEntry>,gBedData<gBedBaseEntry>> results=getBackground<gBedBaseEntry,gNarrowPeak>(convertToGBedDataBase(regions),convertToGBedDataNarrow(excludedRegions),nregions,ngcbins,nrpbins,referenceFileName);
      ret=List::create(_["sel"]=convert(results.first),_["bkg"]=convert(results.second));
    }
  }else if(regionsType=="narrowPeak"){
    if(excludedRegionsType=="bed6"){
      std::pair<gBedData<gNarrowPeak>,gBedData<gBedBaseEntry>> results=getBackground<gNarrowPeak,gBedBaseEntry>(convertToGBedDataNarrow(regions),convertToGBedDataBase(excludedRegions),nregions,ngcbins,nrpbins,referenceFileName);
      ret=List::create(_["sel"]=convert(results.first),_["bkg"]=convert(results.second));
    }else if(excludedRegionsType=="narrowPeak"){
      std::pair<gBedData<gNarrowPeak>,gBedData<gBedBaseEntry>> results=getBackground<gNarrowPeak,gNarrowPeak>(convertToGBedDataNarrow(regions),convertToGBedDataNarrow(excludedRegions),nregions,ngcbins,nrpbins,referenceFileName);
      ret=List::create(_["sel"]=convert(results.first),_["bkg"]=convert(results.second));
    }
  }else{
    Rcpp::stop("unsupported bed format");
  }
  return ret;
}



//[[Rcpp::export]]
std::vector<std::string> get_sequences_internal(const DataFrame & bedData,const std::string& type,const std::string fileName){
  if(type=="bed6"){
    return convertToGBedDataBase(bedData).getSequences(fileName);
  }else if(type=="narrowPeak"){
    return convertToGBedDataNarrow(bedData).getSequences(fileName);
  }else{
    Rcpp::stop("unsupported bed format");
  }
  
}

//[[Rcpp::export]]
S4 kdistr_internal(const std::vector<std::string> & sequences,unsigned long kmerlength,unsigned long ngaps,bool doubleStrand,bool estimatekmers=false,bool strict=false,bool use_float=false){

  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmDistr<float>(sequences,kmerlength,ngaps,doubleStrand,estimatekmers,strict));
    
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmDistr<double>(sequences,kmerlength,ngaps,doubleStrand,estimatekmers,strict));
  }
}

//[[Rcpp::export]]
S4 kdistrProfile_internal(const std::string & sequence,unsigned long win,unsigned long kmerlength,unsigned long ngaps,bool doubleStrand,bool estimatekmers=false,bool strict=false,bool use_float=false){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmDistrProfile<float>(sequence,win,kmerlength,ngaps,doubleStrand,estimatekmers,strict));
    
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmDistrProfile<double>(sequence,win,kmerlength,ngaps,doubleStrand,estimatekmers,strict));
  }
}



//[[Rcpp::export]]
List kdmGkmersString_internal(unsigned long kmerlength, unsigned long ngaps,bool doubleStrand){
  std::pair<std::vector<std::string>,std::vector<std::string>> ret;
  if(ngaps>0){
    geco::methods::kmers::gGappedKmerCounter<float>  counter(kmerlength,kmerlength-ngaps,doubleStrand);
    ret=counter.getStrings();
  }else{
    geco::methods::kmers::gLmerCounter<float>  counter(kmerlength,doubleStrand);
    ret=counter.getStrings();
  }
  return List::create(_["fwd"]=ret.first,_["rev"]=ret.second);
}


//[[Rcpp::export]]
S4 kdmotifs_internal(const S4 & kdistr,unsigned long rank, unsigned long np, std::string initialization_method, const NumericMatrix & W0, const std::string & clustering_method, bool battacharyya, bool init_only, double tolerance, unsigned long maxIterations, bool use_float){
  std::map<std::string,initMethod> methods { 
    {std::string("random"),initMethod::random},
    {std::string("nndsvd"),initMethod::nndsvd},
    {std::string("cluster"),initMethod::cluster},
    {std::string("custom"),initMethod::custom}
  };
  initMethod im=methods.at(initialization_method);
  
  std::map<std::string,cluster::gHClusterMethod> cmethods { 
    {std::string("single"),cluster::gHClusterMethod::METHOD_METR_SINGLE}, 
    {std::string("complete"),cluster::gHClusterMethod::METHOD_METR_COMPLETE}, 
    {std::string("average"),cluster::gHClusterMethod::METHOD_METR_AVERAGE}, 
    {std::string("mcquitty"),cluster::gHClusterMethod::METHOD_METR_WEIGHTED}, 
    {std::string("ward.D"),cluster::gHClusterMethod::METHOD_METR_WARD_D}, 
    {std::string("centroid"),cluster::gHClusterMethod::METHOD_METR_CENTROID}, 
    {std::string("median"),cluster::gHClusterMethod::METHOD_METR_MEDIAN}, 
    {std::string("ward.D2"),cluster::gHClusterMethod::METHOD_METR_WARD_D2} 
  };
  cluster::gHClusterMethod cm=cmethods.at(clustering_method);
  
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmMotifs<float>(convertKDistrFromR<float>(kdistr),rank,np,im,convertFromR<float>(W0),cm,battacharyya,init_only,tolerance,maxIterations));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmMotifs<double>(convertKDistrFromR<double>(kdistr),rank,np,im,convertFromR<double>(W0),cm,battacharyya,init_only,tolerance,maxIterations));
  }
}

//[[Rcpp::export]]
S4 kdmotifs_from_sequences_internal(const std::vector<std::string> & sequences, unsigned long kmerlength,unsigned long ngaps,bool double_strand,bool estimate, unsigned long rank, unsigned long np, std::string initialization_method, const NumericMatrix & W0, const std::string & clustering_method, bool battacharyya, bool init_only, double tolerance, unsigned long maxIterations, bool use_float){
  std::map<std::string,initMethod> methods { 
    {std::string("random"),initMethod::random},
    {std::string("nndsvd"),initMethod::nndsvd},
    {std::string("cluster"),initMethod::cluster},
    {std::string("custom"),initMethod::custom}
  };
  initMethod im=methods.at(initialization_method);
  
  std::map<std::string,cluster::gHClusterMethod> cmethods { 
    {std::string("single"),cluster::gHClusterMethod::METHOD_METR_SINGLE}, 
    {std::string("complete"),cluster::gHClusterMethod::METHOD_METR_COMPLETE}, 
    {std::string("average"),cluster::gHClusterMethod::METHOD_METR_AVERAGE}, 
    {std::string("mcquitty"),cluster::gHClusterMethod::METHOD_METR_WEIGHTED}, 
    {std::string("ward.D"),cluster::gHClusterMethod::METHOD_METR_WARD_D}, 
    {std::string("centroid"),cluster::gHClusterMethod::METHOD_METR_CENTROID}, 
    {std::string("median"),cluster::gHClusterMethod::METHOD_METR_MEDIAN}, 
    {std::string("ward.D2"),cluster::gHClusterMethod::METHOD_METR_WARD_D2} 
  };
  cluster::gHClusterMethod cm=cmethods.at(clustering_method);
  
  geco::methods::kdmotifs::gKDMInfo info(kmerlength,ngaps,double_strand,estimate);
  
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmMotifs2<float>(sequences,info,rank,np,im,convertFromR<float>(W0),cm,battacharyya,init_only,tolerance,maxIterations));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmMotifs2<double>(sequences,info,rank,np,im,convertFromR<double>(W0),cm,battacharyya,init_only,tolerance,maxIterations));
  }
}


//[[Rcpp::export]]
S4 kdmotifs_estimate_internal(const S4 & motifs, bool use_float){
  if(use_float){
    gKDMotifs<float> A=convertKDMotifsFromR<float>(motifs);
    std::cout << "-1" << std::endl;
    gKDMotifs<float> B=geco::methods::kdmotifs::kdmEstimateMotifs(A);
    std::cout << "-2" << std::endl;
    return  convertFromG<float>(B);
    
    //return convertFromG<float>(geco::methods::kdmotifs::kdmEstimateMotifs(convertKDMotifsFromR<float>(motifs)));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmEstimateMotifs(convertKDMotifsFromR<double>(motifs)));
  }
}


//[[Rcpp::export]]
S4 kdmotifs_internal_from_motifs(const S4 & kdmotifs,unsigned long rank, unsigned long np, std::string initialization_method, const NumericMatrix & W0, const std::string & clustering_method, bool battacharyya, bool init_only, double tolerance, unsigned long maxIterations, bool use_float){
  std::map<std::string,initMethod> methods { {std::string("random"),initMethod::random},{std::string("nndsvd"),initMethod::nndsvd},{std::string("cluster"),initMethod::cluster},{std::string("custom"),initMethod::custom}};
  initMethod im=methods.at(initialization_method);
  std::map<std::string,cluster::gHClusterMethod> cmethods { {std::string("single"),cluster::gHClusterMethod::METHOD_METR_SINGLE}, {std::string("complete"),cluster::gHClusterMethod::METHOD_METR_COMPLETE}, {std::string("average"),cluster::gHClusterMethod::METHOD_METR_AVERAGE}, {std::string("mcquitty"),cluster::gHClusterMethod::METHOD_METR_WEIGHTED}, {std::string("ward.D"),cluster::gHClusterMethod::METHOD_METR_WARD_D}, {std::string("centroid"),cluster::gHClusterMethod::METHOD_METR_CENTROID}, {std::string("median"),cluster::gHClusterMethod::METHOD_METR_MEDIAN}, {std::string("ward.D2"),cluster::gHClusterMethod::METHOD_METR_WARD_D2} };
  cluster::gHClusterMethod cm=cmethods.at(clustering_method);
  
    if(use_float){
      return convertFromG<float>(geco::methods::kdmotifs::kdmMotifs<float>(convertKDMotifsFromR<float>(kdmotifs),rank,np,im,convertFromR<float>(W0),cm,battacharyya,init_only,tolerance,maxIterations));
      
    }else{
      return convertFromG<double>(geco::methods::kdmotifs::kdmMotifs<double>(convertKDMotifsFromR<double>(kdmotifs),rank,np,im,convertFromR<double>(W0),cm,battacharyya,init_only,tolerance,maxIterations));
    }
}



//[[Rcpp::export]]
S4 kdmotifs_merge_internal(const S4 & kdmotifs,const NumericMatrix & labels,unsigned int rank,bool use_float){
  gMDense<indexType> clustering=convertFromR<indexType>(labels);
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmMergeKDMotifs<float>(convertKDMotifsFromR<float>(kdmotifs),clustering,rank));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmMergeKDMotifs<double>(convertKDMotifsFromR<double>(kdmotifs),clustering,rank));
  }
}

//[[Rcpp::export]]
List kdmSilhouette_internal(NumericMatrix & distance,const NumericMatrix & clusterLabels,bool use_float){
  if(use_float){
    std::pair<gMDense<float>,std::pair<gMDense<indexType>,gMDense<float>>> res=geco::methods::kdmotifs::kdmSilhouette<float>(convertFromR<float>(distance),convertFromR<indexType>(clusterLabels));
    return List::create(_["silhouette"] = convertFromG<float>(res.first),_["label"] = convertFromG<indexType>(res.second.first),_["sil"] = convertFromG<float>(res.second.second));
  }else{
    std::pair<gMDense<double>,std::pair<gMDense<indexType>,gMDense<double>>> res=geco::methods::kdmotifs::kdmSilhouette<double>(convertFromR<double>(distance),convertFromR<indexType>(clusterLabels));
    return List::create(_["silhouette"] = convertFromG<double>(res.first),_["label"] = convertFromG<indexType>(res.second.first),_["sil"] = convertFromG<double>(res.second.second));
  }
}

//[[Rcpp::export]]
List kdmGetW0_internal(const S4 & counts, const std::string & clustering_method, bool good_clusters, unsigned int minRank, unsigned int maxRank, bool bhattacharyya,bool use_float){
  std::map<std::string,cluster::gHClusterMethod> cmethods { 
    {std::string("single"),cluster::gHClusterMethod::METHOD_METR_SINGLE}, 
    {std::string("complete"),cluster::gHClusterMethod::METHOD_METR_COMPLETE}, 
    {std::string("average"),cluster::gHClusterMethod::METHOD_METR_AVERAGE}, 
    {std::string("mcquitty"),cluster::gHClusterMethod::METHOD_METR_WEIGHTED}, 
    {std::string("ward.D"),cluster::gHClusterMethod::METHOD_METR_WARD_D}, 
    {std::string("centroid"),cluster::gHClusterMethod::METHOD_METR_CENTROID}, 
    {std::string("median"),cluster::gHClusterMethod::METHOD_METR_MEDIAN}, 
    {std::string("ward.D2"),cluster::gHClusterMethod::METHOD_METR_WARD_D2} 
  };
  cluster::gHClusterMethod cm=cmethods.at(clustering_method);
  
  if(use_float){
    std::pair< gMDense<float>, std::pair< gMDense<float>, gMDense<indexType>>> res=geco::methods::kdmotifs::kdmGetW0(convertKDistrFromR<float>(counts),cm,good_clusters,minRank,maxRank,bhattacharyya);
    return List::create(_["W0"] = convertFromG<float>(res.first),_["silhouette"] = convertFromG<float>(res.second.first),_["ngood"] = convertFromG<indexType>(res.second.second));
  }else{
    std::pair< gMDense<double>, std::pair< gMDense<double>, gMDense<indexType>>> res=geco::methods::kdmotifs::kdmGetW0(convertKDistrFromR<double>(counts),cm,good_clusters,minRank,maxRank,bhattacharyya);
    return List::create(_["W0"] = convertFromG<double>(res.first),_["silhouette"] = convertFromG<double>(res.second.first),_["ngood"] = convertFromG<indexType>(res.second.second));
  }

}

//[[Rcpp::export]]
List kdmRankInfo_internal(NumericMatrix & distance, const std::string & clustering_method, bool use_float){
  
  std::map<std::string,cluster::gHClusterMethod> cmethods { 
    {std::string("single"),cluster::gHClusterMethod::METHOD_METR_SINGLE}, 
    {std::string("complete"),cluster::gHClusterMethod::METHOD_METR_COMPLETE}, 
    {std::string("average"),cluster::gHClusterMethod::METHOD_METR_AVERAGE}, 
    {std::string("mcquitty"),cluster::gHClusterMethod::METHOD_METR_WEIGHTED}, 
    {std::string("ward.D"),cluster::gHClusterMethod::METHOD_METR_WARD_D}, 
    {std::string("centroid"),cluster::gHClusterMethod::METHOD_METR_CENTROID}, 
    {std::string("median"),cluster::gHClusterMethod::METHOD_METR_MEDIAN}, 
    {std::string("ward.D2"),cluster::gHClusterMethod::METHOD_METR_WARD_D2} 
  };
  cluster::gHClusterMethod cm=cmethods.at(clustering_method);
  
  if(use_float){
    std::pair<gMDense<float>,gMDense<indexType>> res=geco::methods::kdmotifs::kdmRankInfo<float>(convertFromR<float>(distance),cm);
    return List::create(_["silhouette"] = convertFromG<float>(res.first),_["ngood"] = convertFromG<indexType>(res.second));
  }else{
    std::pair<gMDense<double>,gMDense<indexType>> res=geco::methods::kdmotifs::kdmRankInfo<double>(convertFromR<double>(distance),cm);
    return List::create(_["silhouette"] = convertFromG<double>(res.first),_["ngood"] = convertFromG<indexType>(res.second));
  }
}



//[[Rcpp::export]]
S4 kdmotifsfrompwm_internal(const S4 & pwmset,unsigned long kmerlength,unsigned long ngaps, bool doubleStrand,unsigned long pad_size, bool use_float){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmMotifsFromPWM<float>(convertToPWMSet(pwmset),kmerlength,ngaps,doubleStrand,pad_size));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmMotifsFromPWM<double>(convertToPWMSet(pwmset),kmerlength,ngaps,doubleStrand,pad_size));
  }
}

//[[Rcpp::export]]
S4 kdmotifsfrompwm2_internal(const S4 & pwmset,unsigned long kmerlength,unsigned long ngaps, bool doubleStrand,bool use_float){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmMotifsFromPWM2<float>(convertToPWMSet(pwmset),kmerlength,ngaps,doubleStrand));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmMotifsFromPWM2<double>(convertToPWMSet(pwmset),kmerlength,ngaps,doubleStrand));
  }
}

//[[Rcpp::export]]
std::vector<double> kdmGetPWMThresholds_internal(const S4 & pwmset,double alpha,double initialGranularity,bool forcedGranularity,double maxGranularity){
    return kdmGetPWMThresholds<double>(convertToPWMSet(pwmset),alpha,initialGranularity,forcedGranularity,maxGranularity);
}


//[[Rcpp::export]]
S4 kdmotifsReduce_internal(const S4 & kmotifs,unsigned long rank, double tolerance=1e-6,unsigned long maxIterations=20000, bool use_float=false){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmMotifsReduce<float>(convertKDMotifsFromR<float>(kmotifs),rank,tolerance,maxIterations));
    
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmMotifsReduce<double>(convertKDMotifsFromR<double>(kmotifs),rank,tolerance,maxIterations));
  }
}



//[[Rcpp::export]]
S4 kdfeatures_internal_from_kdistr(const S4  & motifs, const S4 & kdistr, bool use_float){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmFeatures(convertKDMotifsFromR<float>(motifs),convertKDistrFromR<float>(kdistr)));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmFeatures(convertKDMotifsFromR<double>(motifs),convertKDistrFromR<double>(kdistr)));
  }
}

//[[Rcpp::export]]
S4 kdfeatures_internal_from_seqs(const S4  & motifs, const std::vector<std::string> & sequences, bool strict, bool use_float){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmFeatures(convertKDMotifsFromR<float>(motifs),sequences,strict));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmFeatures(convertKDMotifsFromR<double>(motifs),sequences,strict));
  }
}

//[[Rcpp::export]]
S4 kdfeatures_profile_internal(const S4 & motifs, const std::string & sequence, unsigned long win, bool strict, bool use_float){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmFeaturesProfile(convertKDMotifsFromR<float>(motifs),sequence,win,strict));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmFeaturesProfile(convertKDMotifsFromR<double>(motifs),sequence,win,strict));
  }
}


//[[Rcpp::export]]
S4 kdmToKdistr_internal(const S4 & m, bool use_float=true){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::toKDistr(convertKDMotifsFromR<float>(m)));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::toKDistr(convertKDMotifsFromR<double>(m)));
  }
}

//[[Rcpp::export]]
S4 kdmToKDMotifs_internal(const S4 & d, bool use_float=true){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::toKDMotifs(convertKDistrFromR<float>(d)));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::toKDMotifs(convertKDistrFromR<double>(d)));
  }
}



//[[Rcpp::export]]
NumericMatrix kdmSimilarity_kdistr_internal(const S4 & d1,const S4 & d2,bool use_float=false){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmSimilarity(convertKDistrFromR<float>(d1),convertKDistrFromR<float>(d2)));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmSimilarity(convertKDistrFromR<double>(d1),convertKDistrFromR<double>(d2)));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmSimilarity_kdmotifs_internal(const S4 & d1,const S4 & d2,bool use_float=false){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmSimilarity(convertKDMotifsFromR<float>(d1),convertKDMotifsFromR<float>(d2)));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmSimilarity(convertKDMotifsFromR<double>(d1),convertKDMotifsFromR<double>(d2)));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmDistance_kdistr_internal(const S4 & d1,const S4 & d2,bool bhattacharyya,bool use_float=false){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmDistance(convertKDistrFromR<float>(d1),convertKDistrFromR<float>(d2),bhattacharyya));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmDistance(convertKDistrFromR<double>(d1),convertKDistrFromR<double>(d2),bhattacharyya));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmDistance_kdmotifs_internal(const S4 & d1,const S4 & d2,bool bhattacharyya,bool use_float=false){
  if(use_float){
    return convertFromG<float>(geco::methods::kdmotifs::kdmDistance(convertKDMotifsFromR<float>(d1),convertKDMotifsFromR<float>(d2),bhattacharyya));
  }else{
    return convertFromG<double>(geco::methods::kdmotifs::kdmDistance(convertKDMotifsFromR<double>(d1),convertKDMotifsFromR<double>(d2),bhattacharyya));
  }
}


//[[Rcpp::export]]
NumericMatrix kdmSimilarity_kdistr_self_internal(const S4 & d1,bool use_float=false){
  if(use_float){
    gKDistr<float> M=convertKDistrFromR<float>(d1);
    return convertFromG<float>(geco::methods::kdmotifs::kdmSimilarity(M,M));
  }else{
    gKDistr<double> M=convertKDistrFromR<double>(d1);
    return convertFromG<double>(geco::methods::kdmotifs::kdmSimilarity(M,M));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmSimilarity_kdmotifs_self_internal(const S4 & d1,bool use_float=false){
  if(use_float){
    gKDMotifs<float> M=convertKDMotifsFromR<float>(d1);
    return convertFromG<float>(geco::methods::kdmotifs::kdmSimilarity(M,M));
  }else{
    gKDMotifs<double> M=convertKDMotifsFromR<double>(d1);
    return convertFromG<double>(geco::methods::kdmotifs::kdmSimilarity(M,M));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmDistance_kdistr_self_internal(const S4 & d1,bool bhattacharyya,bool use_float=false){
  if(use_float){
    gKDistr<float> M=convertKDistrFromR<float>(d1);
    return convertFromG<float>(geco::methods::kdmotifs::kdmDistance(M,M,bhattacharyya));
  }else{
    gKDistr<double> M=convertKDistrFromR<double>(d1);
    return convertFromG<double>(geco::methods::kdmotifs::kdmDistance(M,M,bhattacharyya));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmDistance_kdmotifs_self_internal(const S4 & d1,bool bhattacharyya,bool use_float=false){
  if(use_float){
    gKDMotifs<float> M=convertKDMotifsFromR<float>(d1);
    return convertFromG<float>(geco::methods::kdmotifs::kdmDistance(M,M,bhattacharyya));
  }else{
    gKDMotifs<double> M=convertKDMotifsFromR<double>(d1);
    return convertFromG<double>(geco::methods::kdmotifs::kdmDistance(M,M,bhattacharyya));
  }
}




//[[Rcpp::export]]
List kdmModelCorrelation_internal(const S4 & a,const S4 & b,unsigned long niter=0,bool use_float=false){
  if(a.is("kdmodel")){
    if(b.is("kdmodel")){
      NumericMatrix cc=b.slot("coefficients");
      if(cc.nrow()>0){
        if(use_float){
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<float>(a),convertKDModelFromR<float>(b),niter);
          return List::create(_["correlation"] = convertFromG<float>(res.first),_["p.values"] = convertFromG<float>(res.second));
        }else{
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<double>(a),convertKDModelFromR<double>(b),niter);    
          return List::create(_["correlation"] = convertFromG<double>(res.first),_["p.values"] = convertFromG<double>(res.second));
        }
      }else{
        if(use_float){
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<float>(a),niter);
          return List::create(_["correlation"] = convertFromG<float>(res.first),_["p.values"] = convertFromG<float>(res.second));
        }else{
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<double>(a),niter);    
          return List::create(_["correlation"] = convertFromG<double>(res.first),_["p.values"] = convertFromG<double>(res.second));
        }
      }
    }else if(b.is("kdmotifs")){
      if(use_float){
        auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<float>(a),convertKDMotifsFromR<float>(b),niter);
        return List::create(_["correlation"] = convertFromG<float>(res.first),_["p.values"] = convertFromG<float>(res.second));
      }else{
        auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<double>(a),convertKDMotifsFromR<double>(b),niter);    
        return List::create(_["correlation"] = convertFromG<double>(res.first),_["p.values"] = convertFromG<double>(res.second));
      }
      
      
    }else{
      Rcpp::stop("b must be either an object of class \"kdmodel\" or of class \"kdmotifs\" or NULL" );
    }
  }else{
    Rcpp::stop("a must be an object of class \"kdmodel\"" );
  }
}


//[[Rcpp::export]]
List kdmModelCorrelation2_internal(const S4 & a,const S4 & b,unsigned long niter=0,bool use_float=false){
  if(a.is("kdmodel")){
    if(b.is("kdmodel")){
      NumericMatrix cc=b.slot("coefficients");
      if(cc.nrow()>0){
        std::cout << "doing 2" << std::endl;
        if(use_float){
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation2(convertKDModelFromR<float>(a),convertKDModelFromR<float>(b),niter);
          return List::create(_["correlation"] = convertFromG<float>(res.first),_["p.values"] = convertFromG<float>(res.second));
        }else{
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation2(convertKDModelFromR<double>(a),convertKDModelFromR<double>(b),niter);
          return List::create(_["correlation"] = convertFromG<double>(res.first),_["p.values"] = convertFromG<double>(res.second));
        }
      }else{
        if(use_float){
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<float>(a),niter);
          return List::create(_["correlation"] = convertFromG<float>(res.first),_["p.values"] = convertFromG<float>(res.second));
        }else{
          auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<double>(a),niter);    
          return List::create(_["correlation"] = convertFromG<double>(res.first),_["p.values"] = convertFromG<double>(res.second));
        }
      }
    }else if(b.is("kdmotifs")){
      if(use_float){
        auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<float>(a),convertKDMotifsFromR<float>(b),niter);
        return List::create(_["correlation"] = convertFromG<float>(res.first),_["p.values"] = convertFromG<float>(res.second));
      }else{
        auto res=geco::methods::kdmotifs::kdmModelsCorrelation(convertKDModelFromR<double>(a),convertKDMotifsFromR<double>(b),niter);    
        return List::create(_["correlation"] = convertFromG<double>(res.first),_["p.values"] = convertFromG<double>(res.second));
      }
      
      
    }else{
      Rcpp::stop("b must be either an object of class \"kdmodel\" or of class \"kdmotifs\" or NULL" );
    }
  }else{
    Rcpp::stop("a must be an object of class \"kdmodel\"" );
  }
}

//[[Rcpp::export]]
NumericMatrix predict_internal_from_features(const S4 model,const S4 features,bool use_float=false){
  if(use_float){
    return convertFromG(predict(convertKDModelFromR<float>(model),convertKDMFeaturesFromR<float>(features)));
  }else{
    return convertFromG(predict(convertKDModelFromR<double>(model),convertKDMFeaturesFromR<double>(features)));
  }
}

//[[Rcpp::export]]
NumericMatrix predict_internal_from_kdistr(const S4 model,const S4 kdistr,bool use_float=false){
  if(use_float){
    return convertFromG(predict(convertKDModelFromR<float>(model),convertKDistrFromR<float>(kdistr)));
  }else{
    return convertFromG(predict(convertKDModelFromR<double>(model),convertKDistrFromR<double>(kdistr)));
  }
}

//[[Rcpp::export]]
NumericMatrix predict_internal_from_sequences(const S4 model,const std::vector<std::string> & sequences,bool strict=false,bool use_float=false){
  if(use_float){
    return convertFromG(predict(convertKDModelFromR<float>(model),sequences,strict));
  }else{
    return convertFromG(predict(convertKDModelFromR<double>(model),sequences,strict));
  }
}

//[[Rcpp::export]]
NumericMatrix predict_internal_from_bed(const S4 model,const DataFrame & bedData, const std::string & bedType, const std::string & filename,bool strict=false,bool use_float=false){
  if(bedType=="bed6"){
    if(use_float){
      return convertFromG(predict(convertKDModelFromR<float>(model),convertToGBedDataBase(bedData),filename,strict));
    }else{
      return convertFromG(predict(convertKDModelFromR<double>(model),convertToGBedDataBase(bedData),filename,strict));
    }
  }else if(bedType=="narrowPeak"){
    if(use_float){
      return convertFromG(predict(convertKDModelFromR<float>(model),convertToGBedDataNarrow(bedData),filename,strict));
    }else{
      return convertFromG(predict(convertKDModelFromR<double>(model),convertToGBedDataNarrow(bedData),filename,strict));
    }
  }else stop("bedType nust be either 'bed6' or 'narrowPeak'");
}


//[[Rcpp::export]]
S4 kdmLoadPWMSet(const std::string filename){
  return convert(readFromMemeFile(filename,gArray<gScore>(0.25,4,false),0.1));
}


//[[Rcpp::export]]
NumericMatrix kdmFeaturesFromPWMSet_internal(const std::vector<std::string> & sequences,const S4 & pwmset,bool double_strand,const std::string & mode){
  gPWMSet set=convertToPWMSet(pwmset);
  gMatrix<gScore> res;
  if(mode=="max"){
    res=set.getMaxScore(sequences,double_strand,gPWM::PSSMScore);
  }else if(mode=="avg"){
    res=set.getAvgScore(sequences,double_strand,gPWM::PSSMScore);
  }
  return convert(res);
}

//[[Rcpp::export]]
NumericMatrix kdmFeaturesProfileFromPWMSet_internal(const std::string & sequence,const S4 & pwmset){
  gPWMSet set=convertToPWMSet(pwmset);
  gMatrix<gScore> res=set.getScores(sequence,gPWM::PSSMScore);
  return convert(res);
}




//[[Rcpp::export]]
List kdmMotifProfiles2_internal(const DataFrame & bedData,const std::vector<unsigned int> & centers, const std::vector<unsigned int> & labels,const S4 & motifs, unsigned int halfInterval,unsigned int halfWin,std::string genomeFile,bool strict,bool use_float){
  gBedData<gBedBaseEntry> bed=convertToGBedDataBase(bedData);
  if(use_float){
    std::pair<gMDense<float>,std::vector<float>> res=geco::methods::kdmotifs::getMotifsProfile2(convertKDMotifsFromR<float>(motifs),bed,centers,labels,genomeFile,halfInterval,halfWin,strict);
    NumericMatrix P=convertFromG<float>(res.first);
    return List::create(_["Y"] = P,_["x"] = res.second);
  }else{
    std::pair<gMDense<double>,std::vector<double>> res=geco::methods::kdmotifs::getMotifsProfile2(convertKDMotifsFromR<double>(motifs),bed,centers,labels,genomeFile,halfInterval,halfWin,strict);
    NumericMatrix P=convertFromG<double>(res.first);
    return List::create(_["Y"] = P,_["x"] = res.second);
  }
}

//[[Rcpp::export]]
List kdmMotifProfiles3_internal(const DataFrame & bedData,const std::vector<unsigned int> & centers, const std::vector<unsigned int> & labels,const S4 & motifs, unsigned int halfInterval,unsigned int halfWin,std::string genomeFile,bool strict,bool use_float){
  gBedData<gBedBaseEntry> bed=convertToGBedDataBase(bedData);
  if(use_float){
    auto res=geco::methods::kdmotifs::getMotifsProfile3(convertKDMotifsFromR<float>(motifs),bed,centers,labels,genomeFile,halfInterval,halfWin,strict);
    List pfeat(res.first.first.size());
    List nfeat(res.first.second.size());
    for(indexType i=0;i<res.first.first.size();i++){
      // NumericMatrix rm=convertFromG<float>(res.first.first[i]);
      // pfeat[i]=rm;
      pfeat[i]=convertFromG<float>(res.first.first[i]);
    }
    for(indexType i=0;i<res.first.second.size();i++){
      // NumericMatrix rm=convertFromG<float>(res.first.second[i]);
      // nfeat[i]=rm;
      nfeat[i]=convertFromG<float>(res.first.second[i]);
    }

    return List::create(_["pfeat"] = pfeat,_["nfeat"] = nfeat,_["x"] = res.second);
  }else{
    auto res=geco::methods::kdmotifs::getMotifsProfile3(convertKDMotifsFromR<double>(motifs),bed,centers,labels,genomeFile,halfInterval,halfWin,strict);
    List pfeat(res.first.first.size());
    List nfeat(res.first.second.size());
    for(indexType i=0;i<res.first.first.size();i++){
      // NumericMatrix rm=convertFromG<double>(res.first.first[i]);
      // pfeat[i]=rm;
      pfeat[i]=convertFromG<double>(res.first.first[i]);
    }
    res.first.first.clear();
    for(indexType i=0;i<res.first.second.size();i++){
      // NumericMatrix rm=convertFromG<double>(res.first.second[i]);
      // nfeat[i]=rm;
      nfeat[i]=convertFromG<double>(res.first.second[i]);
    }
    res.first.second.clear();
    
    return List::create(_["pfeat"] = pfeat,_["nfeat"] = nfeat,_["x"] = res.second);
  }
}

//[[Rcpp::export]]
List kdmMotifProfiles3_seq_internal(const std::vector<std::string> & pseq,const std::vector<std::string> & nseq, const S4 & motifs, unsigned int halfInterval,unsigned int halfWin, bool strict,bool use_float){
  if(use_float){
    auto res=geco::methods::kdmotifs::getMotifsProfile3(convertKDMotifsFromR<float>(motifs),pseq,nseq,halfInterval,halfWin,strict);
    List pfeat(res.first.first.size());
    List nfeat(res.first.second.size());
    for(indexType i=0;i<res.first.first.size();i++){
      NumericMatrix rm=convertFromG<float>(res.first.first[i]);
      pfeat[i]=rm;
    }
    for(indexType i=0;i<res.first.second.size();i++){
      NumericMatrix rm=convertFromG<float>(res.first.second[i]);
      nfeat[i]=rm;
    }
    
    return List::create(_["pfeat"] = pfeat,_["nfeat"] = nfeat,_["x"] = res.second);
  }else{
    auto res=geco::methods::kdmotifs::getMotifsProfile3(convertKDMotifsFromR<double>(motifs),pseq,nseq,halfInterval,halfWin,strict);
    List pfeat(res.first.first.size());
    List nfeat(res.first.second.size());
    for(indexType i=0;i<res.first.first.size();i++){
      NumericMatrix rm=convertFromG<double>(res.first.first[i]);
      pfeat[i]=rm;
    }
    for(indexType i=0;i<res.first.second.size();i++){
      NumericMatrix rm=convertFromG<double>(res.first.second[i]);
      nfeat[i]=rm;
    }
    
    return List::create(_["pfeat"] = pfeat,_["nfeat"] = nfeat,_["x"] = res.second);
  }
}


//[[Rcpp::export]]
NumericMatrix kdmGetProfilesParameters_internal(const List & profile,unsigned int halfInterval, unsigned int halfWin, unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90Interval,bool use_float){
  if(use_float){
    NumericMatrix Y=profile["Y"];
    std::vector<float> x=profile["x"];
    std::pair<gMDense<float>,std::vector<float>> p=std::make_pair(convertFromR<float>(Y),x);
    gMDense<float> A=getMotifsProfileParameters(p,halfInterval,halfWin,nrep,maxSelectedPeak,half90Interval);
    return convertFromG<float>(A);
  }else{
    NumericMatrix Y=profile["Y"];
    std::vector<double> x=profile["x"];
    std::pair<gMDense<double>,std::vector<double>> p=std::make_pair(convertFromR<double>(Y),x);
    gMDense<double> A=getMotifsProfileParameters(p,halfInterval,halfWin,nrep,maxSelectedPeak,half90Interval);
    return convertFromG<double>(A);
  }
}

//[[Rcpp::export]]
List kdmMotifProfiles_internal(const DataFrame & bedData,const std::vector<unsigned int> & centers, const std::vector<unsigned int> & labels,const S4 & motifs, unsigned int interval,unsigned int win,unsigned int nrep,unsigned int maxSelectedPeak,unsigned int half90interval,std::string genomeFile,bool use_float){
  gBedData<gBedBaseEntry> bed=convertToGBedDataBase(bedData);
  if(use_float){
    auto res=geco::methods::kdmotifs::getMotifsProfile(convertKDMotifsFromR<float>(motifs),bed,centers,labels,genomeFile,interval,win,nrep,maxSelectedPeak,half90interval,false);
    NumericMatrix A=convertFromG<float>(res.first);
    NumericMatrix P=convertFromG<float>(res.second.first);
    return List::create(_["parameters"] = A,_["profile"] = List::create(_["Y"] = P,_["x"] = res.second.second));
  }else{
    auto res=geco::methods::kdmotifs::getMotifsProfile(convertKDMotifsFromR<double>(motifs),bed,centers,labels,genomeFile,interval,win,nrep,maxSelectedPeak,half90interval,false);
    NumericMatrix A=convertFromG<double>(res.first);
    NumericMatrix P=convertFromG<double>(res.second.first);
    return List::create(_["parameters"] = A,_["profile"] = List::create(_["Y"] = P,_["x"] = res.second.second));
  }
}

//[[Rcpp::export]]
NumericMatrix kdmGetPWM_internal(const S4 & W,unsigned len,const std::vector<double> & pwm0,double tol,unsigned maxIter,bool use_float){
  if(use_float){
    gKDMotifs<float> w=convertKDMotifsFromR<float>(W);
    if(w.i_k==0){
      geco::methods::kmers::gLmerCounter<float> counter(w.i_l,w.i_doubleStrand);
      return convertFromG<float>(counter.getPWM(w,len,pwm0,tol,maxIter));
    }else{
      geco::methods::kmers::gGappedKmerCounter<float> counter(w.i_l,w.i_l-w.i_k,w.i_doubleStrand);
      return convertFromG<float>(counter.getPWM(w,len,pwm0,tol,maxIter));
    }
  }else{
    gKDMotifs<double> w=convertKDMotifsFromR<double>(W);
    if(w.i_k==0){
      geco::methods::kmers::gLmerCounter<double> counter(w.i_l,w.i_doubleStrand);
      return convertFromG<double>(counter.getPWM(w,len,pwm0,tol,maxIter));
    }else{
      geco::methods::kmers::gGappedKmerCounter<double> counter(w.i_l,w.i_l-w.i_k,w.i_doubleStrand);
      return convertFromG<double>(counter.getPWM(w,len,pwm0,tol,maxIter));
    }
  } 
}


//[[Rcpp::export]]
List kdmProfileInfo_internal(const S4 & motifs,const DataFrame & regions,const std::vector<unsigned int> & centers,const std::vector<unsigned int> & labels,const std::string & genomeFile,unsigned int halfInterval,unsigned int halfWin,double tolerance,bool strict,bool use_float){
  gBedData<gBedBaseEntry> bed=convertToGBedDataBase(regions);
  if(use_float){
    auto W = convertKDMotifsFromR<float>(motifs);
    geco::methods::kdmotifs::motifProfileInfo<float> res=geco::methods::kdmotifs::kdmProfileInfo<float>(W,bed,centers,labels,genomeFile,halfInterval,halfWin,(float) tolerance,strict);
    List pfeat(res.pfeat.size());
    List nfeat(res.nfeat.size());
    IntegerMatrix phits=convertFromG2<int,IntegerMatrix>(res.phits);
    IntegerMatrix nhits=convertFromG2<int,IntegerMatrix>(res.nhits);
    
    IntegerMatrix pcounts=convertFromG2<int,IntegerMatrix>(res.pcounts);
    IntegerMatrix ncounts=convertFromG2<int,IntegerMatrix>(res.ncounts);
    
    for(indexType i=0;i<res.pfeat.size();i++){
      NumericMatrix rm=convertFromG<float>(res.pfeat[i]);
      pfeat[i]=rm;
    }
    for(indexType i=0;i<res.nfeat.size();i++){
      NumericMatrix rm=convertFromG<float>(res.nfeat[i]);
      nfeat[i]=rm;
    }
    return List::create(_["x"] = res.x, _["pfeat"] = pfeat,_["nfeat"] = nfeat, _["thresholds"] = res.thresholds, _["phits"]=phits,_["nhits"]=nhits,_["pcounts"]=pcounts,_["ncounts"]=ncounts);
  }else{
    auto W = convertKDMotifsFromR<double>(motifs);
    geco::methods::kdmotifs::motifProfileInfo<double> res=geco::methods::kdmotifs::kdmProfileInfo<double>(W,bed,centers,labels,genomeFile,halfInterval,halfWin,tolerance,strict);
    List pfeat(res.pfeat.size());
    List nfeat(res.nfeat.size());
    IntegerMatrix phits=convertFromG2<int,IntegerMatrix>(res.phits);
    IntegerMatrix nhits=convertFromG2<int,IntegerMatrix>(res.nhits);
    
    IntegerMatrix pcounts=convertFromG2<int,IntegerMatrix>(res.pcounts);
    IntegerMatrix ncounts=convertFromG2<int,IntegerMatrix>(res.ncounts);
  

    for(indexType i=0;i<res.pfeat.size();i++){
      NumericMatrix rm=convertFromG<double>(res.pfeat[i]);
      pfeat[i]=rm;
    }
    for(indexType i=0;i<res.nfeat.size();i++){
      NumericMatrix rm=convertFromG<double>(res.nfeat[i]);
      nfeat[i]=rm;
    }
    return List::create(_["x"] = res.x, _["pfeat"] = pfeat,_["nfeat"] = nfeat, _["thresholds"] = res.thresholds, _["aurocs"] = res.aurocs, _["phits"]=phits,_["nhits"]=nhits,_["pcounts"]=pcounts,_["ncounts"]=ncounts);
  }
}

//[[Rcpp::export]]
List kdmProfileInfoSeq_internal(const S4 & motifs,const std::vector<std::string> & pseq, const std::vector<std::string> & nseq, unsigned int halfInterval,unsigned int halfWin,double tolerance,bool strict,bool use_float){
  if(use_float){
    auto W = convertKDMotifsFromR<float>(motifs);
    geco::methods::kdmotifs::motifProfileInfo<float> res=geco::methods::kdmotifs::kdmProfileInfo<float>(W,pseq,nseq,halfInterval,halfWin,(float) tolerance,strict);
    List pfeat(res.pfeat.size());
    List nfeat(res.nfeat.size());
    IntegerMatrix phits=convertFromG2<int,IntegerMatrix>(res.phits);
    IntegerMatrix nhits=convertFromG2<int,IntegerMatrix>(res.nhits);
    
    IntegerMatrix pcounts=convertFromG2<int,IntegerMatrix>(res.pcounts);
    IntegerMatrix ncounts=convertFromG2<int,IntegerMatrix>(res.ncounts);
    
    for(indexType i=0;i<res.pfeat.size();i++){
      NumericMatrix rm=convertFromG<float>(res.pfeat[i]);
      pfeat[i]=rm;
    }
    for(indexType i=0;i<res.nfeat.size();i++){
      NumericMatrix rm=convertFromG<float>(res.nfeat[i]);
      nfeat[i]=rm;
    }
    return List::create(_["x"] = res.x, _["pfeat"] = pfeat,_["nfeat"] = nfeat, _["thresholds"] = res.thresholds, _["phits"]=phits,_["nhits"]=nhits,_["pcounts"]=pcounts,_["ncounts"]=ncounts);
  }else{
    auto W = convertKDMotifsFromR<double>(motifs);
    geco::methods::kdmotifs::motifProfileInfo<double> res=geco::methods::kdmotifs::kdmProfileInfo<double>(W,pseq,nseq,halfInterval,halfWin,tolerance,strict);
    List pfeat(res.pfeat.size());
    List nfeat(res.nfeat.size());
    IntegerMatrix phits=convertFromG2<int,IntegerMatrix>(res.phits);
    IntegerMatrix nhits=convertFromG2<int,IntegerMatrix>(res.nhits);
    
    IntegerMatrix pcounts=convertFromG2<int,IntegerMatrix>(res.pcounts);
    IntegerMatrix ncounts=convertFromG2<int,IntegerMatrix>(res.ncounts);
    
    
    for(indexType i=0;i<res.pfeat.size();i++){
      NumericMatrix rm=convertFromG<double>(res.pfeat[i]);
      pfeat[i]=rm;
    }
    for(indexType i=0;i<res.nfeat.size();i++){
      NumericMatrix rm=convertFromG<double>(res.nfeat[i]);
      nfeat[i]=rm;
    }
    return List::create(_["x"] = res.x, _["pfeat"] = pfeat,_["nfeat"] = nfeat, _["thresholds"] = res.thresholds, _["aurocs"] = res.aurocs, _["phits"]=phits,_["nhits"]=nhits,_["pcounts"]=pcounts,_["ncounts"]=ncounts);
  }
}


//[[Rcpp::export]]
List kdmCentrimo_internal(const IntegerMatrix & pcounts,const IntegerMatrix & ncounts,unsigned int npseqs,unsigned int nnseqs, std::string tail,bool symmetric,bool use_float){
  std::map<std::string,geco::methods::kdmotifs::tailType> ttypes { 
    {std::string("two"),geco::methods::kdmotifs::tailType::TWO_TAILED},
    {std::string("lower"),geco::methods::kdmotifs::tailType::LOWER}, 
    {std::string("upper"),geco::methods::kdmotifs::tailType::UPPER} 

  };
  geco::methods::kdmotifs::tailType tt=ttypes.at(tail);

  gMDense<int> cpcounts = convertFromR2<int,IntegerMatrix>(pcounts);
  gMDense<int> cncounts = convertFromR2<int,IntegerMatrix>(ncounts);
  
  if(use_float){
    auto res=geco::methods::kdmotifs::motifsEnrichment<float>(cpcounts,cncounts,npseqs,nnseqs,tt,symmetric);
    List retval=List::create();
    retval.assign(res.begin(),res.end());
    return retval;
  }else{
    auto  res=geco::methods::kdmotifs::motifsEnrichment<double>(cpcounts,cncounts,npseqs,nnseqs,tt,symmetric);
    List retval=List::create();
    retval.assign(res.begin(),res.end());
    return retval;
  }
}
