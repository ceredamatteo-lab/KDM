kdistr <- setClass(
  # Set the name for the class
  "kdistr",
  
  # Define the slots
  slots = c(
    l = "integer",
    k = "integer",
    doubleStrand = "logical",
    estimated = "logical"
  ),
  contains = "dgCMatrix"
)

setMethod("[", c("kdistr"),function(x, i, j, ..., drop=TRUE){
  cat("MYMETHOD kdistr\n")
  if(!missing(i)){
    stop("cannot subset kdistr rows")
  }
  x<-kdistr(x@.Data[i,j,drop=FALSE],l=x@l,k=x@k,doubleStrand=x@doubleStrand,estimated=x@estimated)
  x
})



kdmotif <- setClass(
  # Set the name for the class
  "kdmotifs",
  
  # Define the slots
  slots = c(
    l = "integer",
    k = "integer",
    doubleStrand = "logical",
    estimated = "logical"
  ),
  contains = "matrix"
)

setMethod("show", "kdmotifs",function(object){
  cat("kdm motifs matrix\n")
  cat("\t#N motifs: ",ncol(object@.Data),"\n",sep="")
  cat("\tKmer length: ",object@l,"\n",sep="")
  cat("\tN# gaps: ",object@k,"\n",sep="")
  cat("\tDouble strand: ",object@doubleStrand,"\n",sep="")
})


setMethod("[", c("kdmotifs"),function(x, i, j, ..., drop=TRUE){
  if(!missing(i)){
    stop("cannot subset kdmotifs rows")
    }
  x@.Data<-x@.Data[i,j,drop=FALSE]
  ## less clever: update slot, return instance
  ## x@slt = x@slt[i]
  ## x
  ## clever: by default initialize is a copy constructor, too
  x
})

kdfeatures <- setClass(
  # Set the name for the class
  "kdfeatures",
  
  # Define the slots
  slots = c(
    l = "integer",
    k = "integer",
    doubleStrand = "logical",
    estimated = "logical"
  ),
  contains = "matrix"
)

setMethod("[", c("kdfeatures"),
          ## we won't support subsetting on j; dispatching on 'drop' doesn't
          ## make sense (to me), so in rebellion we'll quietly ignore it.
          function(x, i, j, ..., drop=TRUE)
          {
            x@.Data<-x@.Data[i,j,drop=FALSE]
            ## less clever: update slot, return instance
            ## x@slt = x@slt[i]
            ## x
            ## clever: by default initialize is a copy constructor, too
            x
          })





kdmodel <- setClass(
  # Set the name for the class
  "kdmodel",
  
  # Define the slots
  slots = c(
    coefficients = "matrix",
    intercept = "numeric"
  ),
  contains = "kdmotifs"
)


setMethod("show", c("kdmodel"),function(object){
  cat("kdm Model\n")
  cat("\tOutcomes: ",ncol(object@coefficients),"\n",sep="")
  cat("\tPredictors: ")
  show(as.kdmotifs(object))
})

setAs("kdmodel", "kdmotifs",function(from){
  m<-new("kdmotifs",from@.Data,l=from@l,k=from@k,doubleStrand=from@doubleStrand,estimated=from@estimated)
  m
})


setMethod("[", c("kdmodel"),
          ## we won't support subsetting on j; dispatching on 'drop' doesn't
          ## make sense (to me), so in rebellion we'll quietly ignore it.
          function(x, i, j, ..., drop=TRUE)
          {
            stop("cannot subset kdmodel")
          })

as.kdmotifs<-function(x){
  m<-new("kdmotifs",x@.Data,l=x@l,k=x@k,doubleStrand=x@doubleStrand,estimated=x@estimated)
  m
}


kpwm<-setClass(
  "kpwm",
  slots=c(
    background="numeric"
  ),
  contains = "matrix"
)

kpwmset<-setClass(
  "kpwmset",
  contains = "list"
)

kbed<-setClass(
  "kbed",
  slots=c(
    type = "character"
  ),
  contains = "data.frame"
)

setMethod("[", c("kbed"),
          ## we won't support subsetting on j; dispatching on 'drop' doesn't
          ## make sense (to me), so in rebellion we'll quietly ignore it.
          function(x, i, j, ..., drop=TRUE)
          {
            if(!missing(j)){
              stop("cannot subset kbed columns")
            }
            kbed(data.frame(x)[i,],type=x@type)
          })




.onLoad <- function(libname, pkgname){
  initPackage();
}


kdmReadBedFile<-function(filename,type,force.names=T){
  #df=readBedFile_internal(filename,type);
  kbed(readBedFile_internal(filename,type),type=type)
}

kdmTestConvertNarrowPeak<-function(filename){
  bed<-kdmReadBedFile(filename,type="narrowPeak")
  rr<-test_convert_internal(bed)
  rr
}


kdmRemoveOverlappingEntries<-function(bedData){
    kbed(remove_overlapping_internal(bedData,bedData@type),type=bedData@type)
  }

kdmGetSequence<-function(bedData,file.name){
  get_sequences_internal(bedData,bedData@type,file.name)
}

kdmBackground<-function(selected,excluded,nregions,ref,n.gcbins=20,n.rpbins=20){
  tret=background_internal(selected,selected@type,excluded,excluded@type,n.gcbins,n.rpbins,nregions,ref);
  list(sel=kbed(tret[[1]],type=selected@type),bkg=kbed(tret[[2]],type="bed6"))
}

kdmDistr<-function(data,kmerlength,ngaps,double.strand,estimate=FALSE,strict=FALSE,use.float=FALSE){
  res<-c();
  if(is.character(data)){
    res<-kdistr_internal(data,kmerlength,ngaps,double.strand,estimate,strict,use.float)
  }else{
    stop("data must be a sequence array")
  }
  res
}

kdmDistrProfile<-function(sequence,win,kmerlength,ngaps,double.strand,estimate=FALSE,strict=FALSE,use.float=FALSE){
  res<-c();
  if(is.character(sequence)){
    res<-kdistrProfile_internal(sequence,win,kmerlength,ngaps,double.strand,estimate,strict,use.float)
  }else{
    stop("data must be a sequence array")
  }
  res
}


kdmGkmersString<-function(kmerlength,ngaps,double.strand){
  ss<-kdmGkmersString_internal(kmerlength,ngaps,double.strand)
  if(double.strand){
    as.data.frame(ss)
  }else{
    as.data.frame(ss[[1]])
  }
}

kdmMotifs<-function(kdistr,rank,np,initialization,W0=NULL,clustering.method="complete",bhattacharyya=FALSE,init.only=FALSE,tolerance=1e-6,max.iter=20000,use.float=FALSE){
  res<-c();
  if(class(kdistr)=="kdistr"){
    if(!is.null(W0)){
      res<-kdmotifs_internal(kdistr,rank,np,initialization,W0,clustering.method,bhattacharyya,init.only,tolerance,max.iter,use.float)
    }else{
      res<-kdmotifs_internal(kdistr,rank,np,initialization,matrix(0,0,0),clustering.method,bhattacharyya,init.only,tolerance,max.iter,use.float)
    }
  }else{
    if(!is.null(W0)){
      res<-kdmotifs_internal_from_motifs(kdistr,rank,np,initialization,W0,clustering.method,bhattacharyya,init.only,tolerance,max.iter,use.float)
    }else{
      res<-kdmotifs_internal_from_motifs(kdistr,rank,np,initialization,matrix(0,0,0),clustering.method,bhattacharyya,init.only,tolerance,max.iter,use.float)
    }
    
  }
    res
}

kdmMotifsFromSequences<-function(sequences,kmerlength,ngaps,double.strand,estimate=FALSE,rank,np,initialization,W0=NULL,clustering.method="complete",bhattacharyya=FALSE,init.only=FALSE,tolerance=1e-6,max.iter=20000,use.float=FALSE){
  res<-c();
  if(!is.null(W0)){
    res<-kdmotifs_from_sequences_internal(sequences,kmerlength,ngaps,double.strand,estimate=FALSE,rank,np,initialization,W0,clustering.method,bhattacharyya,init.only,tolerance,max.iter,use.float)
  }else{
    res<-kdmotifs_from_sequences_internal(sequences,kmerlength,ngaps,double.strand,estimate=FALSE,rank,np,initialization,matrix(0,0,0),clustering.method,bhattacharyya,init.only,tolerance,max.iter,use.float)
  }
  res
}




kdmEstimateMotifs<-function(motifs,use.float=FALSE){
  kdmotifs_estimate_internal(motifs,use.float)
}

kdmMotifsFromPWM<-function(pwmSet,kmerlength,ngaps,double.strand=FALSE,pad.size=kmerlength,use.float=FALSE){
  res<-kdmotifsfrompwm_internal(pwmSet,kmerlength,ngaps,double.strand,pad.size,use.float)
  res
}

kdmMotifsFromPWM2<-function(pwmSet,kmerlength,ngaps,double.strand=FALSE,use.float=FALSE){
  res<-kdmotifsfrompwm2_internal(pwmSet,kmerlength,ngaps,double.strand,use.float)
  res
}

kdmGetPWMThresholds<-function(pwmSet,alpha,initialGranularity=1e-4,forcedGranularity=TRUE,maxGranularity=1e-10){
  thresholds<-kdmGetPWMThresholds_internal(pwmSet,alpha,initialGranularity,forcedGranularity,maxGranularity)
  thresholds
}


kdmMotifsReduce<-function(motifs,rank,tolerance=1e-6,max.iter=20000,use.float=FALSE){
  res<-kdmotifsReduce_internal(motifs,rank,tolerance,max.iter,use.float)
  res
}

kdmBuildModel<-function(motifs,coeffmatrix,intercept,remove.null.coeff=TRUE){
  if(!inherits(motifs,"kdmotifs")){
    stop("motifs must be coercible to a 'kdmotifs' object")
  }
  model<-NULL
  dc<-dim(coeffmatrix)
  dm<-dim(motifs)
  if(dc[1]!=dm[2]){
    stop("number of coefficients must equal the number of motifs")
  }
  if(missing(intercept)){
    intercept<-rep(0,dc[2])
  }else{
    if(length(intercept)!=dc[2]){
      stop("number of intercepts must be equal to the coedfficients columns")
    }
  }
  if(remove.null.coeff){
    sel<-which(apply(coeffmatrix,1,function(a){sum(a!=0)})>0)
    model<-new("kdmodel",motifs[,sel],coefficients=coeffmatrix[sel,,drop=F],intercept=intercept);
  }else{
    model<-new("kdmodel",motifs,coefficients=coeffmatrix,intercept=intercept);
  }
  model
}

kdmPredict<-function(model,data,predtype,bedType,genomeFile,strict=FALSE,use.float=FALSE){
  if(class(model)!="kdmodel") stop("model must be an objecct of class 'kdmodel'")
  ret<-NULL;
  if(class(data)=="character" & length(data)==1){
    stop("not yet implemented")
  }else if(class(data)=="kbed"){
    ret<-predict_internal_from_bed(model,data,bedType,genomeFile,strict,use.float);
  }else if(class(data)=="character" & length(data)>1){
    ret<-predict_internal_from_sequences(model,data,strict,use.float);
  }else if(class(data)=="kdistr"){
    ret<-predict_internal_from_kdistr(model,data,use.float);
  }else stop("data must be a string vector or a kbed oject or a kdistr object")
  
}

kdmFeatures<-function(motifs, data, strict=FALSE, use.float=FALSE){
  res<-c()
  dclass<-class(data)
  if(dclass=="kdistr"){
    res<-kdfeatures_internal_from_kdistr(motifs,data,use.float)
  }else if(dclass=="character" & length(data)>1){
    res<-kdfeatures_internal_from_seqs(motifs,data,strict,use.float)
  }
  res
}

kdmFeaturesProfile<-function(motifs, sequence, win, strict=FALSE, use.float=FALSE){
  kdfeatures_profile_internal(motifs,sequence,win,strict,use.float)
}


kdmMergeKdmotifs<-function(d1,d2){
  if(class(d1)=="list"){
    R<-matrix(0,nrow(d1[[1]]),sum(unlist(lapply(d1,ncol))))
    cc<-0;
    for(i in 1:length(d1)){
      nc<-ncol(d1[[i]]);
      R[,cc+(1:nc)]<-d1[[i]]@.Data;
      cc<-cc+nc
    }
    res<-new("kdmotifs",R,l=d1[[1]]@l,k=d1[[1]]@k,doubleStrand=d1[[1]]@doubleStrand,estimated=d1[[1]]@estimated)
  }else{
    res<-d1;
    if(!inherits(d1,"kdmotifs") ){
      stop("Objects d1 and d2 must both be either kdistr or kdmotifs")
    }else if( d1@l==d2@l & d1@k==d2@k & d1@doubleStrand==d2@doubleStrand ){
      res@.Data<-cbind(d1,d2)
    }else{
      stop("d1 and d2 are incopatible (i.e. differnt kmer length, differnt gap number or different strandness")
    }
  }
  res
}



kdmGetW0<-function(counts,clustering.method="complete",good.clusters=TRUE,min.rank,max.rank,bhattacharyya=FALSE,use.float=FALSE){
  kdmGetW0_internal(counts,clustering.method,good.clusters,min.rank,max.rank,bhattacharyya,use.float);
}


kdmRankInfo<-function(distance,clustering.method="complete",use.float=FALSE){
  kdmRankInfo_internal(distance,clustering.method,use.float)
}

kdmSilhouette<-function(distance,cluster.labels,use.float=FALSE){
  kdmSilhouette_internal(distance,cluster.labels,use.float)
}

kdmClusterMotifs<-function(motifs,labels,use.float=FALSE){
  rank<-max(labels)
  mmot<-kdmotifs_merge_internal(motifs,matrix(labels-1,length(labels),1),rank,use.float)
  mmot
}

kdmMergeKdistr<-function(d1,d2,use.float=FALSE){
}

kdmSimilarity<-function(d1,d2,use.float=FALSE){
  res<-NULL;
  if((class(d1)!=class(d2))|(class(d1)!="kdistr" & class(d1)!="kdmotifs")){
    error("Objects d1 and d2 must both be either kdistr or kdmotifs")
  }
  if(class(d1)=="kdistr"){
    res<-kdmSimilarity_kdistr_internal(d1,d2,use.float)
  }else if(class(d2)=="kdmotifs"){
    res<-kdmSimilarity_kdmotifs_internal(d1,d2,use.float)
  }
  res
}

kdmDistance<-function(d1,d2,bhattacharyya=FALSE,use.float=FALSE){
  res<-NULL;
  if(missing(d2)){
    if((class(d1)!="kdistr" & class(d1)!="kdmotifs")){    
      error("Objects d1 must be either kdistr or kdmotifs")
    }
    if(class(d1)=="kdistr"){
      res<-kdmDistance_kdistr_self_internal(d1,bhattacharyya,use.float)
    }else if(class(d1)=="kdmotifs"){
      res<-kdmDistance_kdmotifs_self_internal(d1,bhattacharyya,use.float)
    }
    
  }else{
    if((class(d1)!=class(d2))|(class(d1)!="kdistr" & class(d1)!="kdmotifs")){
      error("Objects d1 and d2 must both be either kdistr or kdmotifs")
    }
    if(class(d1)=="kdistr"){
      res<-kdmDistance_kdistr_internal(d1,d2,bhattacharyya,use.float)
    }else if(class(d2)=="kdmotifs"){
      res<-kdmDistance_kdmotifs_internal(d1,d2,bhattacharyya,use.float)
    }
  }
  res
}

kdmModelsCorrelation<-function(a,b=NULL,niter=0){
  if(is.null(b)){
    b=new("kdmodel")
  }
  if(class(a)=="list"){
    nmod<-length(a)
    res<-list(correlation=diag(nmod),p.values=diag(nmod))
    for(i in 1:(nmod-1)){
      for(j in (i+1):nmod){
        tmp<-kdmModelCorrelation_internal(a[[i]],a[[j]],niter);
        res$correlation[i,j]<-res$correlation[i,j]<-tmp$correlation
        res$p.values[i,j]<-res$p.values[i,j]<-tmp$p.values
        cat(i,j,res$correlation[i,j],"\n")
      }
    }
  }else{
    res<-kdmModelCorrelation_internal(a,b,niter);
  }
  res
}

kdmModelsCorrelation2<-function(a,b=NULL,niter=0){
  if(is.null(b)){
    b=new("kdmodel")
  }
  if(class(a)=="list"){
    nmod<-length(a)
    res<-list(correlation=diag(nmod),p.values=diag(nmod))
    for(i in 1:(nmod-1)){
      for(j in (i+1):nmod){
        tmp<-kdmModelCorrelation2_internal(a[[i]],a[[j]],niter);
        res$correlation[i,j]<-res$correlation[i,j]<-tmp$correlation
        res$p.values[i,j]<-res$p.values[i,j]<-tmp$p.values
        cat(i,j,res$correlation[i,j],"\n")
      }
    }
  }else{
    res<-kdmModelCorrelation2_internal(a,b,niter);
  }
  res
}



kdmFeaturesFromPWMSet<-function(sequences,pwmSet,double.strand=TRUE,mode="max"){
  kdmFeaturesFromPWMSet_internal(sequences,pwmSet,double.strand,mode);
}

kdmFeaturesProfileFromPWMSet<-function(sequence,pwmSet){
  kdmFeaturesProfileFromPWMSet_internal(sequence,pwmSet);
}





kdmMotifProfiles2<-function(data,centers,labels,motifs,half.interval,half.win,genome.file,strict=FALSE,use.float=FALSE){
  if(length(centers)==1){
    RES<-kdmMotifProfiles2_internal(data,rep(centers[1],nrow(data)),labels,motifs,half.interval,half.win,genome.file,strict,use.float)
  }else{
    RES<-kdmMotifProfiles2_internal(data,centers,labels,motifs,half.interval,hfal.win,genome.file,strict,use.float)
  }
  RES
}

kdmMotifProfiles3<-function(data,centers,labels,motifs,half.interval,half.win,genome.file,strict=FALSE,use.float=FALSE){
  if(length(centers)==1){
    RES<-kdmMotifProfiles3_internal(data,rep(centers[1],nrow(data)),labels,motifs,half.interval,half.win,genome.file,strict,use.float)
  }else{
    RES<-kdmMotifProfiles3_internal(data,centers,labels,motifs,half.interval,half.win,genome.file,strict,use.float)
  }
  RES
}

kdmMotifProfiles3_seq<-function(pseq,nseq,motifs,half.interval,half.win,strict=FALSE,use.float=FALSE){
  RES<-kdmMotifProfiles3_seq_internal(pseq,nseq,motifs,half.interval,half.win,strict,use.float)
  RES
}


kdmGetProfilesParameters<-function(profile,half.interval,half.win,nrep=100,max.selected.peak,half90.interval,use.float=FALSE){
  Am=kdmGetProfilesParameters_internal(profile,half.interval,half.win,nrep,max.selected.peak,half90.interval,use.float)
  colnames(Am)<-c("baseline","delta","peak","width","cor","aic","fmin","rsd","x95","class","selected","max")
  A<-as.data.frame(Am)
  classnames<-c("CENTRALLY ENRICHED","TWO MERGED PEAKS","TWO OVERLAPPING PEAKS","TWO DISTINCT PEAKS","CONSTANTLY ENRICHED","CENTRALLY DEPLETED","CONSTANTLY DEPLETED")    
  A$class<-as.factor(classnames[A$class])
  A$class<-ordered(A$class,levels=classnames)
  A
}

kdmMotifProfiles<-function(data,centers,labels,motifs,interval,win,maxSelectedPeak,half90interval,genome.file,nrep=100,use.float=FALSE){
  if(length(centers)==1){
    RES<-kdmMotifProfiles_internal(data,rep(centers[1],nrow(data)),labels,motifs,interval,win,nrep,maxSelectedPeak,half90interval,genome.file,use.float)
  }else{
    RES<-kdmMotifProfiles_internal(data,centers,labels,motifs,interval,win,nrep,maxSelectedPeak,half90interval,genome.file,use.float)
  }
  colnames(RES$parameters)<-c("baseline","delta","peak","width","cor","aic","fmin","rsd","x95","class","selected","max")
  RES$parameters<-as.data.frame(RES$parameters)
  
  classnames<-c("CENTRALLY ENRICHED","TWO MERGED PEAKS","TWO OVERLAPPING PEAKS","TWO DISTINCT PEAKS","CONSTANTLY ENRICHED","CENTRALLY DEPLETED","CONSTANTLY DEPLETED")  
  RES$parameters$class=as.factor(classnames[RES$parameters$class])
  RES$parameters$class<-ordered(RES$parameters$class,levels=classnames)
  
  RES
}

kdmAnnotateMotifs<-function(motifs,REF,nhits,normalizeSimilarity=FALSE){
  W<-as(motifs,"kdmotifs")
  refMotifs<-REF$W
  ANNOT<-REF$names
  
  S<-kdmSimilarity(W,refMotifs)
  if(normalizeSimilarity){
    SUP<-kdmSimilarity(refMotifs,refMotifs)
    diag(SUP)<-NA
    a<-apply(SUP,1,mean,na.rm=TRUE)
    for(i in 1:ncol(S)){
      S[,i]<-S[,i]/a[i]
    }
  }
  
  TFnames<-unique(unlist(strsplit(unlist(ANNOT),split=",")))
  tt<-matrix(0,ncol(W),length(TFnames))
  gg<-matrix(0,ncol(W),length(TFnames))
  sco<-matrix(0,ncol(W),length(TFnames))
  
  SR<--(t(apply(S,1,rank))-length(ANNOT))+1
  thresh<-max(S[upper.tri(S)])*nhits
  for(i in 1:ncol(W)){
    ss<-which(S[i,]>thresh)
    if(length(ss)>0){
      for(j in 1:length(ss)){
        cc<-match(ANNOT[[ss[j]]],TFnames)
        gg[i,cc]<-gg[i,cc]+SR[i,ss[j]]
        sco[i,cc]<-sco[i,cc]+S[i,ss[j]]
        tt[i,cc]<-tt[i,cc]+1
      }
    }
  }
  s1<-which(tt>0)
  gg[s1]<-gg[s1]/tt[s1]
  sco[s1]<-sco[s1]/tt[s1]
  
  annotations<-list()
  for(i in 1:nrow(tt)){
    aa<-max(tt[i,])
    if(aa>0){
      sel<-which(tt[i,]==max(tt[i,]))
      oo<-sel[order(gg[i,sel],decreasing=FALSE)]
      annotations[[i]]<-data.frame(name=TFnames[oo],num=tt[i,oo],rank=gg[i,oo],score=sco[i,oo])
    }else{
      annotations[[i]]<-data.frame(name=c())
    }
  }
  annotations
}


kdmPlotMotifProfiles<-function(motifs,P,REF,highlights=NULL,nhits=0.5,scale=FALSE,useModel=FALSE,normalizeSimilarity=FALSE,ecolors=c("darkblue","gray50","gold")){
  
  require(ComplexHeatmap)
  require(circlize)
  ismodel = (class(motifs)=="kdmodel")
  
  A<-P$parameters
  Y<-P$profile$Y
  pos<- P$profile$x
  
  sFF<-if(useModel){t(apply(as.matrix(A[,1:4]),1,function(a){a[1]+a[2]*(exp(-((pos-a[3])/a[4])^2) + exp(-((pos+a[3])/a[4])^2))}))}else{Y}
  sFF<-if(scale){t(scale(t(sFF)))}else{log2(sFF)}
  
  annotations<-kdmAnnotateMotifs(motifs,REF,nhits)
  rownames(sFF)<-unlist(lapply(annotations,function(a){paste(a$name,collapse=",")}))
  
  tindex2<-integer(0)
  labels2<-integer(0)
  if(!is.null(highlights) & length(highlights)>0){
    for(o in 1:length(highlights)){
      g<-grep(highlights[o],rownames(sFF))
      present<-which((g %in% tindex2))
      toadd<-which(!(g %in% tindex2))
      tindex2<-c(tindex2,g[toadd])
      labels2<-c(labels2,rep(highlights[o],length(toadd)))
      labels2[present]<-paste(labels2[present],highlights[o],sep=",")
    }
  }
  oo<-order(A$selected,A$class)
  A<-A[oo,]
  tindex2<-match(tindex2,oo)
  rsf<-range(sFF)
  if(rsf[1]<0 & rsf[2]>0){
    cfun<-colorRamp2(c(rsf[1], 0, rsf[2]), ecolors)
  }else if(rsf[2]>0){
    cfun<-colorRamp2(c(0, rsf[2]), ecolors[2:3])
  }else{
    cfun<-colorRamp2(c(rsf[1], 0), ecolors[1:2])
  }
  classnames<-c("CENTRALLY ENRICHED","TWO MERGED PEAKS","TWO OVERLAPPING PEAKS","TWO DISTINCT PEAKS","CONSTANTLY ENRICHED","CENTRALLY DEPLETED","CONSTANTLY DEPLETED")  
  classcolors<-setNames(RColorBrewer::brewer.pal(name = "Dark2", n = 7), classnames)	
  
  if(ismodel){
    row_ha<-rowAnnotation(
      Highlights=anno_mark(at=tindex2,labels=labels2,labels_gp=grid::gpar(fontsize = 8)),
      Model=anno_barplot(motifs@coefficients[oo],width=unit(2,"cm")),
      Baseline=anno_barplot(A$baseline,baseline=1,width=unit(2,"cm"),ylim=range(A$baseline)),
      Peak=anno_barplot(A$peak,width=unit(2,"cm"),ylim=c(0,max(max(A$peak),1))),
      Delta=anno_barplot(A$delta,baseline=0,ylim=c(0,max(max(A$delta),1))),
      Width=anno_barplot(A$width,ylim=c(0,max(max(A$width),1))),
      Cor=anno_barplot(A$cor,ylim=c(0,1))
    )
  }else{
    row_ha<-rowAnnotation(
      #       Highlights=anno_mark(at=tindex2,labels=labels2),
      #       Baseline=anno_barplot(A$baseline,baseline=1,width=unit(2,"cm"),ylim=range(A$baseline)),
      #       Peak=anno_barplot(A$peak,width=unit(2,"cm"),ylim=c(0,max(max(A$peak),1))),
      #       Delta=anno_barplot(A$delta,baseline=0,width=unit(2,"cm"),ylim=range(A$delta)),
      #       Width=anno_barplot(A$width,width=unit(2,"cm"),ylim=range(A$width)),
      #       Cor=anno_barplot(A$cor,width=unit(2,"cm")),
      #       width=unit(11,"cm")
      Highlights=anno_mark(at=tindex2,labels=labels2,labels_gp=grid::gpar(fontsize = 8)),
      Baseline=anno_barplot(A$baseline,baseline=1,ylim=range(A$baseline)),
      Peak=anno_barplot(A$peak,ylim=c(0,max(max(A$peak),1))),
      Delta=anno_barplot(A$delta,baseline=0,ylim=c(0,max(max(A$delta),1e-6))),
      Width=anno_barplot(A$width,ylim=c(0,max(max(A$width),1e-6))),
      Cor=anno_barplot(A$cor,ylim=c(0,1))
    )
  }
  row_ha_left<-rowAnnotation(
    Classification=A$class,
    annotation_legend_param = list(Classification = list(direction = "horizontal",nrow=4)),
    col=list(Classification=classcolors)
  )
  nn<-c("unselected","selected")
  dd<-data.frame(spec=nn[A$selected+1])
  h1<-Heatmap(sFF[oo,],name="Profile",row_split=dd,cluster_row_slices=FALSE,cluster_rows=FALSE,cluster_columns=FALSE,right_annotation=row_ha,left_annotation=row_ha_left,row_names_gp = grid::gpar(fontsize = 8),col=cfun,row_gap=unit(5,"mm"),heatmap_legend_param = list(direction = "horizontal"))
  draw(h1,heatmap_legend_side="bottom",padding = unit(c(10, 10, 10, 10), "mm"))
  h1
}

kdmGetPWM<-function(W,pwm0,tol=1e-4,max.iter=10000,use.float=FALSE){
  if(nrow(pwm0)==4){
    ret<-kdmGetPWM_internal(W,ncol(pwm0),as.numeric(pwm0),tol,max.iter,use.float)
    rownames(ret)<-c("A","C","G","T")
    ret
  }else{
    stop("pwm0 must be a 4 rows matrix")
  }
}

kdmGetProfileInfo<-function(motifs,regions,centers,labels,genomeFile,halfInterval,halfWin,tolerance=1e-10,strict,use_float){
  result<-kdmProfileInfo_internal(motifs,regions,centers,labels,genomeFile,halfInterval,halfWin,tolerance,strict,use_float)  
  result
}

kdmGetProfileInfoSeq<-function(motifs,pseq,nseq,halfInterval,halfWin,tolerance=1e-10,strict,use_float){
  result<-kdmProfileInfoSeq_internal(motifs,pseq,nseq,halfInterval,halfWin,tolerance,strict,use_float)  
  result
}

kdmCentrimo<-function(profileInfo,symmetric,tail,use_float=FALSE){
  
  res1<-kdmCentrimo_internal(profileInfo$pcounts,profileInfo$ncounts,ncol(profileInfo$phits),ncol(profileInfo$nhits),tail,symmetric,use_float);
  res2<-list()
  for(i in 1:length(res1)){
    res2[[i]]<-as.data.frame(res1[[i]])
    names(res2[[i]])<-c("Start","End","max_pos","mean_pos","np_in","np_out","nn_in","nn_out","Centrality","Centrality_p.value","Centrality_q.value","Enrichment","Enrichment_p.value","Enrichment_q.value")
  }
  res2
}

kdmIC<-function(W){
  ics<-rep(0,ncol(W))
  for(i in 1:ncol(W)){
    w<-as.numeric(W[,i]^2)
    
    sel<-which(w>0);
    ics[i]<- -log2(1/nrow(W))+sum(w[sel]*log2(w[sel]))
  }
  ics
}

kdmTomTargetInfo<-function(targetMotifs){
  IC<-function(W){
    ics<-rep(0,ncol(W))
    for(i in 1:ncol(W)){
      w<-as.numeric(W[,i]^2)
      
      sel<-which(w>0);
      ics[i]<- -log2(1/nrow(W))+sum(w[sel]*log2(w[sel]))
    }
    ics
  }
  targetInfo<-list(similarity=kdmSimilarity(targetMotifs,targetMotifs),ic=kdmIC(targetMotifs))
  WPAR<-data.frame(mean=rep(NA,ncol(targetMotifs)),sd=NA,Intercept=NA,Slope=NA)	
  MIC<-matrix(NA,ncol(targetMotifs),ncol(targetMotifs));
  for(i in 1:ncol(targetMotifs)){
    MIC[i,]<-sqrt(targetInfo$ic*targetInfo$ic[i])
    mm<-lm(targetInfo$similarity[i,]~MIC[i,])
    d<-targetInfo$similarity[i,]-predict(mm)
    
    # 		ss<-which(abs(d)<3*(sd(d)))
    # 		WPAR$mean[i]=mean(d[ss]);
    # 		WPAR$sd[i]=sd(d[ss]);
    
    WPAR$mean[i]=mean(d);
    WPAR$sd[i]=sd(d);
    
    WPAR$Intercept[i]=mm$coefficients[1];
    WPAR$Slope[i]=mm$coefficients[2];
  }
  targetInfo$WPAR<-WPAR
  targetInfo$mod<-lm(sd~ic,data=data.frame(ic=targetInfo$ic,sd=targetInfo$WPAR$sd))
  
  inf<-MIC[upper.tri(MIC)]
  sim<-targetInfo$similarity[upper.tri(targetInfo$similarity)]
  
  targetInfo$m2<-list()
  targetInfo$m2$mod<-lm(sim~inf)
  res<-resid(targetInfo$m2$mod)
  targetInfo$m2$mean<-mean(res)
  targetInfo$m2$sd<-sd(res)
  targetInfo
}

kdmTom<-function(queryMotifs, targetMotifs, correction.method="BH", threshold=0.5, select.by="qvalue", best.hit=FALSE, sort.by="similarity",targetInfo=NULL,targetSel=NULL){
  IC<-function(W){
    ics<-rep(0,ncol(W))
    for(i in 1:ncol(W)){
      w<-as.numeric(W[,i]^2)
      
      sel<-which(w>0);
      ics[i]<- -log2(1/nrow(W))+sum(w[sel]*log2(w[sel]))
    }
    ics
  }
  
  queryMotifs<-as(queryMotifs,"kdmotifs")
  targetMotifs<-as(targetMotifs,"kdmotifs")
  if(length(colnames(targetMotifs))!=ncol(targetMotifs)){
    colnames(targetMotifs)<-1:ncol(targetMotifs)
  }
  
  if(length(colnames(queryMotifs))!=ncol(queryMotifs)){
    colnames(queryMotifs)<-1:ncol(queryMotifs)
  }
  
  tinfo<-if(is.null(targetInfo)){
    tinfo<-kdmTomTargetInfo(targetMotifs);
  }else{
    targetInfo
  }
  
  
  icst<-tinfo$ic
  tnames<-colnames(targetMotifs)
  
  SQ<-kdmSimilarity(queryMotifs,targetMotifs)
  icsq<-kdmIC(queryMotifs)
  qnames<-colnames(queryMotifs)
  
  WPAR<-tinfo$WPAR
  
  PV<-QV<-EV<-matrix(NA,nrow(SQ),ncol(SQ))
  for(i in 1:nrow(SQ)){
    mm<-lm(SQ[i,]~sqrt(icst*icsq[i]))
    d<-SQ[i,]-predict(mm)
    PV[i,]=pnorm(d,WPAR$mean,WPAR$sd,lower.tail=FALSE)
    QV[i,]=p.adjust(PV[i,],method=correction.method)
    EV[i,]=PV[i,] * (2*ncol(targetMotifs))
  }
  dd<-c()
  if(!best.hit){
    tot<-0;
    OO<-c()
    if(select.by=="qvalue"){
      OO<-QV;
    }else if(select.by=="pvalue"){
      OO<-PV
    }else if(select.by=="evalue"){
      OO<-EV;
    }else if(select.by=="similarity"){
      OO<-SQ;
    }
    
    for(i in 1:nrow(PV)){
      ss<-which(OO[i,]<threshold)		
      tot<-tot+length(ss)
    }
    if(tot>0){
      dd<-data.frame(query=rep("",tot),target=rep("",tot),p.value=NA,q.value=NA,e.value=NA,similarity=NA)		
      ccount<-1
      for(i in 1:nrow(OO)){
        sel<-which(OO[i,]<threshold)
        if(sort.by=="similarity"){
          oo<-sel[order(SQ[i,sel],decreasing=TRUE)]
        }else if(sort.by=="pvalue"){
          oo<-sel[order(PV[i,sel])]
        }else if(sort.by=="qvalue"){
          oo<-sel[order(QV[i,sel])]
        }else if(sort.by=="evalue"){
          oo<-sel[order(EV[i,sel])]
        }
        if(length(oo)>0){
          dd[ccount:(ccount+length(oo)-1),]<-data.frame(query=rep(qnames[i],length(oo)),target=tnames[oo],p.value=PV[i,oo],q.value=QV[i,oo],e.value=EV[i,oo],similarity=SQ[i,oo]);
        }
        ccount<-ccount+length(oo)
      }
    }else{
      dd<-c()
    }
  }else{
    dd<-data.frame(query=qnames,target="none",p.value=NA,q.value=NA,e.value=NA,similarity=NA)
    for(i in 1:nrow(PV)){
      ss<-ncol(QV)
      if(select.by=="qvalue"){
        ss<-which(QV[i,]<threshold)		
      }else if(select.by=="pvalue"){
        ss<-which(PV[i,]<threshold)		
      }else if(select.by=="evalue"){
        ss<-which(EV[i,]<threshold)		
      }else if(select.by=="similarity"){
        ss<-which(SQ[i,]<threshold)		
      }
      if(length(ss)>0){
        if(sort.by=="similarity"){
          sel<-ss[which.max(SQ[i,ss,drop=FALSE])]
        }else if(sort.by=="qvalue"){
          sel<-ss[which(QV[i,ss]==min(QV[i,ss]))]
          if(length(sel)>1){
            sel<-sel[which.max(SQ[i,sel])]
          }
        }else if(sort.by=="evalue"){
          sel<-ss[which(EV[i,ss]==min(EV[i,ss]))]
          if(length(sel)>1){
            sel<-sel[which.max(SQ[i,sel])]
          }
        }else if(sort.by=="pvalue"){
          sel<-ss[which.min(PV[i,ss])]
        }
        dd$target[i]<-tnames[sel]
        dd$p.value[i]<-PV[i,sel]
        dd$q.value[i]<-QV[i,sel]
        dd$e.value[i]<-EV[i,sel]
        dd$similarity[i]<-SQ[i,sel]
      }
    }
  }
  nn<-colnames(queryMotifs)
  dd$rank<-NA
  for(r in 1:ncol(queryMotifs)){
    s1<-which(dd$query==nn[r])
    if(length(s1)>0){
      dd$rank[s1]<-1:length(s1)
    }
  }
  dd
}











