library(rkdMotifs)
library(universalmotif)
library(ggplot2)
library(reshape2)
library(dplyr)
#options(width=260)
setwd("/adat/Progetti/RBP/MF")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("DATASET.Rdata")
load("mCrossBase.Rdata")

queryfile="/adat/Progetti/RBP/mCrossBase.meme"

read_centrimo<-function(file){
		ctm<-read.delim(file,header=TRUE,stringsAsFactors=FALSE,row.names=NULL,comment.char="#")
		ctm2<-ctm[,1:(ncol(ctm)-1)]
		colnames(ctm2)<-colnames(ctm[,2:ncol(ctm)])
		ctm2
}

center <- function(s, l) {
  n <- length(s)
  start <- floor((n - l) / 2) + 1
  end <- start + l - 1
  return(s[start:end])
}

bdir<-"/tmp/tempres2/"
#system("mkdir /tmp/tempres/")

exp="ENCSR001KKZ"
u=which(exp_info$Experiment==exp)

#=1
total=rbind()
for(u in 1:nrow(exp_info)){
  print(u)
exp=exp_info$Experiment[u]

train0=DATASET[[u]]$train
trainLabels=DATASET[[u]]$trainLabels

# SELEZIONE DEL MOTIVO CENTRALE E ARRICCHITO
train=train0
#train$start=train$start-100
#train$end=train$end+100
seqs<-kdmGetSequence(train,file.name="/adat/database/2bit/hg38.2bit")
spos<-seqs[which(trainLabels==1)]
sneg<-seqs[which(trainLabels==0)]
cat(file=paste(bdir,"seqs.fa",sep=""),rbind(paste(">S",1:length(seqs),sep=""),seqs),sep="\n")
cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos),sep=""),spos),sep="\n")
cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg),sep=""),sneg),sep="\n")

cmd_centrimo<-paste("centrimo --oc ",bdir,"centrimo/ --norc --neg ",bdir,"nseqs.fa --local ",bdir,"pseqs.fa ",queryfile,sep="")
system(cmd_centrimo,intern=TRUE,ignore.stderr=TRUE)
centrimo<-read_centrimo(paste(bdir,"centrimo/centrimo.tsv",sep=""))

#all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05 & centrimo$bin_location<=0)]
all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05)]
if(length(all_selected_motifs)==0){next}
top_centrimo=all_selected_motifs[1:min(length(all_selected_motifs),5)]

pp=which(colnames(mcross_W)%in%top_centrimo)
new_W=mcross_W[,pp]

pn=kdmGetProfileInfo(motifs = new_W,regions = train0,
    centers = rep(36, nrow(DATASET[[u]]$train)),
    labels = trainLabels,
    genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
    halfInterval = 70,
    halfWin = 6,
    tolerance = 1e-10,
    strict = FALSE,
    use_float = FALSE
)

ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
res<-c()
for(k in 1:ncol(new_W)){
   E<-ne[[k]]
   cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
   cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(mcross_W)) < 0.05 & (E$Enrichment_q.value  * ncol(mcross_W))<0.05 & E$Centrality>1 & E$Enrichment>1)
   if(length(cond)>0){
     ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
     res<-rbind(res,data.frame(motif=k,ss[1,]))
   }
}

if(is.null(res)){
  total=rbind(total,data.frame(Experiment=exp,Motif=top_centrimo,Selected=FALSE))
  next
}

total=rbind(total,data.frame(Experiment=exp,Motif=top_centrimo,Selected=top_centrimo%in%colnames(new_W)[res$motif]))
}






selected_kdm=c()
if(!is.null(res)){
    res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
    selected_kdm<-colnames(mcross_W)[res$motif[1:min(5,nrow(res))]]
    result$sel_kdm[u]=length(selected_kdm)
    result$kdm.target[u]=target%in%do.call(rbind,strsplit(selected_kdm,"\\."))[,2]
}

top_centrimo=readRDS("profiles/centrimo_top5.rds")


result=data.frame(Experiment=exp_info$Experiment,target=exp_info$target,
  sel_centrimo=0,sel_kdm=0,centrimo.target=FALSE,kdm.target=FALSE,Intersect=0)


top_kdm=list()
for(u in 1:nrow(exp_info)){
  print(u)
  if(length(top_centrimo[[u]])==0){
    top_kdm[[u]]=NA
    next}
  exp=exp_info$Experiment[u]
  target=exp_info$target[u]

  train0=DATASET[[u]]$train
  trainLabels=DATASET[[u]]$trainLabels
  
  pn=kdmGetProfileInfo(motifs = mcross_W,regions = train0,
    centers = rep(36, nrow(DATASET[[u]]$train)),
    labels = trainLabels,
    genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
    halfInterval = 70,
    halfWin = 6,
    tolerance = 1e-10,
    strict = FALSE,
    use_float = FALSE
  )
  saveRDS(pn,paste0("profiles/",exp,".rds"))
  #pn=readRDS(paste0("profiles/",exp,".rds"))


result_top5=rbind()

for(exp in exp_info$Experiment[77:250]){
  print(exp)
  #sel_centrimo=top_centrimo[[which(names(top_centrimo)==exp)]]
  #if(length(sel_centrimo)==0){next}
  
  u=which(exp_info$Experiment==exp)
  train=DATASET[[u]]$train
  trainLabels=DATASET[[u]]$trainLabels
  seqs<-kdmGetSequence(train,file.name="/adat/database/2bit/hg38.2bit")
  spos<-seqs[which(trainLabels==1)]
  sneg<-seqs[which(trainLabels==0)]
  cat(file=paste(bdir,"seqs.fa",sep=""),rbind(paste(">S",1:length(seqs),sep=""),seqs),sep="\n")
  cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos),sep=""),spos),sep="\n")
  cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg),sep=""),sneg),sep="\n")

  cmd_centrimo<-paste("centrimo --oc ",bdir,"centrimo/ --norc --neg ",bdir,"nseqs.fa --local ",bdir,"pseqs.fa ",queryfile,sep="")
  system(cmd_centrimo,intern=TRUE,ignore.stderr=TRUE)
  centrimo<-read_centrimo(paste(bdir,"centrimo/centrimo.tsv",sep=""))

  #all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05 & centrimo$bin_location<=0)]
  all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05)]
  if(length(all_selected_motifs)==0){next}
  sel_centrimo=all_selected_motifs[1:min(length(all_selected_motifs),5)]


  result_top5=rbind(result_top5,data.frame(Experiment=exp,Method="Centrimo",Motif=sel_centrimo))
  #pn=readRDS(paste0("profiles/",exp,".rds"))
  file_path <- paste0("profiles/", exp, ".rds")

  pn <- tryCatch(
    readRDS(file_path),
    error = function(e) {
      message("Errore su ", file_path, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(pn)){
    pn=kdmGetProfileInfo(motifs = mcross_W,regions = train,
    centers = rep(36, nrow(DATASET[[u]]$train)),
    labels = trainLabels,
    genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
    halfInterval = 70,
    halfWin = 6,
    tolerance = 1e-10,
    strict = FALSE,
    use_float = FALSE
  )
  saveRDS(pn,paste0("profiles/",exp,".rds"))
  }


  ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
  res<-c()
  for(k in 1:ncol(mcross_W)){
    E<-ne[[k]]
    cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
    #cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(mcross_W)) < 0.05 & (E$Enrichment_q.value  * ncol(mcross_W))<0.05 & E$Centrality>1 & E$Enrichment>1)
    cond<-which((E$Centrality_q.value * ncol(mcross_W)) < 0.05 & (E$Enrichment_q.value  * ncol(mcross_W))<0.05 & E$Centrality>1 & E$Enrichment>1)
    if(length(cond)>0){
      #ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
      ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
      res<-rbind(res,data.frame(motif=k,ss[1,]))
        
    }
  }

  #selected_kdm=c()
  if(!is.null(res)){
    #res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
    res<-res[order( rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value)),res$Centrality,decreasing=c(TRUE,TRUE)),]
    selected_kdm<-colnames(mcross_W)[res$motif[1:min(5,nrow(res))]]
    result_top5=rbind(result_top5,data.frame(Experiment=exp,Method="KDM.Centrimo",Motif=selected_kdm))
    #top_kdm[[u]]=selected_kdm
    #result$sel_kdm[u]=length(selected_kdm)
    #result$kdm.target[u]=target%in%do.call(rbind,strsplit(selected_kdm,"\\."))[,2]
  }
}



  if(is.null(res)){top_kdm[[u]]=NA}
  #t=which(names(top_centrimo)==exp)
  #sel_centrimo=top_centrimo[[t]] 
  #if(length(sel_centrimo)>0){
    #result$sel_centrimo[u]=length(sel_centrimo)
    #result$centrimo.target[u]=target%in%do.call(rbind,strsplit(sel_centrimo,"\\."))[,2]
  #}

  #if(length(selected_kdm)>0&length(sel_centrimo)>0){result$Intersect[u]=sum(selected_kdm%in%sel_centrimo)}
}
