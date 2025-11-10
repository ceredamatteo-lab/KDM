# P10.DATA ---- 
library(rkdMotifs)
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")

read_centrimo<-function(file){
  ctm<-read.delim(file,header=TRUE,stringsAsFactors=FALSE,row.names=NULL,comment.char="#")
  ctm2<-ctm[,1:(ncol(ctm)-1)]
  colnames(ctm2)<-colnames(ctm[,2:ncol(ctm)])
  ctm2
  if(nrow(ctm2)==0){
    ctm2[1,]<-NA
    ctm2$enrichment<-NA
    ctm2$enrichment_p.value<-NA
    ctm2<-ctm2[0,]
  }else{
    ctm2$enrichment<-NA
    ctm2$enrichment_p.value<-NA
  }
  ctm2
}

kdm_centrimo_RBP_nopos<-function(pinfo,sig=0.05,mnames){
  tmp<-data.frame(motif=NA,Start=NA,End=NA,max_pos=NA,np_in=NA,np_out=NA,nn_in=NA,nn_out=NA,Centrality=NA,Centrality_p.value=NA,Centrality_q.value=NA,Enrichment=NA,Enrichment_p.value=NA,Enrichment_q.value=NA,bin_location=NA)
  edf<-tmp[which(!is.na(tmp$motif)),]
  nmotifs<-length(pinfo$pfeat)
  
  centrimo<-kdmCentrimo(pinfo, symmetric = FALSE, tail = "upper")
  
  
  res1<-c() ## internal selection basedo on Centrality, external selection based on Centrality and Enrichment, final ordering based on Centrality and enrichment
  res2<-c() ## internal selection basedo on Centrality, external selection based on Centrality and Enrichment, final ordering based on Centrality
  res3<-c() ## internal selection basedo on Centrality, external selection based on Centrality, final ordering based on Centrality
  res4<-c() ## internal selection basedo on Centrality and enrichment, external selection based on Centrality and Enrichment, final ordering based on Centrality and enrichment
  res5<-c() ## internal selection basedo on Centrality and enrichment, external selection based on Centrality and Enrichment, final ordering based on Centrality
  
  ####INTERNAL SELECTION BASED ON CENTRALITY
  res<-c()	
  for(mn in 1:nmotifs){
    E<-centrimo[[mn]]
    E$Start<-E$Start+pinfo$x[1]
    E$End<-E$End+pinfo$x[1]
    E$max_pos<-E$max_pos+pinfo$x[1]
    E$mean_pos<-E$mean_pos+pinfo$x[1]
    E$bin_location<-(E$End+E$Start)/2
    
    E$Centrality_q.value<-p.adjust(E$Centrality_p.value,method="fdr")
    #E$Centrality_q.value<-E$Centrality_p.value * nrow(E)
    
    sel<-which((E$Centrality_q.value * nmotifs) < sig & E$Centrality>1)
    if(length(sel)>0){
      oo<-sel[order(rank(-log10(E$Centrality_q.value[sel])),decreasing=TRUE)]
      res<-rbind(res,data.frame(motif=mn,E[oo[1],]))
    }
  }
  if(!is.null(res)){
    res$Enrichment_q.value<-res$Enrichment_p.value * nrow(res)
    res1<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value),-log10(res$Centrality_q.value)),decreasing=TRUE),]
    res1<-res1[which((res1$Enrichment_q.value * nmotifs)<sig & res1$Enrichment>1),]
    
    res2<-res[order(rank(-log10(res$Centrality_q.value)),decreasing=TRUE),]
    res2<-res2[which((res2$Enrichment_q.value * nmotifs)<sig & res2$Enrichment>1),]
    
    res3<-res[order(rank(-log10(res$Centrality_q.value)),decreasing=TRUE),]
    
  }else{
    res1<-edf
    res2<-edf
    res3<-edf
  }
  
  ####INTERNAL SELECTION BASED ON CENTRALITY AND ENRICHMENT
  res<-c()		
  for(mn in 1:nmotifs){
    E<-centrimo[[mn]]
    E$Start<-E$Start+pinfo$x[1]
    E$End<-E$End+pinfo$x[1]
    E$max_pos<-E$max_pos+pinfo$x[1]
    E$mean_pos<-E$mean_pos+pinfo$x[1]
    E$bin_location<-(E$End+E$Start)/2
    
    E$Centrality_q.value<-p.adjust(E$Centrality_p.value,method="fdr")
    E$Enrichment_q.value<-p.adjust(E$Enrichment_p.value,method="fdr")
    #E$Centrality_q.value<-E$Centrality_p.value * nrow(E)
    
    sel<-which((E$Centrality_q.value * nmotifs) < sig & E$Centrality>1 & (E$Enrichment_q.value * nmotifs) < sig & E$Enrichment >1)
    if(length(sel)>0){
      oo<-sel[order(rank(-log10(E$Centrality_q.value[sel]))+rank(-log10(E$Enrichment_q.value[sel])),-log10(E$Centrality_q.value[sel]),decreasing=TRUE)]
      res<-rbind(res,data.frame(motif=mn,E[oo[1],]))
    }
  }
  if(!is.null(res)){
    res4<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value),-log10(res$Centrality_q.value)),decreasing=TRUE),]
    #res4<-res4[which((res4$Enrichment_q.value * nmotifs)<sig & res4$Enrichment>1),]
    res5<-res[order(rank(-log10(res$Centrality_q.value)),decreasing=TRUE),]
    #res5<-res5[which((res5$Enrichment_q.value * nmotifs)<sig & res5$Enrichment>1),]
  }else{
    res4<-edf
    res5<-edf
  }
  res1$mname<-mnames[res1$motif]
  res2$mname<-mnames[res2$motif]
  res3$mname<-mnames[res3$motif]
  res4$mname<-mnames[res4$motif]
  res5$mname<-mnames[res5$motif]
  list(res1,res2,res3,res4,res5)
}

bdir<-"/tmp/tempresRBP/"
if (!dir.exists(bdir)) {dir.create(bdir, recursive = TRUE)}
genomeFile<-"/adat/database/2bit/hg38.2bit"
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")
queryfile="RBP/mCrossBase.meme"
selexp<-which(exp_info$target %in% mcross_info$target)

original_length=35
breaks<-c(0,5,10000,Inf)
labels<-c("1-5",">5","Not selected")

selcen<-list()
jaccard<-list()
for(i in 1:5){
  selcen[[i]]<-list()
  jaccard[[i]]<-data.frame(target=exp_info$target[selexp],cell=exp_info$cell[selexp],ncentrimo=NA,nkdm=NA,int=NA,union=NA,j=NA,rank_centrimo=NA,rank_kdm=NA)
}

tmp_centrimo=list()
for(u in 1:length(selexp)){
  print(u)
  exp<-selexp[u]
  
  tmat<-mcross_info$name[which(mcross_info$target==exp_info$target[exp])]
  train0=DATASET[[exp]]$train
  trainLabels=DATASET[[exp]]$trainLabels
  
  train=train0
  train$start=train$start-original_length
  train$end=train$end+original_length
  
  seqs<-kdmGetSequence(train,file.name=genomeFile)
  spos<-seqs[which(trainLabels==1)]
  sneg<-seqs[which(trainLabels==0)]
  cat(file=paste(bdir,"seqs.fa",sep=""),rbind(paste(">S",1:length(seqs),sep=""),seqs),sep="\n")
  cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos),sep=""),spos),sep="\n")
  cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg),sep=""),sneg),sep="\n")
  cmd_centrimo<-paste("centrimo --oc ",bdir,"centrimo/ --local --norc --neg ",bdir,"nseqs.fa ",bdir,"pseqs.fa ",queryfile,sep="")
  system(cmd_centrimo,intern=TRUE,ignore.stderr=TRUE)
  tmpcentrimo<-read_centrimo(paste(bdir,"centrimo/centrimo.tsv",sep=""))
  tmpcentrimo$Start=tmpcentrimo$bin_location - tmpcentrimo$bin_width/2
  tmpcentrimo$End=tmpcentrimo$bin_location + tmpcentrimo$bin_width/2
  
  tmp_centrimo[[u]]=tmpcentrimo
}
saveRDS(tmp_centrimo,"kdmCentrimo/ENCODE_eCLIP_centrimo.rds")

tmpcentrimo=readRDS("kdmCentrimo/ENCODE_eCLIP_centrimo.rds")
for(u in 1:length(selexp)){
  exp<-selexp[u]
  trainLabels<-DATASET[[exp]]$trainLabels
  if(nrow(tmpcentrimo[[u]])>0){
    npos<-sum(trainLabels==1)
    nneg<-sum(trainLabels==0)
    for(i in 1:nrow(tmpcentrimo[[u]])){
      tmpcentrimo[[u]]$enrichment[i]<-(tmpcentrimo[[u]]$sites_in_bin[i]/npos)/(tmpcentrimo[[u]]$neg_sites_in_bin[i]/nneg)
      tmpcentrimo[[u]]$enrichment_p.value[i]<-binom.test(x=tmpcentrimo[[u]]$sites_in_bin[i],n=npos,p=tmpcentrimo[[u]]$neg_sites_in_bin[i]/nneg,alternative="gr")$p.value
      tmpcentrimo[[u]]$enrichment_q.value[i]<-tmpcentrimo[[u]]$enrichment_p.value[i] * nrow(tmpcentrimo[[u]])
    }
  }else{
    tmpcentrimo[[u]][1,]<-NA
    tmpcentrimo[[u]]$enrichment<-NA
    tmpcentrimo[[u]]$enrichment_p.value<-NA
    tmpcentrimo[[u]]$enrichment_q.value<-NA
    tmpcentrimo[[u]]<-tmpcentrimo[[u]][0,]
  }
}

for(u in 1:length(selexp)){
  print(u)
  exp<-selexp[u]
  tmat<-mcross_info$name[which(mcross_info$target==exp_info$target[exp])]
  train0=DATASET[[exp]]$train
  trainLabels=DATASET[[exp]]$trainLabels
  
  centrimo<-list()
  #centrimo[[1]]<-tmpcentrimo[[u]][which(tmpcentrimo[[u]]$E.value<0.05 & (tmpcentrimo[[u]]$fisher_adj_pvalue * ncol(h13_W))<0.05),]
  centrimo[[1]]<-tmpcentrimo[[u]][which(tmpcentrimo[[u]]$E.value<0.05 & (tmpcentrimo[[u]]$enrichment_q.value)<0.05 & tmpcentrimo[[u]]$enrichment>1),]
  centrimo[[1]]<-centrimo[[1]][order(rank(-log10(centrimo[[1]]$E.value))+rank(-log10(centrimo[[1]]$enrichment_q.value)),-log10(centrimo[[1]]$E.value),decreasing=TRUE),]
  
  #centrimo[[2]]<-tmpcentrimo[[u]][which(tmpcentrimo[[u]]$E.value<0.05 & (tmpcentrimo[[u]]$fisher_adj_pvalue * ncol(h13_W))<0.05),]
  centrimo[[2]]<-tmpcentrimo[[u]][which(tmpcentrimo[[u]]$E.value<0.05 & (tmpcentrimo[[u]]$enrichment_q.value)<0.05 & tmpcentrimo[[u]]$enrichment>1),]
  centrimo[[2]]<-centrimo[[2]][order(-log10(centrimo[[2]]$E.value),decreasing=TRUE),]
  
  centrimo[[3]]<-tmpcentrimo[[u]][which(tmpcentrimo[[u]]$E.value<0.05),]
  centrimo[[3]]<-centrimo[[3]][order(-log10(centrimo[[3]]$E.value),decreasing=TRUE),]
  
  centrimo[[4]]<-centrimo[[1]]
  centrimo[[5]]<-centrimo[[2]]
  
  
  pinfo<-kdmGetProfileInfo(motifs=mcross_W,train0,centers=rep(36,nrow(train0)),labels = trainLabels,genomeFile = genomeFile,halfInterval=70,halfWin=5,tolerance = 1e-10,strict=FALSE,use_float=TRUE)
  kdm<-kdm_centrimo_RBP_nopos(pinfo=pinfo,mnames=mcross_info$name,sig=0.05)
  
  for(i in 1:5){
    a<-kdm[[i]]$mname
    b<-centrimo[[i]]$motif_id
    selcen[[i]][[u]]<-list(kdm=a,centrimo=b)
    jaccard[[i]]$ncentrimo[u]<-length(b)
    jaccard[[i]]$nkdm[u]<-length(a)
    jaccard[[i]]$int[u]<-length(intersect(a,b))
    jaccard[[i]]$union[u]<-length(union(a,b))
    jaccard[[i]]$j[u]<-jaccard[[i]]$int[u]/jaccard[[i]]$union[u]
    if(length(tmat)>0){
      jaccard[[i]]$rank_centrimo[u]<-min(match(tmat,b),na.rm=TRUE)
      jaccard[[i]]$rank_kdm[u]<-min(match(tmat,a),na.rm=TRUE)
    }
  }
  
}
save(file="RBP/Paper_Figure/p10.Rdata",jaccard,selcen)

# P10.PLOT ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(irr)

if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
}

load(paste0(data_folder,"P10.Rdata"))
breaks=c(0,5,10000,Inf)
labels=c("1-5",">5","Not selected")
plots2=list()
for(i in c(3,2,1,5,4)){
  d=jaccard[[i]]%>%
    mutate(class.centrimo=cut(rank_centrimo,breaks = breaks,labels=labels),
           class.kdm=cut(rank_kdm,breaks = breaks,labels=labels))
  all_combos <- expand.grid(class.centrimo = unique(d$class.centrimo),class.kdm = unique(d$class.kdm))
  tmp <- d %>%
    count(class.centrimo,class.kdm) %>%
    right_join(all_combos, by = c("class.centrimo","class.kdm")) %>%
    mutate(Freq = ifelse(is.na(n), 0, n)) %>%
    group_by(class.centrimo) %>% mutate(tot.Var1 = sum(Freq)) %>% ungroup() %>%
    group_by(class.kdm) %>% mutate(tot.Var2 = sum(Freq)) %>% ungroup() %>%
    mutate(Jaccard = Freq / (tot.Var1 + tot.Var2 - Freq))
  kk=kappa2(data.frame(d$class.centrimo, d$class.kdm))$value
  plots2[[length(plots2)+1]]=ggplot(tmp,aes(x=class.kdm,y=class.centrimo,fill=Jaccard))+
    scale_y_discrete(limits=rev)+
    geom_point(col="black",shape=21,size=9)+
    coord_equal()+
    geom_text(aes(label=Freq))+scale_fill_gradient(low = "white",high = "purple")+theme_bw()+
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1),plot.title = element_text(size=8))+
    ylab("Centrimo")+xlab("KDM.Centrimo")+ggtitle(paste0("K.Cohen coefficient: ",round(kk,3)))
}

pdf(paste0(plot_folder,"P10_RBP.pdf"),height = 5,width = 5)
plots2[[5]]
dev.off()
