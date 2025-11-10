library(rkdMotifs)
library(universalmotif)

# P22 DATA ----
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")

hw=6
half_interval=70
original_length=35
genomeFile<-"/adat/database/2bit/hg38.2bit"

### BEGIN KDM WALL
basedir<-"/adat/Progetti/RBP/MF/W/W_POS_NEG_ALL/70/"
wall_kdm_info<-c()
wtmp<-c()
for(u in 1:nrow(exp_info)){
  print(u)
  modelfile<-paste(basedir,exp_info$Experiment[u],"_train.rds",sep="")
  W<-as(readRDS(modelfile),"kdmotifs")
  colnames(W)<-1:ncol(W)
  pn<-kdmGetProfileInfo(motifs=W,regions=DATASET[[u]]$train,centers=rep(original_length+1,nrow(DATASET[[u]]$train)),
                        labels=DATASET[[u]]$trainLabels,genomeFile=genomeFile,halfInterval=half_interval,halfWin=hw,tolerance=1e-10,strict=FALSE,use_float=FALSE)
  ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
  res<-c()
  for(k in 1:ncol(W)){
    E<-ne[[k]]
    cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
    cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
    #cond<-which((E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
    if(length(cond)>0){
      #ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
      ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),-log10(E$Centrality_q.value[cond]),decreasing=TRUE)],]
      res<-rbind(res,data.frame(motif=k,ss[1,]))
      
    }
  }
  if(length(res)>1){
    #res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
    res<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value)),-log10(res$Centrality_q.value),decreasing=TRUE),]
    mindex<-1:min(5,nrow(res))
    selected_motifs<-res$motif[mindex]
    wall_kdm_info<-rbind(wall_kdm_info,data.frame(motif=selected_motifs,threshold=pn$thresholds[selected_motifs],auroc=pn$aurocs[selected_motifs],centrality=res$Centrality[mindex],enrichment=res$Enrichment[mindex],experiment=exp_info$Experiment[u],target=exp_info$target[u],cell=exp_info$cell[u]))
    wtmp<-cbind(wtmp,W[,selected_motifs])
  }
  #print(wall_kdm_info)
}
wall_kdm<-new("kdmotifs",wtmp,l=as.integer(6),k=as.integer(2),doubleStrand=FALSE,estimated=FALSE)
colnames(wall_kdm)<-1:ncol(wall_kdm)
wall_kdm_tinfo<-kdmTomTargetInfo(wall_kdm)
RBPs<-unique(exp_info$target)
tom1<-kdmTom(wall_kdm,wall_kdm,best.hit=FALSE,threshold=0.05,targetInfo=wall_kdm_tinfo)
for(i in 1:length(RBPs)){
  sel<-which(wall_kdm_info$target==RBPs[i])
  which((tom1$query %in% sel) & (tom1$target %in% sel))
}
#save(file=paste(basedir,"wall_kdm.Rdata",sep=""),wall_kdm,wall_kdm_info,wall_kdm_tinfo)
save(file=paste("RBP/Paper_Figure/","wall_kdm_",hw,".Rdata",sep=""),wall_kdm,wall_kdm_info,wall_kdm_tinfo)

### END   KDM WALL


### BEGIN STREME WALL
basedir<-"/adat/Progetti/RBP/MF/STREME/POS_NEG_ALL/"
pwmall_streme_info<-c()
pwmall_streme<-list()
for(u in 1:nrow(exp_info)){
  print(u)
  load(paste(basedir,exp_info$Experiment[u],"_train.Rdata",sep=""))
  if(res$nmotifs>0){
    npos=exp_info$nTrainPeaks[u]
    nneg=npos
    centrimo=res$centrimo
    all_selected_motifs=c()
    if(nrow(centrimo)>0){
      centrimo$enrichment<-(centrimo$sites_in_bin/npos)/(centrimo$neg_sites_in_bin/nneg)
      for(i in 1:nrow(centrimo)){
        centrimo$enrichment_p.value[i]<-binom.test(x=centrimo$sites_in_bin[i],n=npos,p=centrimo$neg_sites_in_bin[i]/nneg,alternative="gr")$p.value
        centrimo$enrichment_q.value[i]<-centrimo$enrichment_p.value[i] * nrow(centrimo)
      }
      centrimo<-centrimo[order(rank(-log10(centrimo$E.value))+rank(-log10(centrimo$enrichment_q.value)),-log10(centrimo$E.value),decreasing=TRUE),]
      res$centrimo=centrimo	
      all_selected_motifs<-res$centrimo$motif_id[which(res$centrimo$E.value<0.05 & res$centrimo$enrichment_q.value <0.05 & res$centrimo$enrichment >1 & res$centrimo$bin_location<=0)]
    }
    if(length(all_selected_motifs)>0){
      mindex<-1:min(5,length(all_selected_motifs))
      selected_motifs<-all_selected_motifs[mindex]
      pwm_names<-unlist(lapply(strsplit(names(res$pwms),split=" "),function(x){x[1]}))
      ids<-match(selected_motifs,pwm_names)
      for(i in 1:length(ids)){
        pwmall_streme<-c(pwmall_streme,create_motif(res$pwms[[ ids[i] ]]@.Data,name=pwm_names[ids[i]],alphabet="ACGT"))
      }
      pwmall_streme_info<-rbind(pwmall_streme_info,data.frame(name=pwm_names[ids],experiment=exp_info$Experiment[u],target=exp_info$target[u],cell=exp_info$cell[u]))
    }
  }
}
for(i in 1:length(pwmall_streme)){
  pwmall_streme[[i]]@name<-as.character(i)
}
pwmall_streme_info$name<-1:nrow(pwmall_streme_info)
write_meme(motifs=pwmall_streme,paste("RBP/Paper_Figure/","pwmall_streme_new.meme",sep=""),strand="+",overwrite=TRUE)
save(file=paste("RBP/Paper_Figure/","pwmall_streme_new.Rdata",sep=""),pwmall_streme,pwmall_streme_info)
### END   STREME WALL

### BEGIN KDM POS_NEG_ALL 70 CLASSIFICATION VS WALL_TRAIN SELECTION 
wall_kdmTF<-unique(wall_kdm_info$target)
basedir<-"W/W_POS_NEG_ALL/70/"
classification<-data.frame(target=exp_info$target,cell=exp_info$cell,class="No match",targetcell=NA,nmotifs=0)
for(u in 1:nrow(exp_info)){
  print(u)
  if(exp_info$target[u] %in% wall_kdmTF){
    modelfile<-paste(basedir,exp_info$Experiment[u],"_test.rds",sep="")
    W<-as(readRDS(modelfile),"kdmotifs")
    colnames(W)<-1:ncol(W)
    pn<-kdmGetProfileInfo(motifs=W,regions=DATASET[[u]]$test,centers=rep(original_length+1,nrow(DATASET[[u]]$test)),
                          labels=DATASET[[u]]$testLabels,genomeFile=genomeFile,halfInterval=half_interval,halfWin=hw,
                          tolerance=1e-10,strict=FALSE,use_float=FALSE)
    ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
    res<-c()
    for(k in 1:ncol(W)){
      E<-ne[[k]]
      cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
      cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
      #cond<-which((E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
      if(length(cond)>0){
        #ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
        ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),-log10(E$Centrality_q.value[cond]),decreasing=TRUE)],]
        res<-rbind(res,data.frame(motif=k,ss[1,]))
        
      }
    }
    if(length(res)>1){
      #res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
      res<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value)),-log10(res$Centrality_q.value),decreasing=TRUE),]
      selected_motifs<-res$motif[1:min(5,nrow(res))]
      classification$nmotifs[u]<-length(selected_motifs)
      
      tom<-kdmTom(W,wall_kdm,best.hit=FALSE,threshold=0.05,targetInfo=wall_kdm_tinfo)
      tom$target2<-wall_kdm_info$target[as.numeric(tom$target)]
      tom$cell<-wall_kdm_info$cell[as.numeric(tom$target)]
      
      tom<-tom[which( (tom$q.value * length(selected_motifs) < 0.05) & (tom$query %in% selected_motifs) ),]
      ss<-which(tom$target2==exp_info$target[u])
      if(length(ss)>0){
        classification$class[u]<-as.character(cut(min(tom$rank[ss]),c(0,1,3,10,Inf)))
        classification$targetcell[u]<-tom$cell[ss[which.min(tom$rank[ss])]]
      }
    }else{
      classification$class[u]<-"No central motif"
    }
  }else{
    classification$class[u]<-"Not in reference"
  }
  #print(rbind(class_kdmr_train[u,],classification[u,]))
}
class_kdm_wall_train<-classification
### END   KDM POS_NEG_ALL 70 CLASSIFICATION VS MCROSSBASE SELECTION TRAIN		

### BEGIN STREME POS_NEG_ALL 70 CLASSIFICATION VS PWMALL_STREME SELECTION 
basedir<-"STREME/POS_NEG_ALL/"
load(file=paste(basedir,"pwmall_streme_new.Rdata",sep=""))
pwmall_stremeTF<-unique(pwmall_streme_info$target)
classification<-data.frame(target=exp_info$target,cell=exp_info$cell,class="No match",targetcell=NA,nmotifs=0)
for(u in 1:nrow(exp_info)){
  #print(u)
  if(exp_info$target[u] %in% pwmall_stremeTF){
    load(paste(basedir,exp_info$Experiment[u],"_test.Rdata",sep=""))
    if(res$nmotifs>0){
      npos=exp_info$nTrainPeaks[u]
      nneg=npos
      centrimo=res$centrimo
      all_selected_motifs=c()
      if(nrow(centrimo)>0){
        centrimo$enrichment<-(centrimo$sites_in_bin/npos)/(centrimo$neg_sites_in_bin/nneg)
        for(i in 1:nrow(centrimo)){
          centrimo$enrichment_p.value[i]<-binom.test(x=centrimo$sites_in_bin[i],n=npos,p=centrimo$neg_sites_in_bin[i]/nneg,alternative="gr")$p.value
          centrimo$enrichment_q.value[i]<-centrimo$enrichment_p.value[i] * nrow(centrimo)
        }
        centrimo<-centrimo[order(rank(-log10(centrimo$E.value))+rank(-log10(centrimo$enrichment_q.value)),-log10(centrimo$E.value),decreasing=TRUE),]
        res$centrimo=centrimo	
        all_selected_motifs<-res$centrimo$motif_id[which(res$centrimo$E.value<0.05 & res$centrimo$enrichment_q.value <0.05 & res$centrimo$enrichment >1 & res$centrimo$bin_location<=0)]
      }
      if(length(all_selected_motifs)>0){
        system("mkdir /tmp/tempres/")
        mindex<-1:min(5,length(all_selected_motifs))
        selected_motifs<-all_selected_motifs[mindex]
        pwm_names<-unlist(lapply(strsplit(names(res$pwms),split=" "),function(x){x[1]}))
        ids<-match(selected_motifs,pwm_names)
        tmp_pwms<-list()
        for(i in 1:length(ids)){
          tmp_pwms<-c(tmp_pwms,create_motif(res$pwms[[ ids[i] ]]@.Data,name=pwm_names[ids[i]],alphabet="ACGT"))
        }
        write_meme(motifs=tmp_pwms,"/tmp/tempres/exp.meme",strand="+")
        cmd_tomtom<-"tomtom -oc /tmp/tempres/tomtom/ -xalph -thresh 0.05 -norc /tmp/tempres/exp.meme RBP/Paper_Figure/pwmall_streme_new.meme"
        rr<-system(cmd_tomtom,intern=TRUE,ignore.stderr=TRUE)
        tom<-read.delim("/tmp/tempres/tomtom/tomtom.tsv",comment.char="#")
        if(nrow(tom)>0){
          tom$target2<-pwmall_streme_info$target[as.numeric(tom$Target_ID)]
          tom$cell<-pwmall_streme_info$cell[as.numeric(tom$Target_ID)]
          tom$rank<-NA
          un<-unique(tom$Query_ID)
          for(i in 1:length(un)){
            ss<-which(tom$Query_ID==un[i])
            tom$rank[ss]<-1:length(ss)
          }
        }
        mm<-which(tom$Query_ID %in% selected_motifs)
        if(length(mm)>0){
          tom<-tom[mm,]
          tom<-tom[which(tom$q.value * length(selected_motifs) < 0.05),]
          if(nrow(tom)>0){
            ss<-which(tom$target2==exp_info$target[u])
            if(length(ss)>0){
              classification$class[u]<-as.character(cut(min(tom$rank[ss]),c(0,1,3,10,Inf)))
              classification$targetcell[u]<-tom$cell[ss[which.min(tom$rank[ss])]]
            }
          }
        }
        system("rm -Rf /tmp/tempres/")
      }else{
        classification$class[u]<-"No central motif"
      }
    }else{
      stop("non dovevi finire qua")
    }
  }else{
    classification$class[u]<-"Not in reference"
  }
  #print(rbind(class_pwmall_train[u,],classification[u,]))
}
class_streme_pwmall_train<-classification
### END   STREME POS_NEG_ALL 70 CLASSIFICATION VS PWMALL_STREME SELECTION

tt3=cbind(class_kdm_wall_train$class,class_streme_pwmall_train$class)
colnames(tt3)=c("KDM_WALL_KDMtom","STREME_STREMEALL_Tomtom")
tt3=data.frame(tt3)
tt3[tt3=="No match"]="No match to target"
tt3[tt3=="No central motif"]="No central motifs"
tt3[tt3=="Not in reference"]="Target not in reference"
saveRDS(tt3,"RBP/Paper_Figure/P22.rds")





# P22 PLOT ----

if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
  script_folder="Scripts/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
  script_folder="PARTI_FIGURE/"
}
source(paste0(script_folder,"00_plot_functions.R"))
result=readRDS(paste0(data_folder,"P22.rds"))[,2:1]

p3<-plotClassification(result,"KDM",7)

pdf(paste0(plot_folder,"P22_RBP.pdf"),height = 8,width = 10)
ggarrange(plotlist=p3,nrow=1)
dev.off()
