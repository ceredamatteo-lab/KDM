library(rkdMotifs)
library(universalmotif)
library(ggplot2)
library(reshape2)
library(dplyr)

# P4B.DATA ----
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")
queryfile="RBP/mCrossBase.meme"

pwmProfile<-function(seqs,refFile,pwmnames,pwms){
  bdir<-tempdir()
  slen<-nchar(seqs[1])
  snames<-paste("S",1:length(seqs),sep="")
  cat(file=paste(bdir,"seqs.fa",sep=""),rbind(paste(">",snames,sep=""),seqs),sep="\n")
  #cmd_fimo<-paste("fimo --oc ",bdir,"/fimo/ --max-strand --thresh ",threshold," --max-stored-scores 5000000 --skip-matched-sequence ",paste("--motif ",pwmnames,sep="",collapse=" ")," --bfile --uniform-- ",refFile," ",bdir,"seqs.fa > ",bdir,"/fimo.tsv",sep="")
  cmd_fimo<-paste("fimo --norc --thresh 2 --max-stored-scores 5000000 --skip-matched-sequence ",paste("--motif ",pwmnames,sep="",collapse=" ")," --bfile --uniform-- ",refFile," ",bdir,"seqs.fa > ",bdir,"/fimo.tsv",sep="")
  system(cmd_fimo,ignore.stderr=TRUE)
  fimo<-read.delim(paste(bdir,"/fimo.tsv",sep=""))
  fimo$pos<-NA
  res<-list()
  res$score<-list()
  res$pvalue<-list()
  for(i in 1:length(pwms)){
    sel<-which(fimo$motif_id==pwmnames[i])
    lenpwm<-ncol(pwms[[i]])   
    res$score[[i]]<-t(matrix(fimo$score[sel],slen-lenpwm+1,length(seqs)))
    res$p.value[[i]]<-t(matrix(fimo$p.value[sel],slen-lenpwm+1,length(seqs)))
    hlenpwm<-floor(lenpwm/2)
    res$score[[i]]<-cbind(matrix(NA,length(seqs),hlenpwm),res$score[[i]],matrix(NA,length(seqs),lenpwm-1-hlenpwm))
    res$p.value[[i]]<-cbind(matrix(NA,length(seqs),hlenpwm),res$p.value[[i]],matrix(NA,length(seqs),lenpwm-1-hlenpwm))
  }
  res
}


pwmWin<-function(profile,hwin){
  win<-2*hwin+1
  res<-matrix(NA,nrow(profile),ncol(profile))
  for(s in 1:nrow(profile)){
    for(i in 1:(ncol(res)-win+1)){
      res[s,i+hwin]<-max(profile[s,i:(i+win-1)])
    }
  }
  res
}


bdir<-"/tmp/tempres/"
if (!dir.exists(bdir)) {dir.create(bdir, recursive = TRUE)}

pwms<-kdmLoadPWMSet(queryfile)
all_names=do.call(rbind,strsplit(names(pwms)," "))[,1]
all_names_info <- all_names %>% strsplit("\\.") %>% do.call(rbind, .) %>% 
  as.data.frame() %>% setNames(c("CellLine","RBP","n")) %>% mutate(Corr = NA,Corr_W=NA)

half_window=6 # KDM half_window
half_interval=140 # Distance to extend upstream/downstream from the CL site for profile visualization.
total_length=half_interval*2+1
original_length=35 # Half length of the original bound peaks
expand=half_interval-original_length 

for(j in 1:nrow(all_names_info)){
  print(j)
  name=all_names[j]
  u=which(exp_info$target==all_names_info$RBP[j] & exp_info$cell==all_names_info$CellLine[j])
  if(length(u)<1){next}
  uu=u
  corrs=c()
  corrs2=c()
  for(u in uu){
    train0=DATASET[[u]]$train
    trainLabels=DATASET[[u]]$trainLabels
    seqs<-kdmGetSequence(train0,file.name="/adat/database/2bit/hg38.2bit")
    
    pp=which(colnames(mcross_W)==name)
    length_pwm=ncol(read_meme(queryfile)[pp][[1]]@motif)
    half_pwm=floor(length_pwm/2)
    new_W=mcross_W[,pp]
    pn=kdmGetProfileInfo(motifs = new_W,regions = train0,
                         centers = rep(original_length+1, nrow(DATASET[[u]]$train)),
                         labels = trainLabels,
                         genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
                         halfInterval = half_interval+half_window,
                         halfWin = half_window,
                         tolerance = 1e-10,
                         strict = FALSE,
                         use_float = FALSE
    )
    
    train2=train0
    train2$start=train2$start-expand-half_pwm
    train2$end=train2$end+expand+half_pwm
    seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
    prof<-pwmProfile(seqs2,queryfile,name,pwms[j])
    #prof<-pwmProfile(seqs2,queryfile,h13_info$name[pp],h13_u[pp])
    if(length_pwm%%2==0){prof$score[[1]][,ncol(prof$score[[1]])-half_pwm+1]=NA}
    
    prof1=c(as.vector(t(pn$pfeat[[1]])),as.vector(t(pn$nfeat[[1]])))
    prof2=as.vector(t(prof$score[[1]]))
    prof2=prof2[!is.na(prof2)]
    corrs=c(corrs,cor.test(prof1,prof2)$estimate)

    R<-pwmWin(prof$score[[1]],half_window)
    prof1=c(as.vector(t(pn$pfeat[[1]][,(1+half_window):(total_length-half_window)])),
            as.vector(t(pn$nfeat[[1]][,(1+half_window):(total_length-half_window)])))
    prof2=as.vector(t(R))
    prof2=prof2[!is.na(prof2)]
    corrs2=c(corrs2,cor.test(prof1,prof2)$estimate)
  }
  all_names_info$Corr[j]=mean(corrs)
  all_names_info$Corr_W[j]=mean(corrs2)
}

saveRDS(all_names_info,"RBP/Paper_Figure/p4B.rds")

# P9.PLOT ----
library(ggplot2)
library(reshape2)
library(dplyr)
if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
}

corr_info=melt(readRDS(paste0(data_folder,"p4B.rds")),id.vars = c("CellLine","RBP","n"))%>%
  filter(!is.na(value))
pA=ggplot(corr_info,aes(x=value,fill=variable,col=variable))+
  geom_histogram(alpha=0.2,position = "identity",bins=100)+
  #geom_density(alpha=0.4)+
  theme_bw()+theme(panel.grid = element_blank())

pdf(paste0(plot_folder,"P4B_RBP.pdf"),height = 5,width = 5)
pA
dev.off()


