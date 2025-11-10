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

get_center <- function(x, n) {
  sapply(x, function(str) {
    len <- nchar(str)
    start <- floor((len - n) / 2) + 1
    substr(str, start, start + n - 1)
  })
}

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
  as.data.frame() %>% setNames(c("CellLine","RBP","n")) %>% mutate(ID=paste0(CellLine,"_",RBP)) 
  #mutate(Corr = NA,Corr_W=NA)

#half_window_1=6 
#half_window_2=35
half_windows=c(6,35)
half_interval=150 # Distance to extend upstream/downstream from the CL site for profile visualization.
total_length=half_interval*2+1
original_length=35 # Half length of the original bound peaks
expand=half_interval-original_length 


Correlations_list <- list()
info=rbind()
for(j in 1:nrow(all_names_info)){
  print(j)
  name=all_names[j]
  u=which(exp_info$target==all_names_info$RBP[j] & exp_info$cell==all_names_info$CellLine[j])
  if(length(u)<1){next}
  uu=u
 #corrs=c()
 #corrs2=c()
  for(u in uu){
    exp=exp_info$Experiment[u]
    info=rbind(info,c(name,exp))
    train0=DATASET[[u]]$train
    trainLabels=DATASET[[u]]$trainLabels
    N=exp_info$nTrainPeaks[u]
    if(N>500){
        idx_sel=c(1:500,(1+N):(500+N))
        train0=train0[idx_sel,]
        trainLabels=trainLabels[idx_sel]
    }
    idx_pos=which(trainLabels==1)
    seqs<-kdmGetSequence(train0,file.name="/adat/database/2bit/hg38.2bit")
    
    pp=which(colnames(mcross_W)==name)
    length_pwm=ncol(read_meme(queryfile)[pp][[1]]@motif)
    half_pwm=floor(length_pwm/2)
    new_W=mcross_W[,pp]

    pns=list(
    kdmGetProfileInfo(motifs = new_W,regions = train0,
                         centers = rep(original_length+1, nrow(train0)),
                         labels = trainLabels,
                         genomeFile = "/adat/database/2bit/hg38.2bit",
                         halfInterval = half_interval+half_windows[1],
                         halfWin = half_windows[1],
                         tolerance = 1e-10,
                         strict = FALSE,
                         use_float = FALSE),
    kdmGetProfileInfo(motifs = new_W,regions = train0,
                         centers = rep(original_length+1, nrow(train0)),
                         labels = trainLabels,
                         genomeFile = "/adat/database/2bit/hg38.2bit",
                         halfInterval = half_interval+half_windows[2],
                         halfWin = half_windows[2],
                         tolerance = 1e-10,
                         strict = FALSE,
                         use_float = FALSE)
    )
    
    # SCORE DELLA PWM 
    train2=train0
    train2$start=train2$start-expand-half_pwm
    train2$end=train2$end+expand+half_pwm
    seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
    prof<-pwmProfile(seqs2,queryfile,name,pwms[j])
    mid=(ncol(prof$score[[1]])-1)/2
    if(length_pwm%%2==0){prof$score[[1]][,ncol(prof$score[[1]])-half_pwm+1]=NA}
    
    max_pwm=apply(prof$score[[1]],1,max,na.rm=T)
    max_pwm_pos=max_pwm[idx_pos]
    
    prof_pwm=as.vector(t(prof$score[[1]]))
    prof_pwm=prof_pwm[!is.na(prof_pwm)]
    
    prof_pwm_pos=as.vector(t(prof$score[[1]][idx_pos,]))
    prof_pwm_pos=prof_pwm_pos[!is.na(prof_pwm_pos)]
    
    all_correlations=c()
    for(i in 1:2){
    #print(i)
        half_window=half_windows[i]
        pn=pns[[i]]
        
        max_kdm<- c(apply(pn$pfeat[[1]], 1, max),apply(pn$nfeat[[1]], 1, max))
        max_kdm_pos<- apply(pn$pfeat[[1]], 1, max)
        prof_kdm=c(as.vector(t(pn$pfeat[[1]])),as.vector(t(pn$nfeat[[1]])))
        prof_kdm_pos=as.vector(t(pn$pfeat[[1]]))
        
        corr_kdm_pwm_all=cor.test(prof_kdm,prof_pwm)$estimate # 1.Correlazione lungo tutta la seqeuenza
        corr_kdm_pwm_max=cor.test(max_kdm,max_pwm)$estimate # 2. Correlazione tra i massimi 
        corr_kdm_pwm_all_pos=cor.test(prof_kdm_pos,prof_pwm_pos)$estimate # 3.Correlazione lungo tutta la seqeuenza positive
        corr_kdm_pwm_max_pos=cor.test(max_kdm_pos,max_pwm_pos)$estimate # 4. Correlazione tra i massimi positive 
        
        R<-pwmWin(prof$score[[1]],half_window)
        prof_kdm=c(as.vector(t(pn$pfeat[[1]][,(1+half_window):(total_length-half_window)])),
                as.vector(t(pn$nfeat[[1]][,(1+half_window):(total_length-half_window)])))
        prof_pwm_w=as.vector(t(R))
        prof_pwm_w=prof_pwm_w[!is.na(prof_pwm_w)]
        prof_kdm_pos=as.vector(t(pn$pfeat[[1]][,(1+half_window):(total_length-half_window)]))
        prof_pwm_w_pos=as.vector(t(R[idx_pos,]))
        prof_pwm_w_pos=prof_pwm_w_pos[!is.na(prof_pwm_w_pos)]
        
        corr_kdm_pwm_w=cor.test(prof_kdm,prof_pwm_w)$estimate # 5. Correlazione lungo tutta la sequenza con PWM win
        corr_kdm_pwm_w_pos=cor.test(prof_kdm_pos,prof_pwm_w_pos)$estimate # 6. Correlazione lungo tutta la sequenza con PWM win positive
        all_correlations=c(all_correlations,c(corr_kdm_pwm_all,corr_kdm_pwm_max,corr_kdm_pwm_all_pos,corr_kdm_pwm_max_pos,corr_kdm_pwm_w,corr_kdm_pwm_w_pos))
    }
    
    
    for(i in 1:2){
        half_window=half_windows[i]
        xx=(mid-half_window):(mid+half_window)
        max_pwm=apply(prof$score[[1]][,xx],1,max,na.rm=T)
        
        seqs_sub=unname(get_center(seqs,2*half_window+1))
        FF=kdmFeatures(new_W,seqs_sub)
        
        corr_feature_pwm_max=cor.test(FF@.Data,max_pwm)$estimate
        corr_feature_pwm_max_pos=cor.test(FF@.Data[which(trainLabels==1)],max_pwm[which(trainLabels==1)],)$estimate
    
        all_correlations=c(all_correlations,c(corr_feature_pwm_max,corr_feature_pwm_max_pos))
    }
    
    names(all_correlations)=paste0("Cor",1:16)
    Correlations_list[[length(Correlations_list) + 1]]=all_correlations
   #Correlations_list[[length(Correlations_list) + 1]] <- data.frame(
      #Name = name,
      #Experiment = exp,
      #Corr1 = as.numeric(corrs),
      #Corr2 = as.numeric(corrs2),
      #Corr3 = as.numeric(corrs3),
      #Corr4 = as.numeric(corrs4),
      #stringsAsFactors = FALSE
    #)
  }
  #all_names_info$Corr[j]=mean(corrs)
  #all_names_info$Corr_W[j]=mean(corrs2)
}


Correlations <- do.call(rbind, Correlations_list)

saveRDS(Correlations,"RBP/Paper_Figure/p4B_ALL.rds")
#saveRDS(all_names_info,"RBP/Paper_Figure/p4B.rds")

# P4B.PLOT ----
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

min.peaks=50
corr_info=melt(data.frame(readRDS(paste0(data_folder,"p4B_ALL.rds")))[,c(13,15)])
colnames(corr_info)=c("Method","Correlations")
corr_info$Method <- as.character(corr_info$Method)
corr_info$Method[which(corr_info$Method=="Cor13")]="Half win=6"
corr_info$Method[which(corr_info$Method=="Cor15")]="Half win=35"

info=readRDS(paste0(data_folder,"info_corr.rds"))%>%data.frame()%>%rename(PWM=1,Experiment=2)
load("RBP/ENCODE_eCLIP_DATASET.Rdata")

corr=merge(cbind(info,corr_info),exp_info%>%select(Experiment,nTrainPeaks),by="Experiment")%>%
filter(nTrainPeaks>=min.peaks)
o=corr%>%group_by(Method)%>%summarise(med=median(Correlations))%>%arrange(desc(med))%>%pull(Method)

pA=ggplot(corr%>%mutate(Method=factor(Method,levels=o)),aes(x=Correlations,fill=Method,col=Method))+
	geom_histogram(alpha = 0.5, position = "identity", bins = 100, color = NA) +
  	#geom_density(alpha = 0.4) +
  	theme_bw() +
  	theme(panel.grid = element_blank()) +
  	coord_cartesian(xlim = c(-1, 1))

pdf(paste0(plot_folder,"P4B_RBP.pdf"),height = 5,width = 5)
pA
dev.off()


