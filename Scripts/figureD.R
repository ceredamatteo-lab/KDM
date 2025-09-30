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


bdir<-"/tmp/tempres2/"
#system("mkdir /tmp/tempres/")

u=which(exp_info$Experiment=="ENCSR570WLM")

aucs=data.frame(Experiment=exp_info$Experiment,AUC=0)
for(u in 1:nrow(exp_info)){
    print(u)
    train0=DATASET[[u]]$train
    trainLabels=DATASET[[u]]$trainLabels

# SELEZIONE DEL MOTIVO CENTRALE E ARRICCHITO
    train=train0
    train$start=train$start-100
    train$end=train$end+100
    seqs<-kdmGetSequence(train,file.name="/adat/database/2bit/hg38.2bit")
    spos<-seqs[which(trainLabels==1)]
    sneg<-seqs[which(trainLabels==0)]
    cat(file=paste(bdir,"seqs.fa",sep=""),rbind(paste(">S",1:length(seqs),sep=""),seqs),sep="\n")
    cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos),sep=""),spos),sep="\n")
    cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg),sep=""),sneg),sep="\n")

    cmd_centrimo<-paste("centrimo --oc ",bdir,"centrimo/ --norc --neg ",bdir,"nseqs.fa --local ",bdir,"pseqs.fa ",queryfile,sep="")
    system(cmd_centrimo,intern=TRUE,ignore.stderr=TRUE)
    centrimo<-read_centrimo(paste(bdir,"centrimo/centrimo.tsv",sep=""))

    all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05)]
    if(length(all_selected_motifs)==0){
        print("no centrimo central motifs")
        next}
    sel=all_selected_motifs[1]

    # SALVO SOLO LA PWM SELEZIONATA 

    pwms<-kdmLoadPWMSet(queryfile)
    pp=which(do.call(rbind,strsplit(names(pwms)," "))[,1]==sel)
    all_pwms <- read_meme(queryfile)
    sel_pwms <- all_pwms[pp]
    write_meme(sel_pwms, paste0(bdir,"selected.meme"),overwrite = TRUE)
    pwms<-kdmLoadPWMSet(paste0(bdir,"selected.meme"))
    th2=kdmGetPWMThresholds(pwms,0.001)

    
    # SEQUENZE PER I PROFILI DELLE PWM con FIMO
    train2=train0
    ex=105
    train2$start=train2$start-ex-5
    train2$end=train2$end+ex+5
    seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
    spos2<-seqs2[which(trainLabels==1)]
    #sneg2<-seqs2[which(trainLabels==0)]
    cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos2),sep=""),spos2),sep="\n")
    #cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg2),sep=""),sneg2),sep="\n")

    cmd_fimo<-paste("fimo --oc ",bdir,"fimo/ --norc --thresh 2 --max-stored-scores 5000000 --bfile --uniform-- ",
    bdir,"selected.meme ",bdir,"pseqs.fa ",sep="")
    system(cmd_fimo,intern=TRUE,ignore.stderr=TRUE)
    fimo_pos<-read.table(paste(bdir,"fimo/fimo.tsv",sep=""),header=T)%>%mutate(start=start-146,stop=stop-146,Pos=start+5)

    n1=length(unique(fimo_pos$sequence_name))
    n2=length(unique(fimo_pos$Pos))
    if(nrow(fimo_pos)!=n1*n2){
        print("not all fimo scores")
        next
    }


    # PROFILI DELLE KDM A DIVERSA HW
    train=train0
    pp=which(colnames(mcross_W)==sel)
    new_W=mcross_W[,pp]

    pn=kdmGetProfileInfo(motifs = new_W,regions = train,
        centers = rep(36, nrow(DATASET[[u]]$train)),
        labels = trainLabels,
        genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
        halfInterval = 140+5,
        halfWin = 5,
        tolerance = 1e-10,
        strict = FALSE,
        use_float = FALSE
    )

    ref=fimo_pos%>%mutate(Seq=as.numeric(gsub("P","",sequence_name)))%>%arrange(Seq,Pos)%>%mutate(Ref=ifelse(score>=th2,1,0))%>%pull(Ref)

    pred=c(t(pn$pfeat[[1]]))
    aucs$AUC[u]=auc(roc(ref,pred))
}




fimo_pos=fimo_pos%>%filter(score>=th2)


total_pwm_peaks=0
kdm_good_peaks=0
for(ss in 1:length(spos2)){
peaks=fimo_pos%>%filter(sequence_name==paste0("P",ss))%>%pull(Pos)
if(length(peaks)==0){next}
xx = match(peaks, pn$x)
xx = xx[!is.na(xx)]
total_pwm_peaks=total_pwm_peaks+length(peaks)
kdm_good_peaks=kdm_good_peaks+sum(pn$pfeat[[1]][ss,xx]>pn$threshold)
}



