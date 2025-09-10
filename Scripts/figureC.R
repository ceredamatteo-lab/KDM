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

u=which(exp_info$Experiment=="ENCSR194HZU")
#u=which(exp_info$Experiment==exp)
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

all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05 & centrimo$bin_location<=0)]
if(length(all_selected_motifs)==0){next}
sel=all_selected_motifs[1]

# SALVO SOLO LA PWM SELEZIONATA 

pwms<-kdmLoadPWMSet(queryfile)
pp=which(do.call(rbind,strsplit(names(pwms)," "))[,1]==sel)
all_pwms <- read_meme(queryfile)
sel_pwms <- all_pwms[pp]
write_meme(sel_pwms, paste0(bdir,"selected.meme"),overwrite = TRUE)
pwms<-kdmLoadPWMSet(paste0(bdir,"selected.meme"))
th2=kdmGetPWMThresholds(pwms,0.001)

# PROFILI DELLE KDM A DIVERSA HW
train=train0
pp=which(colnames(mcross_W)==sel)
new_W=mcross_W[,pp]

hw=6
pn=kdmGetProfileInfo(motifs = new_W,regions = train0,
    centers = rep(36, nrow(DATASET[[u]]$train)),
    labels = trainLabels,
    genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
    halfInterval = 140+hw,
    halfWin = hw,
    tolerance = 1e-10,
    strict = FALSE,
    use_float = FALSE
)


# SEQUENZE PER I PROFILI DELLE PWM con FIMO
train2=train0
ex=105
train2$start=train2$start-ex-5
train2$end=train2$end+ex+5
seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
spos2<-seqs2[which(trainLabels==1)]
sneg2<-seqs2[which(trainLabels==0)]
cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos2),sep=""),spos2),sep="\n")
cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg2),sep=""),sneg2),sep="\n")

cmd_fimo<-paste("fimo --oc ",bdir,"fimo/ --norc --thresh 2 --max-stored-scores 5000000 --bfile --uniform-- ",
  bdir,"selected.meme ",bdir,"pseqs.fa ",sep="")
system(cmd_fimo,intern=TRUE,ignore.stderr=TRUE)
fimo_pos<-read.table(paste(bdir,"fimo/fimo.tsv",sep=""),header=T)%>%mutate(start=start-146,stop=stop-146,Pos=start+5)


final=rbind(fimo_pos%>%group_by(Pos)%>%summarise(Feature=median(score))%>%mutate(type="Bound",Method = "PWM"),
data.frame(Pos = pn$x, Feature = apply(pn$pfeat[[1]], 2, median), type = "Bound",Method="KDM")
)
ggplot(final,aes(x=Pos,y=Feature))+geom_line()+facet_wrap(~Method,scales="free_y",ncol=1)


ref=subset(final,Method=="PWM"&type=="Bound")$Feature
s2 <- subset(final, Method == "KDM" & type == "Bound")$Feature
s1 <- center(ref, length(s2))
t <- cor.test(s1, s2, method = "pearson")

exp_info$Correlation[u]=t$estimate
exp_info$p.value[u]=t$p.value
exp_info$selected[u]=sel
}


