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



bdir<-"/tmp/tempres/"
#system("mkdir /tmp/tempres/")

u=which(exp_info$Experiment=="ENCSR570WLM")
#u=which(exp_info$Experiment=="ENCSR366YOG")

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

all_selected_motifs<-centrimo$motif_id[which(centrimo$E.value<0.05 & centrimo$bin_location<=0)]
sel=all_selected_motifs[1]
#selected_motifs<-all_selected_motifs[1:min(5,length(all_selected_motifs))]



pwms<-kdmLoadPWMSet(queryfile)
pp=which(do.call(rbind,strsplit(names(pwms)," "))[,1]==sel)
all_pwms <- read_meme(queryfile)
sel_pwms <- all_pwms[pp]
write_meme(sel_pwms, paste0(bdir,"selected.meme"),overwrite = TRUE)
pwms<-kdmLoadPWMSet(paste0(bdir,"selected.meme"))

ex=105
train2=train
train2$start=train2$start-ex
train2$end=train2$end+ex+10
seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
spos2<-seqs2[which(trainLabels==1)]
sneg2<-seqs2[which(trainLabels==0)]



pp=which(colnames(mcross_W)==sel)
new_W=mcross_W[,pp]

all_pn=list()
hws=c(1,5,15,35)
for(hw in hws){
pn<-kdmGetProfileInfo(motifs=new_W,regions=train,centers=rep(35,nrow(DATASET[[u]]$train)),labels=trainLabels,genomeFile="/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
	halfInterval=70,halfWin=hw,tolerance=1e-10,strict=FALSE,use_float=FALSE)
all_pn[[length(all_pn)+1]]=pn
}

lims=20
good=c()
for(ss in 1:length(spos2)){
	p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]
	peaks=which(p1>0)-141+5
	peaks=peaks[peaks>= -70 & peaks <= 70]
	if(all(peaks >= -lims & peaks <= lims)==FALSE){next}
	#peaks=peaks[peaks >= -lims & peaks <= lims]
	if(length(peaks)>5 | length(peaks)==0){next}

	keep=TRUE
	for(j in 1:4){
		pn=all_pn[[j]]
		th=pn$thresholds
		xx=which(pn$x%in%peaks)
		if(sum(pn$pfeat[[1]][ss,xx]<th)!=0){keep=FALSE}
	}
	if(keep==TRUE){good=c(good,ss)}
}

ss=4288
tmp=list()
ths=c(0)
p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]
tmp[[1]]=data.frame(Pos=seq(1,length(p1))-141+5,Score=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281],Method="PWM")
peaks=which(p1>0)-141+5
#peaks=peaks[peaks >= -lims & peaks <= lims]
for(k in 1:4){
	tmp[[k+1]]=data.frame(Pos=all_pn[[k]]$x,Score=all_pn[[k]]$pfeat[[1]][ss,],Method=paste0("KDM: hw=",hws[k]))
	ths=c(ths,all_pn[[k]]$thresholds)
}
final=do.call(rbind,tmp)%>%mutate(Method=factor(Method,levels=unique(final$Method)))
a=ggplot(final,aes(x=Pos,y=Score))+geom_line()+xlim(c(-69,+69))+geom_vline(xintercept=peaks,col="red")+
facet_wrap(~Method,scales="free",ncol=1)+
geom_hline(data=data.frame(Method=unique(final$Method),val=ths),aes(yintercept=val),col="blue")+
theme_bw()+ggtitle("Profile on a single sequence")


tmp=list()
pwm_feat=matrix(nrow=length(pwms),ncol=nchar(spos2[1]))
pwm_feat <- Reduce("+", lapply(spos2, function(s) {
  kdmFeaturesProfileFromPWMSet(s, pwms)
}))
pwm_feat=pwm_feat/length(spos)
pwm_feat=melt(data.frame(pwm_feat)%>%mutate(Motif=do.call(rbind,strsplit(names(pwms)," "))[,1]))%>%filter(!is.na(value))%>%
mutate(Pos=as.numeric(gsub("X","",variable))-141,type="Bound")%>%rename(Feature=value)%>%select(Pos,Feature,type,Motif)
tmp[[1]]=pwm_feat%>%mutate(Pos=Pos+5,Method="PWM")

pwm_feat_neg=matrix(nrow=length(pwms),ncol=nchar(sneg2[1]))
pwm_feat_neg <- Reduce("+", lapply(sneg2, function(s) {
  kdmFeaturesProfileFromPWMSet(s, pwms)
}))
pwm_feat_neg=pwm_feat_neg/length(sneg)
pwm_feat_neg=melt(data.frame(pwm_feat_neg)%>%mutate(Motif=do.call(rbind,strsplit(names(pwms)," "))[,1]))%>%filter(!is.na(value))%>%
mutate(Pos=as.numeric(gsub("X","",variable))-141,type="Unbound")%>%rename(Feature=value)%>%select(Pos,Feature,type,Motif)
tmp[[2]]=pwm_feat_neg%>%mutate(Pos=Pos+5,Method="PWM")

for(k in 1:4){
	pn=all_pn[[k]]
	profile=rbind()
	pos=colMeans(pn$pfeat[[1]])
	neg=colMeans(pn$nfeat[[1]])

	profile=rbind(profile,data.frame(Pos=pn$x,Feature=pos,type="Bound",Motif=colnames(new_W)[1]))
	profile=rbind(profile,data.frame(Pos=pn$x,Feature=neg,type="Unbound",Motif=colnames(new_W)[1]))
	profile=profile%>%mutate(Motif=factor(Motif,levels=sel),Method=paste0("KDM: hw=",hws[k]))
	tmp[[k+2]]=profile
}
final=do.call(rbind,tmp)%>%mutate(Method=factor(Method,levels=unique(final$Method)))%>%filter(Pos<=64)
b=ggplot(final, aes(x = Pos, y = Feature, col = type)) +
  geom_line() +
  facet_wrap(~ Method, scales = "free",ncol=1) +
  scale_color_manual(values = c("gold", "grey")) +
  theme_bw()+scale_x_continuous(limits = c(-68, 68))+ggtitle("Average profile")



tmp=list()
pwm_feat=matrix(nrow=length(spos2),ncol=nchar(spos2[1]))
for(i in 1:length(spos2)){
	pwm_feat[i,]=kdmFeaturesProfileFromPWMSet(spos2[i], pwms)
}
pos_counts <- colSums(pwm_feat > 0)[1:281]
pwm_feat=matrix(nrow=length(sneg2),ncol=nchar(sneg2[1]))
for(i in 1:length(sneg2)){
	pwm_feat[i,]=kdmFeaturesProfileFromPWMSet(sneg2[i], pwms)
}
neg_counts <- colSums(pwm_feat > 0)[1:281]

tmp[[1]]=rbind(
	data.frame(Pos=seq(1,length(pos_counts))-141+5,Counts=pos_counts,type="Bound",Motif=colnames(new_W)[1],Method="PWM"),
	data.frame(Pos=seq(1,length(neg_counts))-141+5,Counts=neg_counts,type="Unbound",Motif=colnames(new_W)[1],Method="PWM"))

for(k in 1:4){
	pn=all_pn[[k]]
	profile=rbind()
	profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$pcounts),type="Bound",Motif=colnames(new_W)[1]))
	profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$ncounts),type="Unbound",Motif=colnames(new_W)[1]))
	tmp[[k+1]]=profile%>%mutate(Motif=factor(Motif,levels=sel),Method=paste0("KDM: hw=",hws[k]))
}
final=do.call(rbind,tmp)
final=final%>%mutate(Method=factor(Method,levels=unique(final$Method)))%>%filter(abs(Pos)<=64)


c=ggplot(final,aes(x=Pos,y=Counts,fill=type))+geom_area(alpha=0.5)+
	facet_wrap(~ Method, scales = "free",ncol=1) +
  scale_fill_manual(values = c("gold", "grey")) +
  theme_bw()+scale_x_continuous(limits = c(-68, 68))+ggtitle("Hits")


good=c()
for(ss in 1:length(spos2)){
p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]
if(length(which(p1>0))>5 | length(which(p1>0))==0){next}
if(sum(pn$pfeat[[1]][ss,]>=th)==length(which(p1>0))){good=c(good,ss)}
}

ss=sample(good,1)

p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]
peaks=which(p1>0)-141+5
#3p2=kdmFeaturesProfile(motifs=new_W,sequence=spos[ss],win=1)
#xx=seq(1:length(p2))-length(p2)/2

#plot(xx,p2,type="l")

plot(pn$x,pn$pfeat[[1]][ss,],type="l")
abline(v = peaks, col = "red", lty = 2)
abline(h = th, col = "blue", lty = 2)

#ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")

pwms<-kdmLoadPWMSet(queryfile)
pp=which(do.call(rbind,strsplit(names(pwms)," "))[,1]%in%selected_motifs)
all_pwms <- read_meme(queryfile)
sel_pwms <- all_pwms[pp]
write_meme(sel_pwms, paste0(bdir,"selected.meme"),overwrite = TRUE)
pwms<-kdmLoadPWMSet(paste0(bdir,"selected.meme"))

train2=train
train2$start=train2$start-34
train2$end=train2$end+39
seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
spos2<-seqs2[which(trainLabels==1)]
sneg2<-seqs2[which(trainLabels==0)]

pwm_feat=matrix(nrow=length(pwms),ncol=nchar(spos2[1]))
pwm_feat <- Reduce("+", lapply(spos2, function(s) {
  kdmFeaturesProfileFromPWMSet(s, pwms)
}))
pwm_feat=pwm_feat/length(spos)
pwm_feat=melt(data.frame(pwm_feat)%>%mutate(Motif=do.call(rbind,strsplit(names(pwms)," "))[,1]))%>%filter(!is.na(value))%>%
mutate(Pos=as.numeric(gsub("X","",variable))-70,type="Bound")%>%rename(Feature=value)%>%select(Pos,Feature,type,Motif)

pwm_feat_neg=matrix(nrow=length(pwms),ncol=nchar(sneg2[1]))
pwm_feat_neg <- Reduce("+", lapply(sneg2, function(s) {
  kdmFeaturesProfileFromPWMSet(s, pwms)
}))
pwm_feat_neg=pwm_feat_neg/length(sneg)
pwm_feat_neg=melt(data.frame(pwm_feat_neg)%>%mutate(Motif=do.call(rbind,strsplit(names(pwms)," "))[,1]))%>%filter(!is.na(value))%>%
mutate(Pos=as.numeric(gsub("X","",variable))-70,type="Unbound")%>%rename(Feature=value)%>%select(Pos,Feature,type,Motif)


profile=rbind()
for(i in 1:length(selected_motifs)){
pos=colMeans(pn$pfeat[[i]])
neg=colMeans(pn$nfeat[[i]])

profile=rbind(profile,data.frame(Pos=pn$x,Feature=pos,type="Bound",Motif=colnames(new_W)[i]))
profile=rbind(profile,data.frame(Pos=pn$x,Feature=neg,type="Unbound",Motif=colnames(new_W)[i]))
}
profile=profile%>%mutate(Motif=factor(Motif,levels=selected_motifs))


final=do.call(rbind,list(profile%>%mutate(Model="KDM"),pwm_feat%>%mutate(Pos=Pos+5,Model="PWM"),pwm_feat_neg%>%mutate(Pos=Pos+5,Model="PWM")))%>%
mutate(Model=factor(Model,levels=c("PWM","KDM")))%>%filter(Pos<=64)
