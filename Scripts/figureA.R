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

# SELEZIONE DEL MOTIVO CENTRALE E ARRICCHITO

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

# CONVERSIONE DEL MOTIVO IN PWM

pwms<-kdmLoadPWMSet(queryfile)
pp=which(do.call(rbind,strsplit(names(pwms)," "))[,1]==sel)
all_pwms <- read_meme(queryfile)
sel_pwms <- all_pwms[pp]
write_meme(sel_pwms, paste0(bdir,"selected.meme"),overwrite = TRUE)
pwms<-kdmLoadPWMSet(paste0(bdir,"selected.meme"))

# PROFILI DELLE KDM A DIVERSA HW

pp=which(colnames(mcross_W)==sel)
new_W=mcross_W[,pp]

hws <- c(3, 6, 15, 35)

all_pn <- lapply(hws, function(hw) {
    message("Profilo con hw = ", hw)
  kdmGetProfileInfo(motifs = new_W,regions = train,
    centers = rep(36, nrow(DATASET[[u]]$train)),
    labels = trainLabels,
    genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
    halfInterval = 140,
    halfWin = hw,
    tolerance = 1e-10,
    strict = FALSE,
    use_float = FALSE
  )
})


# SEQUENZE PER I PROFILI DELLE PWM

ex=105
train2=train
train2$start=train2$start-ex-5
train2$end=train2$end+ex+10
seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
spos2<-seqs2[which(trainLabels==1)]
sneg2<-seqs2[which(trainLabels==0)]

# SCELTA DELLA SEQUENZA

lims=30
good=c()
for(ss in 1:length(spos2)){
	p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]
	peaks=which(p1>0)-146+5
	#peaks=peaks[peaks>= -100 & peaks <= 100]
	#if(all(peaks >= -lims & peaks <= lims)==FALSE){next}
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


# A: PROFILO SULLA SINGOLA SEQUENZA

ss=4842
p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]

tmp=list(data.frame(Pos   = seq_along(p1) - 146 + 5,Score = p1,Method = "PWM"))
peaks <- which(p1 > 0) - 146 + 5

tmp=c(tmp, lapply(seq_along(hws), function(k) {
  data.frame(Pos   = all_pn[[k]]$x,Score = all_pn[[k]]$pfeat[[1]][ss, ],Method = paste0("KDM: hw=", hws[k]))
}))

ths=c(0, sapply(all_pn, function(x) x$thresholds))

final=do.call(rbind, tmp) %>%mutate(Method = factor(Method, levels = unique(Method)))

pA=ggplot(final,aes(x=Pos,y=Score))+geom_line()+
scale_x_continuous(limits=c(-140,140),expand=c(0,0),breaks=c(-140,-70,0,70,140))+
    geom_vline(xintercept=peaks,col="red",lwd=0.2)+
    facet_wrap(~Method,scales="free_y",ncol=1)+
    geom_hline(data=data.frame(Method=unique(final$Method),val=ths),aes(yintercept=val),col="blue",lwd=0.2)+
    theme_bw()+ggtitle("Profile on a single sequence")



# B: PROFILO MEDIO SULL'INTERO SET DI SEQUENZE

get_pwm_profile=function(seqs, type) {
  mat=Reduce("+", lapply(seqs, function(s) kdmFeaturesProfileFromPWMSet(s, pwms)))
  mat=mat / length(seqs)
  melt(data.frame(mat)) %>%
    filter(!is.na(value)) %>%
    mutate(Pos = as.numeric(gsub("X", "", variable)) - 146 + 5, type = type, Feature = value) %>%
    select(Pos, Feature, type)
}

tmp <- list(
  get_pwm_profile(spos2, "Bound") %>% mutate(Method = "PWM"),
  get_pwm_profile(sneg2, "Unbound") %>% mutate(Method = "PWM")
)

tmp <- c(tmp, lapply(seq_along(hws), function(k) {
  pn  <- all_pn[[k]]
  pos <- colMeans(pn$pfeat[[1]])
  neg <- colMeans(pn$nfeat[[1]])
  bind_rows(
    data.frame(Pos = pn$x, Feature = pos, type = "Bound"),
    data.frame(Pos = pn$x, Feature = neg, type = "Unbound")
  ) %>% mutate(Method = paste0("KDM: hw=", hws[k]))
}))

final <- bind_rows(tmp)
final=final%>%mutate(Method=factor(Method,levels=unique(final$Method)))

center <- function(s, l) {
  n <- length(s)
  start <- floor((n - l) / 2) + 1
  end <- start + l - 1
  return(s[start:end])
}

ref=subset(final,Method=="PWM"&type=="Bound")$Feature
cor <- do.call(rbind, lapply(1:4, function(k) {
  pn <- all_pn[[k]]
  mm <- max(colMeans(pn$pfeat[[1]]))
  s2 <- subset(final, Method == paste0("KDM: hw=", hws[k]) & type == "Bound")$Feature
  s1 <- center(ref, length(s2))
  t <- cor.test(s1, s2, method = "pearson")
  
  c(Method = paste0("KDM: hw=", hws[k]),
    Measure = "Correlation",
    q.val = round(t$p.value,4),
    Score = unname(t$estimate),
    y = mm)
}))

cor <- as.data.frame(cor, stringsAsFactors = FALSE)
cor$q.val <- as.numeric(cor$q.val)
cor$Score <- as.numeric(cor$Score)
cor$y <- as.numeric(cor$y)

cc<-centrimo[which(centrimo$E.value<0.05 & centrimo$bin_location<=0),]
cc=cc[1,]
xmin=floor(cc$bin_location-cc$bin_width/2)
xmax=floor(cc$bin_location+cc$bin_width/2)
rect_df <- data.frame(
  Method = unique(final$Method), 
  xmin=xmin,
  xmax=xmax,
  ymin = -Inf,
  ymax = Inf
)
            
annot_list = lapply(seq_along(all_pn), function(k) {
  pn = all_pn[[k]]
  mm = max(colMeans(pn$pfeat[[1]]))
  
  ne = kdmCentrimo(pn, symmetric = FALSE, tail = "upper")[[1]]
  ad = (max(ne$Start) - 1) / 2
  ne = transform(ne, Start = Start - ad, End = End - ad)
  ne = ne[ne$Start == xmin & ne$End == xmax, ]
  
  data.frame(Method  = rep(paste0("KDM: hw=", hws[k]), 2),Measure = c("Enrichment", "Centrality"),
  q.val   = c(ne$Enrichment_q.value, ne$Centrality_q.value),Score=c(ne$Enrichment,ne$Centrality),y= mm,stringsAsFactors = FALSE)
})
annot = do.call(rbind, annot_list)
annot=rbind(annot,cor)
annot$Score=round(annot$Score,2)
annot$Method   = factor(annot$Method, levels = unique(final$Method))

pB=ggplot(final, aes(x = Pos, y = Feature, col = type)) +
  geom_line(lwd=1) +
  facet_wrap(~ Method, scales = "free",ncol=1) +
  scale_color_manual(values = c("gold", "grey")) +
  theme_bw()+scale_x_continuous(limits=c(-140,140),expand=c(0,0),breaks=c(-140,-70,0,70,140))+
  ggtitle("Average profile")+
  geom_rect(data = rect_df,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "gold", alpha = 0.2)+
    geom_text(data=annot%>%filter(Measure=="Enrichment"),aes(y=y,label=paste0(Measure," q.val= ",q.val," (",Score,")")),
            x=xmax+20,hjust=0,vjust=1,inherit.aes = FALSE,size=3)+
    geom_text(data=annot%>%filter(Measure=="Centrality"),aes(y=y,label=paste0(Measure," q.val= ",q.val," (",Score,")")),
            x=xmax+20,hjust=0,vjust=3,inherit.aes = FALSE,size=3)+
    geom_text(data=annot%>%filter(Measure=="Correlation"),aes(y=y,label=paste0(Measure," q.val= ",q.val," (",Score,")")),
            x=xmax+20,hjust=0,vjust=5,inherit.aes = FALSE,size=3)
            




tmp=list()
pwm_feat=matrix(nrow=length(spos2),ncol=nchar(spos2[1]))
for(i in 1:length(spos2)){
	pwm_feat[i,]=kdmFeaturesProfileFromPWMSet(spos2[i], pwms)
}
pwm_feat=pwm_feat[,1:286]
count = rep(0, ncol(pwm_feat))
max_col = max.col(pwm_feat, ties.method = "first")
max_val = apply(pwm_feat, 1, max)
valid_idx = which(max_val > 0)
tab = table(max_col[valid_idx])
count[as.integer(names(tab))] = tab
pos_counts <- count

pwm_feat=matrix(nrow=length(sneg2),ncol=nchar(sneg2[1]))
for(i in 1:length(sneg2)){
	pwm_feat[i,]=kdmFeaturesProfileFromPWMSet(sneg2[i], pwms)
}
pwm_feat=pwm_feat[,1:286]
count = rep(0, ncol(pwm_feat))
max_col = max.col(pwm_feat, ties.method = "first")
max_val = apply(pwm_feat, 1, max)
valid_idx = which(max_val > 0)
tab = table(max_col[valid_idx])
count[as.integer(names(tab))] = tab
neg_counts <- count



tmp[[1]]=rbind(
	data.frame(Pos=seq(1,length(pos_counts))-146+5,Counts=pos_counts,type="Bound",Motif=colnames(new_W)[1],Method="PWM"),
	data.frame(Pos=seq(1,length(neg_counts))-146+5,Counts=neg_counts,type="Unbound",Motif=colnames(new_W)[1],Method="PWM"))

for(k in 1:4){
	pn=all_pn[[k]]
	profile=rbind()
	profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$pcounts),type="Bound",Motif=colnames(new_W)[1]))
	profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$ncounts),type="Unbound",Motif=colnames(new_W)[1]))
	tmp[[k+1]]=profile%>%mutate(Motif=factor(Motif,levels=sel),Method=paste0("KDM: hw=",hws[k]))
}
final=do.call(rbind,tmp)
final=final%>%mutate(Method=factor(Method,levels=unique(final$Method)))


pC=ggplot(final,aes(x=Pos,y=Counts,fill=type))+geom_area(alpha=0.5)+
	facet_wrap(~ Method,,ncol=1) +
  scale_fill_manual(values = c("gold", "grey")) +
  theme_bw()+ggtitle("Hits")+
  scale_x_continuous(limits=c(-140,140),expand=c(0,0),breaks=c(-140,-70,0,70,140))
  
  
  
library(ggpubr)
ggarrange(pA,pC,pB,nrow=1,align="h",widths=c(3,4,4))
