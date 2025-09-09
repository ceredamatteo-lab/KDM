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

#u=which(exp_info$Experiment=="ENCSR570WLM")
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

ref=subset(final,Method=="PWM"&type=="Bound")$Feature
s2 <- subset(final, Method == "KDM" & type == "Bound")$Feature
s1 <- center(ref, length(s2))
t <- cor.test(s1, s2, method = "pearson")

exp_info$Correlation[u]=t$estimate
exp_info$p.value[u]=t$p.value
exp_info$selected[u]=sel
}



ggplot(profiles, aes(x = Pos, y = Feature)) +
  geom_line(lwd=1) +
  facet_wrap(~ Method, scales = "free_y",ncol=1) 

## B.1: CORRELAZIONE TRA OGNI CURVA E QUELLA DELLA PWM


cor <- do.call(rbind, lapply(1:4, function(k) {
  pn <- all_pn[[k]]
  mm <- max(colMeans(pn$pfeat[[1]]))
  s2 <- subset(final, Method == paste0("KDM: win=", hws[k]*2+1) & type == "Bound")$Feature
  s1 <- center(ref, length(s2))
  t <- cor.test(s1, s2, method = "pearson")
  
  c(Method = paste0("KDM: win=", hws[k]*2+1),
    Measure = "Correlation",
    q.val = round(t$p.value,4),
    Score = unname(t$estimate),
    y = mm)
}))

cor <- as.data.frame(cor, stringsAsFactors = FALSE)
cor$q.val <- as.numeric(cor$q.val)
cor$Score <- as.numeric(cor$Score)
cor$y <- as.numeric(cor$y)

## B.2: REGIONE CENTRALE E ARRICCHITA

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

## B.3: CENTRALITA E ARRICCHIMENTO NELLA REGIONE INDIVIDUATA
            
annot_list = lapply(seq_along(all_pn), function(k) {
  pn = all_pn[[k]]
  mm = max(colMeans(pn$pfeat[[1]]))
  
  ne = kdmCentrimo(pn, symmetric = FALSE, tail = "upper")[[1]]
  ad = (max(ne$Start) - 1) / 2
  ne = transform(ne, Start = Start - ad, End = End - ad)
  
  #ne=subset(ne,End-Start!=0)
  #res<-ne[order(rank(ne$Centrality)+rank(ne$Enrichment),ne$Centrality,decreasing=c(TRUE,TRUE)),]
  #res=res[1,]
  ne = ne[ne$Start == xmin & ne$End == xmax, ]
  
  #data.frame(Method  = rep(paste0("KDM: win=", hws[k]*2+1), 2),Measure = c("Enrichment", "Centrality"),
  #q.val   = c(ne$Enrichment_q.value, ne$Centrality_q.value),Score=c(ne$Enrichment,ne$Centrality),y= mm,stringsAsFactors = FALSE)
  data.frame(Method  = rep(paste0("KDM: win=", hws[k]*2+1), 2),Measure = c("Enrichment", "Centrality"),
  q.val   = c(ne$Enrichment_q.value, ne$Centrality_q.value),Score=c(ne$Enrichment,ne$Centrality),y= mm,stringsAsFactors = FALSE)
})
annot = do.call(rbind, annot_list)
annot=rbind(annot,cor)
annot$q.val=round(annot$q.val,4)
annot$Score=round(annot$Score,2)
annot$Method   = factor(annot$Method, levels = unique(final$Method))


pB=ggplot(final, aes(x = Pos, y = Feature, col = type)) +
  geom_line(lwd=1) +
  facet_wrap(~ Method, scales = "free_y",ncol=1) +
  scale_color_manual(values = c("gold", "grey")) +
  theme_bw()+scale_x_continuous(limits=c(-140,140),expand=c(0,0),breaks=c(-140,-70,0,70,140))+
  ggtitle("Average profile")+
  geom_rect(data = rect_df,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "red", alpha = 0.1)+
    geom_text(data=annot%>%filter(Measure=="Enrichment"),aes(y=y,label=paste0(Measure,": ",Score,"(q.val=",q.val,")")),
            x=xmax+20,hjust=0,vjust=1,inherit.aes = FALSE,size=3)+
    geom_text(data=annot%>%filter(Measure=="Centrality"),aes(y=y,label=paste0(Measure,": ",Score,"(q.val=",q.val,")")),
            x=xmax+20,hjust=0,vjust=3,inherit.aes = FALSE,size=3)+
    geom_text(data=annot%>%filter(Measure=="Correlation"),aes(y=y,label=paste0(Measure,": ",Score,"(p.val=",q.val,")")),
            x=xmax+20,hjust=0,vjust=5,inherit.aes = FALSE,size=3)
            

# C: HITS POSOTIVE E NEGATIVE

tmp=list()
#compute_pwm_counts <- function(seqs, pwms, th2, max_cols = 286) {
#  pwm_feat <- t(sapply(seqs, function(seq) {
#    kdmFeaturesProfileFromPWMSet(seq, pwms)
#  }))
#  pwm_feat <- pwm_feat[, 1:max_cols, drop = FALSE]
#  count <- rep(0, ncol(pwm_feat))
#  max_col <- max.col(pwm_feat, ties.method = "first")
#  max_val <- apply(pwm_feat, 1, max)
#  valid_idx <- which(max_val > th2)
#  tab <- table(max_col[valid_idx])
#  count[as.integer(names(tab))] <- tab
#  
#  return(count)
#}

#pos_counts <- compute_pwm_counts(spos2, pwms, th2)
#neg_counts <- compute_pwm_counts(sneg2, pwms, th2)

tmp[[1]]=rbind(
fimo_pos%>%group_by(sequence_name)%>%slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%select(sequence_name, Pos, p.value)%>%filter(p.value<=0.001)%>%count(Pos)%>%
  rename(Counts=n)%>%mutate(type="Bound",Method="PWM"),
fimo_neg%>%group_by(sequence_name)%>%slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%select(sequence_name, Pos, p.value)%>%filter(p.value<=0.001)%>%count(Pos)%>%
  rename(Counts=n)%>%mutate(type="Unbound",Method="PWM")
)


#tmp[[1]]=rbind(
#	data.frame(Pos=seq(1,length(pos_counts))-146+5,Counts=pos_counts,type="Bound",Method="PWM"),
#	data.frame(Pos=seq(1,length(neg_counts))-146+5,Counts=neg_counts,type="Unbound",Method="PWM"))

for(k in 1:4){
	pn=all_pn[[k]]
	profile=rbind()
	profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$pcounts),type="Bound"))
	profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$ncounts),type="Unbound"))
	tmp[[k+1]]=profile%>%mutate(Method=paste0("KDM: win=",hws[k]*2+1))
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


