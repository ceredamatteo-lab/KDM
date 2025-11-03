library(rkdMotifs)
library(universalmotif)
library(reshape2)
library(dplyr)

# P4.DATA ----
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")
queryfile="RBP/mCrossBase.meme"

read_centrimo<-function(file){
  ctm<-read.delim(file,header=TRUE,stringsAsFactors=FALSE,row.names=NULL,comment.char="#")
  ctm2<-ctm[,1:(ncol(ctm)-1)]
  colnames(ctm2)<-colnames(ctm[,2:ncol(ctm)])
  ctm2
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

half_window=6 # KDM half_window
half_interval=140 # Distance to extend upstream/downstream from the CL site for profile visualization.
total_length=half_interval*2+1
original_length=35 # Half length of the original bound peaks
expand=half_interval-original_length # quanto devo "aggiungere ad entrambi i lati per arrivare a 140"
half_pwm=5 # Half length of the PWMs

u=which(exp_info$Experiment=="ENCSR570WLM")
train0=DATASET[[u]]$train
trainLabels=DATASET[[u]]$trainLabels

## SELECTION OF CENTRAL AND ENRICHED MOTIF ----
train=train0
train$start=train$start-expand
train$end=train$end+expand
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

## SAVING SELECTED PWM ----

pwms<-kdmLoadPWMSet(queryfile)
pp=which(do.call(rbind,strsplit(names(pwms)," "))[,1]==sel)
all_pwms <- read_meme(queryfile)
sel_pwms <- all_pwms[pp]
write_meme(sel_pwms, paste0(bdir,"selected.meme"),overwrite = TRUE)
pwms<-kdmLoadPWMSet(paste0(bdir,"selected.meme"))
th2=kdmGetPWMThresholds(pwms,0.001)

## KDM PROFILE ----
train=train0
pp=which(colnames(mcross_W)==sel)
new_W=mcross_W[,pp]

hws <- c(half_window)
all_pn <- lapply(hws, function(hw) {
  message("Profilo con hw = ", hw)
  kdmGetProfileInfo(motifs = new_W,regions = train,
                    centers = rep(original_length+1, nrow(DATASET[[u]]$train)),
                    labels = trainLabels,
                    genomeFile = "/adat/Progetti/RBP/RNA_exons/Reference/hg38/hg38.2bit",
                    halfInterval = half_interval+hw,
                    halfWin = hw,
                    tolerance = 1e-10,
                    strict = FALSE,
                    use_float = FALSE
  )
})

## PWM PROFILE USING FIMO SCORE ----
train2=train0
train2$start=train2$start-expand-half_pwm
train2$end=train2$end+expand+half_pwm
seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
spos2<-seqs2[which(trainLabels==1)]
sneg2<-seqs2[which(trainLabels==0)]
cat(file=paste(bdir,"pseqs.fa",sep=""),rbind(paste(">P",1:length(spos2),sep=""),spos2),sep="\n")
cat(file=paste(bdir,"nseqs.fa",sep=""),rbind(paste(">N",1:length(sneg2),sep=""),sneg2),sep="\n")

shift=half_interval+half_pwm+1 # in questo modo la posizione 1 diventa la posizione -145 in modo che la prima finestra in cui viene posizionata la PWM abbia il centro in -140

cmd_fimo<-paste("fimo --oc ",bdir,"fimo/ --norc --thresh 2 --max-stored-scores 5000000 --bfile --uniform-- ",
                bdir,"selected.meme ",bdir,"pseqs.fa ",sep="")
system(cmd_fimo,intern=TRUE,ignore.stderr=TRUE)
fimo_pos<-read.table(paste(bdir,"fimo/fimo.tsv",sep=""),header=T)%>%mutate(start=start-shift,stop=stop-shift,Pos=start+half_pwm)

cmd_fimo<-paste("fimo --oc ",bdir,"fimo/ --norc --thresh 2 --max-stored-scores 5000000 --bfile --uniform-- ",
                bdir,"selected.meme ",bdir,"nseqs.fa ",sep="")
system(cmd_fimo,intern=TRUE,ignore.stderr=TRUE)
fimo_neg<-read.table(paste(bdir,"fimo/fimo.tsv",sep=""),header=T)%>%mutate(start=start-shift,stop=stop-shift,Pos=start+half_pwm)

## CHOOSING THE SEQUENCE ----

# train2=train0
# ex=105
# train2$start=train2$start-ex-5
# train2$end=train2$end+ex+10
# seqs2<-kdmGetSequence(train2,file.name="/adat/database/2bit/hg38.2bit")
# spos2<-seqs2[which(trainLabels==1)]
# lims=50

# good=c()
# for(ss in 1:length(spos2)){
# 	p1=kdmFeaturesProfileFromPWMSet(spos2[ss], pwms)[1:281]
# 	peaks=which(p1>th2)-146+5
# 	#peaks=peaks[peaks>= -100 & peaks <= 100]
# 	#if(all(peaks >= -lims & peaks <= lims)==FALSE){next}
# 	peaks=peaks[peaks >= -lims & peaks <= 0]
# 	#if(length(peaks)>5 | length(peaks)==0){next}
#   if(length(peaks)>5 | length(peaks)<=1){next}

# 	keep=TRUE
#   j=1
# 	#for(j in 1:4){
# 		pn=all_pn[[j]]
# 		th=pn$thresholds
# 		xx=which(pn$x%in%peaks)
# 		if(sum(pn$pfeat[[1]][ss,xx]<th)!=0){keep=FALSE}
# 	#}
# 	if(keep==TRUE){good=c(good,ss)}
# }


## PROFILES ON THE SELECTED SEQUENCE ----
pwms<-kdmLoadPWMSet(queryfile)
prof<-pwmProfile(spos2,queryfile,sel,pwms[pp])
R<-pwmWin(prof$score[[1]],half_window)

ss=284
sel_fimo=fimo_pos%>%filter(sequence_name==paste0("P",ss))

tmp=list(
  data.frame(Pos=sel_fimo$Pos,Score = sel_fimo$score,Method = "PWM"),
  data.frame(Pos=sel_fimo$Pos[order(sel_fimo$Pos)],Score = R[ss,][(half_pwm+1):(ncol(R)-half_pwm)],Method = "PWM + window: 13")
)

peaks=sel_fimo$Pos[which(sel_fimo$p.value<=0.001)]
tmp=c(tmp, lapply(seq_along(hws), function(k) {
  data.frame(Pos   = all_pn[[k]]$x,Score = all_pn[[k]]$pfeat[[1]][ss, ],Method = paste0("KDM: win=", hws[k]*2+1))
}))

#ths=c(NA,NA, sapply(all_pn, function(x) x$thresholds))

final=do.call(rbind, tmp) %>%mutate(Method = factor(Method, levels = unique(Method)))
saveRDS(list(final=final,peaks=peaks),"RBP/Paper_Figure/p4A.rds")

# P4.PLOT ----
library(ggplot2)
if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
}

half_interval=140 
data=readRDS(paste0(data_folder,"p4A.rds"))
final=data$final
peaks=data$peaks

pA=ggplot(final,aes(x=Pos,y=Score))+geom_line()+
  scale_x_continuous(limits=c(-half_interval,half_interval),
                     expand=c(0,0),breaks=c(-half_interval,-half_interval/2,0,half_interval/2,half_interval))+
  geom_vline(xintercept=peaks,col="red",lwd=0.2)+
  facet_wrap(~Method,scales="free_y",ncol=1)+
  #geom_hline(data=data.frame(Method=unique(final$Method),val=ths),aes(yintercept=val),col="blue",lwd=0.2)+
  theme_bw()+xlab("Score on  single sequence")+theme( plot.margin = unit(c(1, 2, 1, 1), "cm"))

pdf(paste0(plot_folder,"P4A_RBP.pdf"),height = 5,width = 8)
pA
dev.off()

# P9.DATA ----

## HITS COUNT FOR EACH POSITION ----

tmp=list()
tmp[[1]]=rbind(
  fimo_pos%>%group_by(sequence_name)%>%slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%select(sequence_name, Pos, p.value)%>%filter(p.value<=0.001)%>%count(Pos)%>%
    rename(Counts=n)%>%mutate(type="Bound",Method="PWM"),
  fimo_neg%>%group_by(sequence_name)%>%slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%select(sequence_name, Pos, p.value)%>%filter(p.value<=0.001)%>%count(Pos)%>%
    rename(Counts=n)%>%mutate(type="Unbound",Method="PWM")
)

k=1
pn=all_pn[[k]]
profile=rbind()
profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$pcounts),type="Bound"))
profile=rbind(profile,data.frame(Pos=pn$x[-c(1, length(pn$x))],Counts=c(pn$ncounts),type="Unbound"))
tmp[[k+1]]=profile%>%mutate(Method=paste0("KDM: win=",hws[k]*2+1))
final=do.call(rbind,tmp)
final=final%>%mutate(Method=factor(Method,levels=unique(final$Method)))

## IDENTIFICATION OF THE ENRICHED AND CENTRAL REGION ----

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

## CENTRALITY AND ENRICHMENT IN THE SELECTED REGION ----

annot_list = lapply(seq_along(all_pn), function(k) {
  pn = all_pn[[k]]
  mm=max(pn$pcounts)
  ne = kdmCentrimo(pn, symmetric = FALSE, tail = "upper")[[1]]
  ad = (max(ne$Start) - 1) / 2
  ne = transform(ne, Start = Start - ad, End = End - ad)
  ne = ne[ne$Start == xmin & ne$End == xmax, ]
  data.frame(Method  = rep(paste0("KDM: win=", hws[k]*2+1), 2),Measure = c("Enrichment", "Centrality"),
             q.val   = c(ne$Enrichment_q.value, ne$Centrality_q.value),Score=c(ne$Enrichment,ne$Centrality),y= mm,stringsAsFactors = FALSE)
})
annot = do.call(rbind, annot_list)
annot$q.val=round(annot$q.val,4)
annot$Score=round(annot$Score,2)
annot$Method   = factor(annot$Method, levels = unique(final$Method))

saveRDS(list(final=final,rect_df=rect_df,annot=annot),"RBP/Paper_Figure/p9.rds")

# P9.PLOT -----
library(ggplot2)
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

half_interval=140
data=readRDS(paste0(data_folder,"p9.rds"))
final=data$final
annot=data$annot
rect_df=data$rect_df
xmax=unique(rect_df$xmax)

pB=ggplot(final,aes(x=Pos,y=Counts,fill=type))+geom_area(alpha=0.5)+
  facet_wrap(~ Method,ncol=1) +
  scale_fill_manual(values = c("gold", "grey")) +
  theme_bw()+#ggtitle("Hits")+
  scale_x_continuous(limits=c(-half_interval,half_interval),
                     expand=c(0,0),breaks=c(-half_interval,-half_interval/2,0,half_interval/2,half_interval))+
  geom_rect(data = rect_df,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "red", alpha = 0.1)+
  geom_text(data=annot%>%filter(Measure=="Enrichment"),aes(y=y,label=paste0(Measure,": ",Score,"(q.val=",q.val,")")),
            x=xmax+20,hjust=0,vjust=1,inherit.aes = FALSE,size=3)+
  geom_text(data=annot%>%filter(Measure=="Centrality"),aes(y=y,label=paste0(Measure,": ",Score,"(q.val=",q.val,")")),
            x=xmax+20,hjust=0,vjust=3,inherit.aes = FALSE,size=3)+xlab("Position")+ylab("Hits count")

pdf(paste0(plot_folder,"P9_RBP.pdf"),height = 5,width = 8)
pB
dev.off()

