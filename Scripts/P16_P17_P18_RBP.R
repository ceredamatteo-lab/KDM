library(rkdMotifs)
library(reshape2)
library(dplyr)

# P16+P17 DATA ----
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
c25 <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A","#FF7F00", "gold1","skyblue2", "#FB9A99","palegreen2","#CAB2D6",
"#FDBF6F", "gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown", "black"
)
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")

exp="ENCSR534YOI"

half_window=6
original_length=35
half_interval=70
top_kmer=250
th=0.05
genomeFile<-"/adat/database/2bit/hg38.2bit"

u=which(names(DATASET)==exp)
#modelfile<-paste("RBP/Paper_Figure/",exp,"_train.rds",sep="")
modelfile=paste("/adat/Progetti/RBP/MF/W/W_POS_NEG_ALL/70/",exp,"_train.rds",sep="")
W<-as(readRDS(modelfile),"kdmotifs")
colnames(W)<-1:ncol(W)

pn<-kdmGetProfileInfo(motifs=W,regions=DATASET[[u]]$train,centers=rep(original_length,nrow(DATASET[[u]]$train)),labels=DATASET[[u]]$trainLabels,genomeFile=genomeFile,
halfInterval=half_interval,halfWin=half_window,tolerance=1e-10,strict=FALSE,use_float=FALSE)
ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
res<-c()
for(k in 1:ncol(W)){
    E<-ne[[k]]
    cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
    cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(W)) < th & (E$Enrichment_q.value  * ncol(W))<th & E$Centrality>1 & E$Enrichment>1)
    if(length(cond)>0){
        #ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
        ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),-log10(E$Centrality_q.value[cond]),decreasing=TRUE)],]
        res<-rbind(res,data.frame(motif=k,ss[1,]))
    }
}
#res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
res<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value),-log10(res$Centrality_q.value)),decreasing=TRUE),]
MR<-kdmMotifsReduce(W,pn,motifSelection=res$motif,use_negatives=TRUE,weights=NULL,maxdist=5,siglev=1e-4,bonf=TRUE)
W0=W
W<-MR$WR
colnames(W)<-1:ncol(W)
pn<-kdmGetProfileInfo(motifs=W,regions=DATASET[[u]]$train,centers=rep(original_length,nrow(DATASET[[u]]$train)),labels=DATASET[[u]]$trainLabels,genomeFile=genomeFile,
halfInterval=half_interval,halfWin=half_window,tolerance=1e-10,strict=FALSE,use_float=FALSE)
ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
res<-c()
for(k in 1:ncol(W)){
    E<-ne[[k]]
    cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
    cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(W)) < th & (E$Enrichment_q.value  * ncol(W))<th & E$Centrality>1 & E$Enrichment>1)
    if(length(cond)>0){
        #ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
        ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),-log10(E$Centrality_q.value[cond]),decreasing=TRUE)],]
        res<-rbind(res,data.frame(motif=k,ss[1,]))
    }
}
#res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
res<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value),-log10(res$Centrality_q.value)),decreasing=TRUE),]
selected_motifs<-res$motif[1:min(5,nrow(res))]

tom<-kdmTom(W,mcross_W,best.hit=FALSE,threshold=th,targetInfo=mcross_tinfo)
tom$target2<-mcross_info$target[match(tom$target,mcross_info$name)]
tom$cell<-mcross_info$cell[match(tom$target,mcross_info$name)]
tom<-tom[which( (tom$q.value * length(selected_motifs) < th) & (tom$query %in% selected_motifs) ),]
ss<-which(tom$target2==exp_info$target[u] & tom$rank==1)

good=as.numeric(tom$query[ss])
pwm_sel=tom$target[ss]
W<-W0[,MR$clusters[[good]]]
colnames(W)<-1:ncol(W)

tom<-kdmTom(W,mcross_W,best.hit=FALSE,threshold=th,targetInfo=mcross_tinfo)
tom$target2<-mcross_info$target[match(tom$target,mcross_info$name)]
ranks<-rep(Inf,ncol(W))
buoni<-which(tom$target2==exp_info$target[u] & tom$target==pwm_sel)
ranks[as.numeric(tom$query[buoni])]<-tom$rank[buoni]

KDMM<-MR$WR[,good]
mnames<-c("PWM",paste("KDM-",1:ncol(W)," (",ranks,")",sep=""),"KDMM (1)")
cscale<-c25[1:(ncol(W)+2)]
kPWM<-mcross_W[,which(colnames(mcross_W)==pwm_sel)]
probs<-as.numeric(kPWM)^2
oo<-order(probs,decreasing=TRUE)
sel_kmers<-oo[1:top_kmer]
infoPWM<-data.frame(name="PWM",prob=probs[sel_kmers],x=1:length(sel_kmers))
infoKDM<-c()
for(j in 1:ncol(W)){
    pp<-as.numeric(W[,j]^2)[sel_kmers]
    infoKDM<-rbind(infoKDM,data.frame(name=paste("KDM-",j," (",ranks[j],")",sep=""),prob=pp,x=1:length(sel_kmers)))
}
infoKDMM<-data.frame(name="KDMM (1)",prob=as.numeric(KDMM[,1]^2)[sel_kmers],x=1:length(sel_kmers))
info<-rbind(infoPWM[,1:3],infoKDM,infoKDMM)
info$name<-factor(info$name,levels=c("PWM",paste("KDM-",1:ncol(W)," (",ranks,")",sep=""),"KDMM (1)"))
info2<-info %>% group_by(name) %>% summarise("Cumulated probability"=sum(prob))
infoKDM2<-merge(merge(infoKDM,infoKDM %>% group_by(x) %>% summarise(tot=sum(prob)),by="x"),infoKDMM[,2:3],by="x") %>% mutate(prob=prob.x * prob.y / tot)
infoKDM2<-infoKDM2[,c(2,6,1)]
infoKDM2$plt<-"KDMM (1)"
infoPWM$plt<-"PWM"
info2$name<-factor(info2$name,levels=mnames)
info3<-rbind(infoPWM,infoKDM2)
info3$name<-factor(info3$name,levels=mnames)

save(file="RBP/Paper_Figure/P16_P17.Rdata",info2,info3,cscale)


# P16+P17 PLOT ----

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

load(paste0(data_folder,"P16_P17.Rdata"))

p1<-ggplot(data=info2,aes(x=name,y=`Cumulated probability`,fill=name)) + geom_bar(stat="identity",col="black") + scale_fill_manual(values=cscale)+theme_bw()+theme(panel.grid=element_blank())+
geom_text(aes(y=`Cumulated probability`/2,label=paste0(round(`Cumulated probability`,3)*100,"%")))

p2<-ggplot(data=info3,aes(x=x,y=prob,fill=name)) + geom_bar(position="stack",stat="identity") + scale_fill_manual(values=cscale) + facet_grid(plt~.,scale="free_y")+theme_bw()+
theme(panel.grid=element_blank(),strip.text.y=element_text(angle=0))

pdf(paste0(plot_folder,"P16_RBP.pdf"),height = 5,width = 12)
p2
dev.off()


pdf(paste0(plot_folder,"P17_RBP.pdf"),height = 7,width = 7)
p1
dev.off()

# P18 DATA
label=DATASET[[u]]$trainLabels
anno=read.csv(paste0("RBP/Paper_Figure/",exp,".hg38_multianno.csv"))[which(label==1),]
reg=DATASET[[u]]$train[which(label==1),]
tt<-table(anno$Func.refGene)
tt2<-tt[which(tt>40)]
selpeaks<-which(!is.na(anno$Func.refGene) & (anno$Func.refGene %in% names(tt2)))
anno<-anno[selpeaks,]
seq=kdmGetSequence(reg,genomeFile)
FF=kdmFeatures(W,seq)
ff=data.frame(FF[selpeaks,])
colnames(ff)[1:3]=paste0("Motif-",1:3)
ff$Annovar=anno$Func.refGene
mm=merge(melt(ff,id.vars="Annovar"),ff%>%count(Annovar),by="Annovar")%>%
mutate(Annotation=gsub("_","\n",paste0(Annovar,"\n(",n,")")))

class=unique(mm$Annotation)
res <- data.frame(
  variable = character(),
  group1 = character(),
  group2 = character(),
  pval = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(length(class) - 1)) {
  for (j in (i + 1):length(class)) {
    for (m in unique(mm$variable)) {
      s1 <- subset(mm, Annotation == class[i] & variable == m)$value
      s2 <- subset(mm, Annotation == class[j] & variable == m)$value

      p1 <- wilcox.test(s1, s2, alternative = "greater")$p.value
      p2 <- wilcox.test(s2, s1, alternative = "greater")$p.value

      res <- rbind(
        res,
        data.frame(variable = m, group1 = class[i], group2 = class[j], pval = p1),
        data.frame(variable = m, group1 = class[j], group2 = class[i], pval = p2)
      )
    }
  }
}
res$p.adj=p.adjust(res$pval,method="fdr")
res=res%>%mutate(Score=ifelse(p.adj>0.05,0,-log10(p.adj)),Sign=Score!=0)

res_mod <- res %>%
  mutate(
  group1_mod = factor(gsub("_","\n", group1), levels = unique(gsub("_","\n", group1))),
  group2_mod = factor(gsub("_","\n", group2), levels = unique(gsub("_","\n", group1))))


save(file="RBP/Paper_Figure/P18.Rdata",mm,res_mod)


# P18 PLOT ----

library(ggplot2)
library(ggpubr)
if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
}


load(paste0(data_folder,"P18.Rdata"))
pA=ggplot(mm, aes(x = Annotation, y = value, fill = Annovar)) +
  #geom_boxplot(notch = TRUE, outlier.shape = NA,width=0.5) +
  geom_boxplot(notch = TRUE,width=0.5,outlier.size=0.2,show.legend=F)+
  facet_wrap(~variable, nrow = 1, scales = "free_y") +
  scale_fill_brewer(palette = "Set2")+theme_bw()+theme(panel.grid=element_blank())

group1_levels <- levels(res_mod$group1_mod)
group2_levels <- levels(res_mod$group2_mod)  
colors_group1 <- RColorBrewer::brewer.pal(n = length(group1_levels), name = "Set2")
colors_group2 <- RColorBrewer::brewer.pal(n = length(group2_levels), name = "Set2")


plot_list=list()
for(m in paste0("Motif-",1:3)){
    sub=subset(res_mod,variable==m)
    
    p <- ggplot(sub, aes(x=group2_mod, y=group1_mod)) +
        geom_point(aes(fill=Sign), shape=21, col="white", show.legend=F, size=4) +
        scale_fill_manual(values=c("white","black")) +
        theme_bw() +
        theme(panel.grid=element_blank()) +
        coord_equal()
    
    for(i in seq_along(group2_levels)){
        p <- p + annotate("rect",xmin = -0.5, xmax = 0, 
        ymin = i - 0.5, ymax = i + 0.5, fill = colors_group2[i],col="black")
    }
    
    for(i in seq_along(group1_levels)){
        p <- p + annotate("rect",ymin = -0.5, ymax = 0, 
        xmin = i - 0.5, xmax = i + 0.5,fill = colors_group1[i],col="black")
    }

    plot_list[[length(plot_list)+1]]=p+theme(axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank())
    
    
}


pdf(paste0(plot_folder,"P18_RBP.pdf"),height = 5,width = 10)
ggarrange(pA,ggarrange(plotlist=plot_list,nrow=1),ncol=1,heights=c(7,3))
dev.off()
  
  

