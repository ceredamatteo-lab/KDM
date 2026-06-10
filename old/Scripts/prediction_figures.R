library(glmnet)
library(ROCR)
library(rkdMotifs)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)

basedir<-"/adat/Progetti/KDM/MATERIALE/ALL"

makePlot<-function(data,ntrain,nq=4,yname="Accuracy",toconsider,pval=0.01,npeaks=NULL,alternative="two",avg="mean"){
	pal <- c(
	"dodgerblue2", "#E31A1C", # red
	"green4",
	"#6A3D9A", # purple
	"#FF7F00", # orange
	"gold1",
	"skyblue2", "#FB9A99", # lt pink
	"palegreen2",
	"#CAB2D6", # lt purple
	"#FDBF6F", # lt orange
	"gray70", "khaki2",
	"maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
	"darkturquoise", "green1", "yellow4", "yellow3",
	"darkorange4", "black", "brown"
	)

	if(avg=="mean"){
		avgfun<-mean
	}else{
		avgfun<-median
	}
	apal<-pal[1:length(toconsider)]	
	data<-data[,match(toconsider,names(data))]
	mranks<-apply(apply(data,1,rank),1,avgfun)
	or<-order(mranks)
	data<-data[,or]
	mranks<-mranks[or]
	onames<-names(data)
	names1<-rev(onames)[1:(length(onames)-1)]
	names2<-onames[1:(length(onames)-1)]
	
	
	d1<-melt(cbind(data,Type="DNN based"))
	d1$Type[which(d1$variable %in% c("ls-GKM","KDM-LRLM","SeqGL"))]<-"k-mer based"	
	d1$meanRank<-NA
	for(i in 1:length(onames)){
		sel<-which(d1$variable==onames[i])
		d1$meanRank[sel]<-mranks[i]
	}
	names(d1)<-c("Model type","Method",yname,"Mean rank")
	d1$Method<-factor(d1$Method,levels=onames)
	
	p1<-ggplot(data=d1,aes(x=Method,y=.data[[yname]],fill=`Model type`)) + geom_boxplot(outlier.shape=NA) + scale_fill_manual(values=c("NA",pal[9])) + guides(fill = guide_legend(position="inside")) + theme_bw(base_size=8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))     

	rup<-c()
	for(i in 1:length(names1)){
		for(j in 1:(length(names2)-i+1)){
			x<-data[,match(names1[i],onames)]
			y<-data[,match(names2[j],onames)]
			pptest<-wilcox.test(x,y,paired=TRUE,alternative=alternative)
			pp<-pptest$p.value
			if(is.na(pp)){
				pp<-1
			}
			dd<-(median(x-y))
			rup<-rbind(rup,
			data.frame(
				method1=names1[i],
				method2=names2[j],
				difference=dd,
				p.value=pp
			))
		}
	}
	rup$adjusted.p.value<-p.adjust(rup$p.value,method="bonferroni")
	rup$label<-formatC(rup$adjusted.p.value, format = "e", digits = 2)
	rup$label[which(rup$adjusted.p.value>=pval)]<-"n.s."
# 	for(i in 1:nrow(rup)){
# 		if(rup$adjusted.p.value[i]<pval){
# 			if(rup$difference[i]>0){
# 				rup$label[i]<-"+"
# 			}else{
# 				rup$label[i]<-"-"
# 			}
# 		}else{
# 			rup$label[i]<-""
# 		}
# 	}
  	rup$method1<-factor(rup$method1,levels=rev(names1))
  	rup$method2<-factor(rup$method2,levels=names2)
  	dname<-paste("Median(D",yname,")",sep="")
  	dname2<-sym(dname)
  	names(rup)<-c("Method1","Method2",dname,"p.value","adjusted.p.value","label")
	p2<-ggplot(data=rup,aes(y=Method1,x=Method2,fill=!!dname2,label=label)) + geom_tile(color="black") + geom_text(size=3) + scale_fill_gradient2(low="blue",high="red",mid="white") + theme_bw(base_size=8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	if(!is.null(ntrain)){
		bb<-quantile(ntrain,probs=seq(0,1,length.out=nq+1))
		bb[1]<-0
		bb[nq+1]<-max(range(c(ntrain)))
		class<-cut(ntrain,breaks=bb,include.lowest=TRUE,labels=FALSE)

		qq<-round(quantile(npeaks,seq(0,1,length.out=nq+1)))
		d3<-data.frame(ntrain=npeaks)
		colnames(d3)<-"Number of training sequences"
		uu<-data.frame(x=10^((log10(qq[1:nq])+log10(qq[2:(nq+1)]))/2),y=0.1,label=1:nq)
		pA<-ggplot(data=d3,aes(x=`Number of training sequences`)) + geom_density() + geom_vline(xintercept=qq) + geom_text(data=uu,aes(x=x,y=y,label=label)) + scale_x_continuous(trans="log10",breaks=qq,labels=qq) + theme_bw(base_size=8) + ylab("Density") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		
		aa<-t(apply(data,1,rank))
		MRANK<-as.data.frame(t(sapply(1:nq,function(x){apply(aa[which(class==x),],2,avgfun)})))
		MRANK$class<-1:nq
		d4<-melt(MRANK,id.vars="class")
		names(d4)<-c("Number of training sequences quantiles","Method","rank")
		pB<-ggplot(data=d4,aes(x=`Number of training sequences quantiles`,y=`rank`,color=Method)) + geom_line() + geom_point(size=3,alpha=0.5) + scale_color_manual(values=apal) + scale_x_continuous(breaks=seq(1,nq,1),labels=1:nq) + scale_y_continuous(breaks=1:length(toconsider),labels=length(toconsider):1,limits=c(1,length(toconsider))) + ylab(if(avg=="mean"){"Average rank"}else{"Median rank"}) + guides(color = guide_legend(position="bottom", nrow = 2)) + theme_bw(base_size=8) + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
		legend<-as_ggplot(get_legend(pB))
		pB<-ggplot(data=d4,aes(x=`Number of training sequences quantiles`,y=`rank`,color=Method)) + geom_line() + geom_point(size=3,alpha=0.5) + scale_color_manual(values=apal) + scale_x_continuous(breaks=seq(1,nq,1),labels=1:nq) + scale_y_continuous(breaks=1:length(toconsider),labels=length(toconsider):1,limits=c(1,length(toconsider))) + ylab(if(avg=="mean"){"Average rank"}else{"Median rank"}) + guides(color = "none") + theme_bw(base_size=8) + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

		best_method<-ncol(data)
		best_method_name<-names(data)[best_method]
		data2<-data[,-best_method]-data[,best_method]
		#MAUC<-as.data.frame(t(sapply(1:nq,function(x){apply(data2[which(class==x),],2,function(y){median(y-data[which(class==x),best_method])})})))
		MAUC<-as.data.frame(t(sapply(1:nq,function(x){apply(data2[which(class==x),],2,median)})))
		MAUC$class<-1:nq
		d5<-melt(MAUC,id.vars="class")	
		names(d5)<-c("Number of training sequences quantiles","Method","Average_diff")
		pC<-ggplot(data=d5,aes(x=`Number of training sequences quantiles`,y=Average_diff,color=Method)) + geom_line() + geom_point(size=3,alpha=0.5) + scale_color_manual(values=apal) + scale_x_continuous(breaks=seq(1,nq,1),labels=1:nq) + ylab(paste("Median ",yname," - ", yname,"(",names(data)[best_method],")",sep="")) + guides(color = "none") + theme_bw(base_size=8) + theme(panel.grid.minor.y = element_blank(), panel.grid.minor = element_blank())
		res<-list(p1=p1,p2=p2,pA=pA,pB=pB,pC=pC,legend=legend,data=data,rup=rup)
		res<-list(main=ggarrange(p1,ggarrange(ggarrange(pA,pB,pC,nrow=3,heights=c(1,2,2),align="v"),legend,nrow=2,heights=c(5,1)),ncol=2),sig=p2,data=data,rup=rup)
	}else{
		res<-list(main=p1,sig=p2,data=data,rup=rup)
	
	}
	res
}

### BEGIN TF_BENCHMARK
load(paste(basedir,"/DATASETS/TF_ENCODE/einfo.Rdata",sep=""))
npeaks<-einfo$npeaks
load(paste(basedir,"/DATASETS/TF_BENCHMARK/einfo.Rdata",sep=""))
ntrain<-einfo$nTrainPos


setwd(paste(basedir,"/TFBS_PREDICTION/TF_BENCHMARK",sep=""))


#load("meme.Rdata")
load("seqgl.Rdata")
SeqGL<-seqgl
load("gkm.Rdata")
load("kdm.Rdata")

CRPTS<-read.delim("other_methods_results/CRPTS.csv",sep=",")
DanQ<-read.delim("other_methods_results/DanQ.csv",sep=",")
DeepBind<-read.delim("other_methods_results/DeepBind.csv",sep=",")
DLBSS<-read.delim("other_methods_results/DLBSS.csv",sep=",")
D_SSCA<-read.delim("other_methods_results/D_SSCA.csv",sep=",")
DeepSTF<-read.delim("other_methods_results/DeepSTF.csv",sep=",")
MLSNet<-read.delim("other_methods_results/MLSNet.csv",sep=",")
BERT<-read.delim("other_methods_results/BERT.tsv")

acc<-data.frame(
		CRPTS=CRPTS$ACC,
		DanQ=DanQ$ACC,
		DeepBind=DeepBind$ACC,
		DLBSS=DLBSS$ACC,
		D_SSCA=D_SSCA$ACC,
		DeepSTF=DeepSTF$ACC,
		MLSNet=MLSNet$ACC,
		KDM_LRLM=kdm$ACC,
		SeqGL=SeqGL$ACC,
		lsGKM=gkm$ACC,
		BERT_TFBS=BERT$Accuracy)

auroc<-data.frame(
		CRPTS=CRPTS$ROC.AUC,
		DanQ=DanQ$ROC.AUC,
		DeepBind=DeepBind$ROC.AUC,
		DLBSS=DLBSS$ROC.AUC,
		D_SSCA=D_SSCA$ROC.AUC,
		DeepSTF=DeepSTF$ROC.AUC,
		MLSNet=MLSNet$ROC.AUC,
		KDM_LRLM=kdm$AUROC,
		SeqGL=SeqGL$AUROC,
		lsGKM=gkm$AUROC,
		BERT_TFBS=BERT$ROC.AUC)

auprc<-data.frame(
		CRPTS=CRPTS$PR.AUC,
		DanQ=DanQ$PR.AUC,
		DeepBind=DeepBind$PR.AUC,
		DLBSS=DLBSS$PR.AUC,
		D_SSCA=D_SSCA$PR.AUC,
		DeepSTF=DeepSTF$PR.AUC,
		MLSNet=MLSNet$PR.AUC,
 		KDM_LRLM=kdm$AUPRC,
		SeqGL=SeqGL$AUPRC,
		lsGKM=gkm$AUPRC,
		BERT_TFBS=BERT$PR.AUC)

colnames(auroc)<-c("CRPTS","DanQ","DeepBind","DLBSS","D-SSCA","DeepSTF","MLSNet","KDM-LRLM","SeqGL","ls-GKM","BERT-TFBS")
colnames(auprc)<-c("CRPTS","DanQ","DeepBind","DLBSS","D-SSCA","DeepSTF","MLSNet","KDM-LRLM","SeqGL","ls-GKM","BERT-TFBS")
colnames(acc)<-c("CRPTS","DanQ","DeepBind","DLBSS","D-SSCA","DeepSTF","MLSNet","KDM-LRLM","SeqGL","ls-GKM","BERT-TFBS")

nquantiles<-4
toc<-c("CRPTS","DanQ","DeepBind","DLBSS","D-SSCA","DeepSTF","MLSNet","KDM-LRLM","SeqGL","ls-GKM","BERT-TFBS")
avgfun<-"mean"

tests<-list()
dd<-auroc
res1<-makePlot(dd,ntrain,nq=nquantiles,"auROC",toconsider=toc,pval=0.05,npeaks=npeaks*0.8,alternative="two",avg=avgfun)
pAUROC<-res1$main
tests[[1]]<-res1$sig

dd<-auprc
res2<-makePlot(dd,ntrain,nq=nquantiles,"auPRC",toconsider=toc,pval=0.05,npeaks=npeaks*0.8,alternative="two",avg=avgfun)
pAUPRC<-res2$main
tests[[2]]<-res2$sig

dd<-acc
res3<-makePlot(dd,ntrain,nq=nquantiles,"Accuracy",toconsider=toc,pval=0.05,npeaks=npeaks*0.8,alternative="two",avg=avgfun)
pACC<-res3$main
tests[[3]]<-res3$sig

ggarrange(plotlist=tests,ncol=3)
### END TF_BENCHMARK

### BEGIN TF_ENCODE
setwd(paste(basedir,"/TFBS_PREDICTION/TF_ENCODE",sep=""))
load("gkm.Rdata")
load("seqgl.Rdata")
# load("kdm.Rdata")
# kdm_win_win<-kdm
load("kdm_no_W_win.Rdata")
kdm_win<-kdm


auroc2<-data.frame(SeqGL=seqgl$AUROC,ls_GKM=gkm$AUROC,KDM_LRLM=kdm_win$AUROC)
auprc2<-data.frame(SeqGL=seqgl$AUPRC,ls_GKM=gkm$AUPRC,KDM_LRLM=kdm_win$AUPRC)
acc2<-data.frame(SeqGL=seqgl$ACC,ls_GKM=gkm$ACC,KDM_LRLM=kdm_win$ACC)

names(auroc2)<-names(auprc2)<-names(acc2)<-c("SeqGL","ls-GKM","KDM_LRLM")
toc<-c("SeqGL","ls-GKM","KDM_LRLM")
avgfun<-"mean"
res1<-makePlot(auroc2,ntrain=NULL,nq=4,"auROC",toconsider=toc,pval=0.05,alternative="two",avg="median")
res2<-makePlot(auprc2,ntrain=NULL,nq=4,"auROC",toconsider=toc,pval=0.05,alternative="two",avg="median")
res3<-makePlot(acc2,ntrain=NULL,nq=4,"auROC",toconsider=toc,pval=0.05,alternative="two",avg="median")
tests<-list()
tests[[1]]<-res1$sig
tests[[2]]<-res2$sig
tests[[3]]<-res3$sig
ggarrange(plotlist=tests,ncol=3)
### END TF_ENCODE
