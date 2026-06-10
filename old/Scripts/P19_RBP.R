library(rkdMotifs)
library(universalmotif)
## P19.DATA ----
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")

hw=6
half_interval=70
original_length=35
genomeFile<-"/adat/database/2bit/hg38.2bit"

kdmMotifsReduce<-function(motifs,profile,motifSelection=NULL,use_negatives=TRUE,weights=NULL,maxdist=5,siglev=1e-4,bonf=TRUE){
	require(cluster)
    if(is.null(motifSelection)){
		motifSelection=1:ncol(motifs)
	}
    if(is.null(weights)){
		weights=rep(1,length(motifSelection))
	}
	
	if(length(weights)!=length(motifSelection)){
		error("motifSelection and weights must have the same lenght")
	}
	
	nm<-length(motifSelection)
	npos<-length(profile$x)
	oclusters<-list()
	oweights<-list()
	GG<-EE<-PP<-NULL
	if(nm>1){
		if(use_negatives){
			occurrences<-cbind(profile$phits[motifSelection,],profile$nhits[motifSelection,])
		}else{
			occurrences<-profile$phits[motifSelection,]
		}
		nseqs<-ncol(occurrences)
		
		csig<-siglev
		if(bonf){
			csig<-siglev/((nm^2-nm)/2)
		}
		GG<-EE<-PP<-matrix(0,nm,nm)
		diag(GG)<-1
		diag(EE)<-1
		diag(PP)<-1
		intv<-2:(npos-1)
		ppeak<-apply(occurrences,1,function(a){sum(a>1 & a<npos,na.rm=TRUE)})/nseqs/(npos-2)
		mult<-2*sum(maxdist+1:maxdist) + (2*maxdist+1) * (npos-2-2*maxdist)
		for(ma in 1:(nm-1)){
			for(mb in (ma+1):nm){
# 				if(sign(weights[ma])==sign(weights[mb])){
					prob<-mult*ppeak[ma]*ppeak[mb]
					selp<-which(occurrences[ma,]>1 & occurrences[ma,]<npos & occurrences[mb,]>1 & occurrences[mb,]<npos)
					if(length(selp)>0){
						CCP<-sum(abs(occurrences[ma,selp]-occurrences[mb,selp])<maxdist)
						pv<-binom.test(x=CCP,n=length(selp),p=prob,alternative="gr")$p.value
						GG[ma,mb]<-GG[mb,ma]<-as.numeric(pv<csig)
						EE[ma,mb]<-EE[mb,ma]<-(CCP/length(selp))/prob
						PP[ma,mb]<-PP[mb,ma]<-pv
					}
# 				}
			}
		}

		colnames(GG)<-rownames(GG)<-motifSelection
		colnames(EE)<-rownames(EE)<-motifSelection
		colnames(PP)<-rownames(PP)<-motifSelection
		a<-mean(GG[upper.tri(GG)])
		if(a==0){
			WR<-W[,motifSelection]
			for(i in 1:length(motifSelection)){
				oclusters[[i]]<-motifSelection[i]
				oweights[[i]]<-weights[i]
			}
		}else if(a==1){ 
			oclusters[[1]]<-motifSelection
			oweights[[1]]<-weights
			tw<-matrix(weights,ncol=1)
			WR<-new("kdmotifs",sqrt((W[,oclusters[[1]]]^2 %*% tw)/sum(tw)),l=W@l,k=W@k,doubleStrand=W@doubleStrand,estimated=W@estimated)
		}else{
			D<-dist(GG,method="binary");
			hc<-hclust(D,method="complete")
			sil<-numeric(length(motifSelection))
			for(i in 2:(length(motifSelection)-1)){sil[i]<-mean(silhouette(cutree(hc,k=i),dist=D)[,3])}
			nc<-which.max(sil)
			clusters<-cutree(hc,k=nc)
			clus<-unique(clusters)
			wtmp<-matrix(0,nrow(W),length(clus))
			names<-c();
			for(i in 1:length(clus)){
				sel<-which(clusters==clus[i])
				oclusters[[i]]<-motifSelection[sel]
				oweights[[i]]<-weights[sel]
				if(length(sel)>1){
					tw<-matrix(oweights[[i]],ncol=1)
					wtmp[,i]<-sqrt((W[,oclusters[[i]]]^2 %*% tw)/sum(tw))
					names<-c(names,oclusters[[i]])
				}else{
					wtmp[,i]<-W[,oclusters[[i]]]
					names<-c(names,oclusters[[i]])
				}
				
			}
			WR<-new("kdmotifs",wtmp,l=W@l,k=W@k,doubleStrand=W@doubleStrand,estimated=W@estimated)
		}
	}else{
		GG<-EE<-PP<-NULL
		WR<-motifs[,motifSelection]
		oclusters[[1]]<-motifSelection
		oweights[[1]]<-weights
	}
	list(WR=WR,clusters=oclusters,weights=oweights,GG=GG,EE=EE,PP=PP)
}


### BEGIN CLASSIFICATION MCROSSBASE SELECTION


### BEGIN STREME POS_NEG_ALL CLASSIFICATION VS MCROSSBASE SELECTION TRAIN
mcrossTF<-unique(mcross_info$target)
basedir<-"/adat/Progetti/RBP/MF/STREME/POS_NEG_ALL/"
classification<-data.frame(target=exp_info$target,cell=exp_info$cell,class="No match",targetcell=NA,nmotifs=0)
for(u in 1:nrow(exp_info)){genomeFile<-"/adat/database/2bit/hg38.2bit"
    print(u)
	if(exp_info$target[u] %in% mcrossTF){
		load(paste(basedir,exp_info$Experiment[u],"_train.Rdata",sep=""))
		if(res$nmotifs>0){
		npos=exp_info$nTrainPeaks[u]
		nneg=npos
		centrimo=res$centrimo
		all_selected_motifs=c()
		if(nrow(centrimo)>0){
            centrimo$enrichment<-(centrimo$sites_in_bin/npos)/(centrimo$neg_sites_in_bin/nneg)
            for(i in 1:nrow(centrimo)){
            centrimo$enrichment_p.value[i]<-binom.test(x=centrimo$sites_in_bin[i],n=npos,p=centrimo$neg_sites_in_bin[i]/nneg,alternative="gr")$p.value
            centrimo$enrichment_q.value[i]<-centrimo$enrichment_p.value[i] * nrow(centrimo)
            }
            centrimo<-centrimo[order(rank(-log10(centrimo$E.value))+rank(-log10(centrimo$enrichment_q.value)),-log10(centrimo$E.value),decreasing=TRUE),]
            res$centrimo=centrimo	
            all_selected_motifs<-res$centrimo$motif_id[which(res$centrimo$E.value<0.05 & res$centrimo$enrichment_q.value <0.05 & res$centrimo$enrichment >1 & res$centrimo$bin_location<=0)]
        }
		if(length(all_selected_motifs)>0){
			selected_motifs<-all_selected_motifs[1:min(5,length(all_selected_motifs))]
			classification$nmotifs[u]<-length(selected_motifs)
			tom<-res$tomtom
			if(nrow(tom)>0){
				tom$target2<-mcross_info$target[match(tom$Target_ID,mcross_info$name)]
				tom$cell<-mcross_info$cell[match(tom$Target_ID,mcross_info$name)]
				tom$rank<-NA
				un<-unique(tom$Query_ID)
				for(i in 1:length(un)){
					ss<-which(tom$Query_ID==un[i])
					tom$rank[ss]<-1:length(ss)
				}
			}
			mm<-which(tom$Query_ID %in% selected_motifs)
			if(length(mm)>0){
				tom<-tom[mm,]
 				tom<-tom[which(tom$q.value * length(selected_motifs) < 0.05),]
				if(nrow(tom)>0){
					ss<-which(tom$target2==exp_info$target[u])
					if(length(ss)>0){
						classification$class[u]<-as.character(cut(min(tom$rank[ss]),c(0,1,3,10,Inf)))
						classification$targetcell[u]<-tom$cell[ss[which.min(tom$rank[ss])]]
					}
				}
			}
		}else{
			classification$class[u]<-"No central motif"
		}
	}else{
		classification$class[u]<-"No motif identified"
		}
	}else{
		classification$class[u]<-"Not in reference"
	}
}
class_streme_train<-classification
### END   STREME POS_NEG_ALL CLASSIFICATION VS MCROSSBASE SELECTION TRAIN

### BEGIN KDM_REDUCED POS_NEG_ALL 70 CLASSIFICATION VS MCROSSBASE SELECTION TRAIN
mcrossTF<-unique(mcross_info$target)
basedir<-"/adat/Progetti/RBP/MF/W/W_POS_NEG_ALL/70/"
classification<-data.frame(target=exp_info$target,cell=exp_info$cell,class="No match",targetcell=NA,nmotifs=0)
for(u in 1:nrow(exp_info)){
    print(u)
	if(exp_info$target[u] %in% mcrossTF){
		modelfile<-paste(basedir,exp_info$Experiment[u],"_train.rds",sep="")
		if(file.exists(modelfile)){
			W<-as(readRDS(modelfile),"kdmotifs")
			colnames(W)<-1:ncol(W)
			pn<-kdmGetProfileInfo(motifs=W,regions=DATASET[[u]]$train,centers=rep(original_length+1,nrow(DATASET[[u]]$train)),labels=DATASET[[u]]$trainLabels,genomeFile=genomeFile,
			halfInterval=half_interval,halfWin=hw,tolerance=1e-10,strict=FALSE,use_float=FALSE)
			ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
			res<-c()
			for(k in 1:ncol(W)){
				E<-ne[[k]]
				cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
				cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
				#cond<-which((E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
				if(length(cond)>0){
					#ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
					ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),-log10(E$Centrality_q.value[cond]),decreasing=TRUE)],]
					res<-rbind(res,data.frame(motif=k,ss[1,]))
					
				}
			}
			if(length(res)>1){
				#res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
            res<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value)),-log10(res$Centrality_q.value),decreasing=TRUE),]
            MR<-kdmMotifsReduce(W,pn,motifSelection=res$motif,use_negatives=TRUE,weights=NULL,maxdist=5,siglev=1e-4,bonf=TRUE)
				W<-MR$WR
				colnames(W)<-1:ncol(W)
				pn<-kdmGetProfileInfo(motifs=W,regions=DATASET[[u]]$train,centers=rep(original_length+1,nrow(DATASET[[u]]$train)),labels=DATASET[[u]]$trainLabels,genomeFile=genomeFile,
				halfInterval=half_interval,halfWin=hw,tolerance=1e-10,strict=FALSE,use_float=FALSE)
				ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
				res<-c()
				for(k in 1:ncol(W)){
					E<-ne[[k]]
					cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
					cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
					#cond<-which((E$Centrality_q.value * ncol(W)) < 0.05 & (E$Enrichment_q.value  * ncol(W))<0.05 & E$Centrality>1 & E$Enrichment>1)
					if(length(cond)>0){
						#ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
						ss<-E[cond[order(rank(-log10(E$Centrality_q.value[cond]))+rank(-log10(E$Enrichment_q.value[cond])),-log10(E$Centrality_q.value[cond]),decreasing=TRUE)],]
						res<-rbind(res,data.frame(motif=k,ss[1,]))
						
					}
				}
				if(length(res)>1){
					#res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
		            res<-res[order(rank(-log10(res$Centrality_q.value))+rank(-log10(res$Enrichment_q.value)),-log10(res$Centrality_q.value),decreasing=TRUE),]
					selected_motifs<-res$motif[1:min(5,nrow(res))]
					classification$nmotifs[u]<-length(selected_motifs)
				
					tom<-kdmTom(W,mcross_W,best.hit=FALSE,threshold=0.05,targetInfo=mcross_tinfo)
					if(length(tom)>1){
						tom$target2<-mcross_info$target[match(tom$target,mcross_info$name)]
						tom$cell<-mcross_info$cell[match(tom$target,mcross_info$name)]
						tom<-tom[which( (tom$q.value * length(selected_motifs) < 0.05) & (tom$query %in% selected_motifs) ),]
						ss<-which(tom$target2==exp_info$target[u])
						if(length(ss)>0){
							classification$class[u]<-as.character(cut(min(tom$rank[ss]),c(0,1,3,10,Inf)))
							classification$targetcell[u]<-tom$cell[ss[which.min(tom$rank[ss])]]
						}
					}
				}
			}else{
				classification$class[u]<-"No central motif"
			}
		}else{
			classification$class[u]<-"No motif identified"
		}
		#print(rbind(class_kdmr_train1[u,],classification[u,]))
	}else{
		classification$class[u]<-"Not in reference"
	}
}
class_kdmr_train<-classification
### END   KDM REDUCED POS_NEG_ALL 70 CLASSIFICATION VS MCROSSBASE SELECTION	TRAIN

tt2=cbind(class_kdmr_train$class,class_streme_train$class)
colnames(tt2)=c("KDMM_mCross_KDMtom","STREME_mCross_Tomtom")
tt2=data.frame(tt2)
tt2[tt2=="No match"]="No match to target"
tt2[tt2=="No central motif"]="No central motifs"
tt2[tt2=="Not in reference"]="Target not in reference"
saveRDS(tt2,"RBP/Paper_Figure/P19.rds")

### END CLASSIFICATION MCROSSBASE SELECTION

## P19.PLOT ----

if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
  script_folder="Scripts/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
  script_folder="PARTI_FIGURE/"
}


source(paste0(script_folder,"00_plot_functions.R"))
result=readRDS(paste0(data_folder,"P19.rds"))[,2:1]

p3<-plotClassification(result,"KDMtom",6)

pdf(paste0(plot_folder,"P19_RBP.pdf"),height = 8,width = 10)
ggarrange(plotlist=p3,nrow=1)
dev.off()









