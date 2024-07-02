
if( Sys.info()['nodename']%in%c('aulonocara.bp.lnf.it')){

  DATA_DIR="/adat/Progetti/RBP/profile"
  #CGB_DIR=''
  #GIT_LOCAL_DIR=''
  #CGB_SHARED='~/Dropbox (HuGeF)/SPLICING/CGB_shared/'

  }

setwd(DATA_DIR)



library(rkdMotifs)
library(parallel)
library(dplyr)
library(stringr)
source("/adat/Progetti/KDM//FACTORBOOK2/RTX/funzioni_importanti.rtx")
library(optparse)

option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="Project name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

protein<-opt$project
message(paste0("Project: ",protein))

flanks<-61
hwin<-15

### Obtain motifs info ###

W<-readRDS(paste0("Model/W_",protein,"_6_2_70.rds"))

pos<-kdmReadBedFile(filename=paste0("Sequences/",protein,"_150.bed"),type="bed6")
centers<-floor((pos$end+pos$start)/2)-pos$start
temp<-pos
temp$start<-pos$start+centers-flanks
temp$end<-pos$start+centers+flanks+1
pseq<-kdmGetSequence(temp,"/adat/database/2bit/hg19.2bit")

neg<-kdmReadBedFile(filename=paste0("Sequences/",protein,"_negatives_150_strand.bed"),type="bed6")
centers<-floor((neg$end+neg$start)/2)-neg$start
temp<-neg
temp$start<-neg$start+centers-flanks
temp$end<-neg$start+centers+flanks+1
nseq<-kdmGetSequence(temp,"/adat/database/2bit/hg19.2bit")


pinfo<-kdmSeqProfileInfo(pseq,nseq,W,half.interval=flanks,half.win=hwin)
saveRDS(pinfo,paste0("Results/",protein,"_pinfo.rds"))
message(paste0("Results/",protein,"_pinfo.rds"))
enr<-kdmEnrichment(pinfo)


### Statisitcal analysis of the motifs: centrality and enrichment ###

RES<-list()
for(m in 1:70){
	E<-enr[[m]]
	E$start<-E$start-46
	ss<-cbind(which(E$start==-45 & E$win==30),which(E$start==-15 & E$win==30),which(E$start==+15 & E$win==+30))
	RES[[m]]<-E[cbind(which(E$start==-45 & E$win==30),which(E$start==-15 & E$win==30),which(E$start==+15 & E$win==+30)),c("enr","enr_q.value","cen","cen_q.value")]
}
saveRDS(RES,paste0("Results/",protein,"_test.rds"))
message(paste0("Results/",protein,"_test.rds"))

### Convert kdm motifs in normal PWMs ###


kdmTopwm<-function(motifs,nmat=1,nflank=5,base=2,alpha=1e-6,cutoff=2,bkg=c(0.25,0.25,0.25,0.25),pseudo=1,um=TRUE,fraction=0.9,estimate=TRUE,threshold=0){
	map<-list()
	map[["a"]]<-c(1,0,0,0)
	map[["c"]]<-c(0,1,0,0)
	map[["g"]]<-c(0,0,1,0)
	map[["t"]]<-c(0,0,0,1)
	map[["n"]]<-bkg	
	map[["A"]]<-c(1,0,0,0)
	map[["C"]]<-c(0,1,0,0)
	map[["G"]]<-c(0,0,1,0)
	map[["T"]]<-c(0,0,0,1)
	map[["N"]]<-bkg		

	getPWM<-function(seq){
		matrix(unlist(map[strsplit(seq,split="")[[1]]]),nrow=4)
	}

    update_ss<-function(kmers,PWM,cutoff,pseudo,bkg,base,klength,alpha){
		###ALIGN
		res<-data.frame(seq=kmers[,1],weight=kmers$weight,pos=NA,score=NA)
		PWMM<-matrix(0,4,ncol(PWM))
		scol<-which(apply(PWM,2,sum)>0)
		
		
		PWMM<-apply(PWM*1000+pseudo,2,function(a){log((a/sum(a)/bkg),base)})
		maxpos<-ncol(PWMM)-klength+1
		for(i in 1:nrow(kmers)){
			pwm<-getPWM(kmers[i,1])
			scores<-rep(0,maxpos)
			for(pos in 1:maxpos){
				scores[pos]<-sum(apply(pwm*PWMM[,pos:(pos+klength-1)],2,sum))
			}
			
			res$pos[i]<-which.max(scores)
			res$score[i]<-scores[res$pos[i]]
		}
		
		###UPDATE
		sel<-which(res$score>cutoff)
		newPWM<-matrix(0,4,ncol(PWMM))
		for(i in 1:length(sel)){
			k<-sel[i]
			start<-res$pos[k]
            newPWM[,start:(start+klength-1)]<-newPWM[,start:(start+klength-1)] + getPWM(res$seq[k])*res$weight[k] + alpha
		}
		newPWM<-newPWM[,which(apply(newPWM,2,sum)>0)]
		newPWM<-apply(newPWM,2,function(a){a/sum(a)})
		rownames(newPWM)<-c("A","C","G","T")
		list(PWM=newPWM,selected=sel,res=res)
	}
    
	cbkg<-(bkg*1000+pseudo)/sum(bkg*1000+pseudo)
	PWMS<-list();
	if(motifs@k>0 & estimate){
		motifs<-kdmEstimateMotifs(motifs)
	}
	
	if(estimate){
		kk1<-kdmGkmersString(kmerlength=motifs@l,ngaps=0,double.strand=FALSE)
	}else{
		kk1<-kdmGkmersString(kmerlength=motifs@l,ngaps=motifs@k,double.strand=FALSE)
	}
	if(motifs@doubleStrand==FALSE){
		#if(estimate){
			#kk2<-kdmGkmersString(kmerlength=motifs@l,ngaps=0,double.strand=TRUE)
		#}else{
			#kk2<-kdmGkmersString(kmerlength=motifs@l,ngaps=motifs@k,double.strand=TRUE)
		#}
		for(k in 1:ncol(motifs)){
			PWMS[[k]]<-list()
			kk1$score<-as.numeric(motifs[,k]^2)
# 			kk2$weight<-alpha*as.numeric(motifs[,k]^2)
			kk1$weight<-as.numeric(motifs[,k]^2)
			
			oo<-order(kk1$score,decreasing=TRUE)
			kmers<-kk1[oo,]
			#kmers<-kmers[1:floor(nrow(kk2)*fraction),]
# 			print(cumsum(kk2$score[oo]))
			last<-min(which(cumsum(kk1$score[oo])>fraction))
			cat("PPPPPP---------------------->",last,"\n")
			kmers<-kmers[1:last,]
# 			browser()
			
			
			for(m in 1:nmat){
				#PWM0<-cbind(matrix(0,4,nflank),getPWM(kmers$fwd[1]),matrix(0,4,nflank))
				PWM0<-cbind(matrix(0,4,nflank),getPWM(kmers[,1]),matrix(0,4,nflank))
				r<-1000000
# 				while(r>0.0001){
				while(r>threshold){
					PWM<-update_ss(kmers,PWM0,cutoff=cutoff,pseudo=pseudo,bkg=cbkg,base=base,motifs@l,alpha=alpha)
					if(ncol(PWM$PWM)==ncol(PWM0)){
						r<-sqrt(sum( (PWM0-PWM$PWM)^2))
					}
					PWM0<-PWM$PWM
					print(PWM0)
					cat("matrix=",m,"\tchange=",r,"\tnkmers=",length(PWM$selected),"\n",sep="")
				}
				PWMS[[k]][[m]]<-PWM$PWM
#  				print(PWM$PWM)
# 				print(PWM$res)
				
				kmers<-kmers[-PWM$selected,]
			}
			if(um){
				#aa<-list()
				for(i in 1:length(PWMS[[k]])){
					aa<-create_motif(PWMS[[k]][[i]],alphabet="DNA",type="PPM",name=paste("GKM.",k,sep=""),bkg=bkg)
					aa@alphabet<-"DNA"
				}
				PWMS[[k]]<-aa
			}
		}
	}else{
		stop("Not yet implemented for double strand motifs")
	}
	PWMS
}


pwms<-kdmTopwm(W,nflank=3,cutoff=1,um=FALSE,fraction=0.99,estimate=FALSE)
saveRDS(pwms,paste0("Results/",protein,"_pwms.rds"))
message(paste0("Results/",protein,"_pwms.rds"))

### Calculate coefficients roles ###

saveRDS(read.table(paste0("coeff/",protein,"_6_2_70.txt"),header=TRUE)%>%filter (model!="intercept") %>% mutate(model=factor(model,levels=seq(1,70)))%>%group_by(model,window) %>% summarise(mean_coeff=mean(coeff),sd=sd(coeff,na.rm=TRUE)),paste0("Results/",protein,"_coeff.rds"))
message(paste0("Results/",protein,"_coeff.rds"))


### Compare kdm motifs with mCross reference ###

mcross<-readRDS("reference_k562_ss.rds")
k562<-kdmLoadPWMSet("/adat/Progetti/RBP/mCrossBase_K562_pwms.meme")
colnames(mcross)<-c(names(k562),"5ss","3ss")
kdmtom<-kdmTom(W,mcross,threshold=0.5,best.hit=FALSE)
kdmtom$mcross<-do.call(rbind,str_split(kdmtom$target,"\\."))[,2]
kdmtom<-kdmtom[,c("query","rank","mcross")]%>%filter(rank<=10)%>%group_by(query)%>%summarise(mCross=paste0(unique(mcross),collapse="-"))
saveRDS(kdmtom,paste0("Results/",protein,"_kdmtom.rds"))
message(paste0("Results/",protein,"_kdmtom.rds"))
