library(rkdMotifs)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(universalmotif)

# P15.DATA ----
setwd("/adat/Progetti/KDM/MATERIALE")
source("/adat/Progetti/KDM/FACTORBOOK2/RTX/funzioni_importanti.rtx")
load("RBP/ENCODE_eCLIP_DATASET.Rdata")
load("RBP/mCrossBase.Rdata")

queryfile="RBP/mCrossBase.meme"
target_info=readRDS("RBP/CISBP_ATTRACT_info.rds")
targetfile="RBP/CISBP_ATTRACT.meme"

## TOMTOM ----
bdir<-"/tmp/tempresRBP/"
if (!dir.exists(bdir)) {dir.create(bdir, recursive = TRUE)}

cmd_tomtom<-paste("tomtom -oc ",bdir,"tomtom/ -xalph -thresh 0.05 -norc ",targetfile," ",queryfile,sep="")
system(cmd_tomtom,intern=TRUE,ignore.stderr=TRUE)
#rr<-system(cmd_tomtom,intern=TRUE,ignore.stderr=TRUE)
	 		
tom<-read.delim(paste0(bdir,"/tomtom/tomtom.tsv"),comment.char="#")	
tom$target2<-mcross_info$target[match(tom$Target_ID,mcross_info$name)]
tom$query2<-target_info$RBP_Name[match(tom$Query_ID,target_info$Motif_ID)]
tom$rank<-NA

classification<-data.frame(target=unique(target_info$RBP_Name),Tomtom="No match")
for(u in 1:nrow(classification)){
	rbp=classification$target[u]
	if(!rbp %in% tom$query2){next}
	
	sel=subset(tom,query2==rbp)
	un<-unique(sel$Query_ID)
	for(i in 1:length(un)){
		ss<-which(sel$Query_ID==un[i])
		sel$rank[ss]<-1:length(ss)
		}

	#sel<-sel[which(sel$q.value * length(un) < 0.05),]
	ss<-subset(sel,target2==query2)
	if(nrow(ss)==0){next}
	classification$Tomtom[u]<-as.character(cut(min(ss$rank),c(0,1,3,10,Inf)))
}
	
## KDMTOM ----
pwm<-kdmLoadPWMSet(targetfile)
kdm_pwm<-kdmMotifsFromPWM(pwm,6,2,FALSE,5)						
colnames(kdm_pwm)<-gsub(" ","",names(pwm))

tom<-kdmTom(kdm_pwm,mcross_W,best.hit=FALSE,threshold=0.05,targetInfo=mcross_tinfo)
tom$target2<-mcross_info$target[match(tom$target,mcross_info$name)]
tom$query2<-target_info$RBP_Name[match(tom$query,target_info$Motif_ID)]
classification$Kdmtom="No match"
for(u in 1:nrow(classification)){
	rbp=classification$target[u]
	if(!rbp %in% tom$query2){next}
	
	sel=subset(tom,query2==rbp)
	un<-unique(sel$query)
	for(i in 1:length(un)){
		ss<-which(sel$query==un[i])
		sel$rank[ss]<-1:length(ss)
		}

	#sel<-sel[which(sel$q.value * length(un) < 0.05),]
	ss<-subset(sel,target2==query2)
	if(nrow(ss)==0){next}
	classification$Kdmtom[u]<-as.character(cut(min(ss$rank),c(0,1,3,10,Inf)))
}
classification[classification=="No match"]="No match to target"
colnames(classification)[2:3]=c("CISBP-ATTRACT_mCross_Tomtom","CISBP-ATTRACT_mCross_KDMtom")
saveRDS(classification,"RBP/Paper_Figure/p15.rds")

## P15.PLOT ----

if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local" ){
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
result=readRDS(paste0(data_folder,"p15.rds"))

p1<-plotClassification(result[,-1],"KDMtom",5)

pdf(paste0(plot_folder,"P15_RBP.pdf"),height = 8,width = 10)
ggarrange(plotlist=p1,nrow=1)
dev.off()

