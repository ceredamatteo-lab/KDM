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

# P15.PLOT ----

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
result=readRDS(paste0(data_folder,"p15.rds"))

p1<-plotClassification(result[,-1],"KDMtom",5)

pdf(paste0(plot_folder,"P15_RBP.pdf"),height = 8,width = 10)
ggarrange(plotlist=p1,nrow=1)
dev.off()

classification=readRDS("tomtom_kdmtom/Classification_comparison.rds")
colors=c("grey90","grey80",brewer.pal(5,"Reds"))
names(colors)=rev(c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference"))


class_1_summary=melt(classification,id.vars="target")%>%
  mutate(value=factor(value,levels=rev(c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference")))) %>%count(variable, value)%>%
  group_by(variable)%>%mutate(tot=sum(n),prop = n / tot)%>%
  arrange(variable, desc(value))%>%mutate(cum=cumsum(n),cprop=cumsum(prop),label = paste0(n,"(",round(prop * 100), "%)"),ypos = cumsum(prop) - prop / 2)

class_1_summary$test=""
r="Tomtom"
m="Kdmtom"

for(i in 1:nrow(class_1_summary)){
	if(class_1_summary$variable[i]==r){next}
	class=class_1_summary$value[i]
	j=which(class_1_summary$value==class & class_1_summary$variable==r)
	if(length(j)==0){j=old}
	test=binom.test(x=class_1_summary$cum[i],n=class_1_summary$tot[i],
	p=class_1_summary$cum[j]/class_1_summary$tot[j])
	p.val=test$p.value
	if(p.val>0.05){class_1_summary$test[i]="="
	} else if (test$estimate>class_1_summary$cum[j]/class_1_summary$tot[j]){class_1_summary$test[i]="+"
	} else {class_1_summary$test[i]="-"}
	old=j
}

p1=ggplot(class_1_summary, aes(x = variable, y = prop, fill = value)) +
    geom_bar(stat = "identity", col = "black", position = "stack") +
    geom_text(aes(label = paste0(label," ",test), y = ypos), size = 3, color = "black",fontface="bold") +
    theme_minimal() +
    theme(axis.title = element_blank(),panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("CISBP+ATTRACT against mCross")



tmp=classification
n1=colnames(tmp)[2]
n2=colnames(tmp)[3]
tmp=data.frame(table(tmp[,2],tmp[,3]))%>%
  mutate(Var1=factor(Var1,levels=c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference")),
         Var2=factor(Var2,levels=c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference")))
tmp=merge(tmp,tmp%>%group_by(Var1)%>%summarise(tot.Var1=sum(Freq)),by="Var1") 
tmp=merge(tmp,tmp%>%group_by(Var2)%>%summarise(tot.Var2=sum(Freq)),by="Var2") 
tmp=tmp%>%mutate(Jaccard=Freq/(tot.Var1+tot.Var2-Freq))

p2=ggplot(tmp,aes(x=Var1,y=Var2,fill=Jaccard))+
  #geom_tile(col="black")+
  scale_y_discrete(limits=rev)+
  geom_point(col="black",shape=21,size=9)+
  coord_equal()+
  geom_text(aes(label=Freq))+scale_fill_gradient(low = "white",high = "forestgreen")+theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1))+xlab(n1)+ylab(n2)
  
  
ggarrange(p1,p2,ncol=1)
	
