library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)

## FUNCTIONS ## 
get_binom_test=function(data,condition,adjust=TRUE){
	data$p.val=NA
	data$sign=NA
	if(length(unique(data$variable))==2){r=unique(data$variable)[!grepl(condition,unique(data$variable))]}
	if(length(unique(data$variable))==3){r=unique(data$variable)[grep("Factorbook",unique(data$variable))]}
	m=unique(data$variable)[grepl(condition,unique(data$variable))]
	
	ref=subset(data,variable==r)%>%mutate(p.adj=NA)
	final=list(ref)
	for(var in unique(data$variable)){
		if(var==r){next}
		m=subset(data,variable==var)
		for(i in 1:nrow(m)){
			if(m$cum[i]==m$tot[i]){next}
			class=m$value[i]
			j=which(ref$value==class)
			if(length(j)==0){j=old}
			test=binom.test(x=m$cum[i],n=m$tot[i],p=ref$cum[j]/ref$tot[j])
			m$p.val[i]=test$p.value
			m$sign[i]=log10(test$estimate/(ref$cum[j]/ref$tot[j]))
			old=j
		}
		
		if(adjust){m$p.adj=p.adjust(m$p.val,method="fdr")
		} else {m$p.adj=m$p.val}
		
		final[[length(final)+1]]<-m
	}
	
	data=do.call(rbind,final)%>%mutate(test=ifelse(is.na(p.adj),"",ifelse(p.adj>0.05,"=",ifelse(sign>0,"+","-"))))
	return(data)
}

get_summary=function(data,motifs,reference,tool,classes_sel){
	mm=melt(data%>%mutate(ID=row_number()),id.vars="ID")%>%
	filter(value%in%classes_sel)%>%
  mutate(value=factor(value,levels=rev(classes_sel))) %>%count(variable, value)%>%
  group_by(variable)%>%mutate(tot=sum(n),prop = n / tot)%>%
  arrange(variable, desc(value))%>%mutate(cum=cumsum(n),cprop=cumsum(prop),label = paste0(n,"(",round(prop * 100), "%)"),ypos = cumsum(prop) - prop / 2)

	mm$Motifs<-NA
	mm$Reference<-NA
	mm$Tool<-NA
	for(i in 1:ncol(data)){
		ss<-which(mm$variable==colnames(data)[i])
		mm$Motifs[ss]<-paste(i,motifs[i],sep="")
		mm$Reference[ss]<-reference[i]
		mm$Tool[ss]<-tool[i]
	}
	
	return(mm)
}

get_jaccard=function(data,classes_sel){
	tmp=data.frame(table(data[,1],data[,2]))%>%
	  mutate(Var1=factor(Var1,levels=classes_sel),
		 Var2=factor(Var2,levels=classes_sel))
	tmp=merge(tmp,tmp%>%group_by(Var1)%>%summarise(tot.Var1=sum(Freq)),by="Var1") 
	tmp=merge(tmp,tmp%>%group_by(Var2)%>%summarise(tot.Var2=sum(Freq)),by="Var2") 
	tmp=tmp%>%mutate(Jaccard=Freq/(tot.Var1+tot.Var2-Freq))
	
	return(tmp)
}

plotClassification<-function(data,condition,sel,adjust=TRUE){
	classes=c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match to target","No central motifs","Target not in reference")
	classes_sel=classes[1:sel]
	colors=c("grey90","grey80",brewer.pal(5,"Reds"))
	names(colors)=rev(classes)
	data=data[which(data[,1]%in%classes_sel),]
	
	a<-strsplit(colnames(data),split="_")
	motifs<-unlist(lapply(a,function(x){x[1]}))
	reference<-unlist(lapply(a,function(x){x[2]}))
	tool<-unlist(lapply(a,function(x){x[3]}))
	xlabels<-motifs
	names(xlabels)<-paste(1:length(motifs),motifs,sep="")

	mm=get_summary(data,motifs,reference,tool,classes_sel)
	mm=get_binom_test(mm,condition,adjust)

	mm2<-data.frame(Motifs=paste(1:length(motifs),motifs,sep=""),Reference=as.factor(reference),Classification_method=as.factor(tool),y=-0.03)
	
	p1=ggplot(mm, aes(x = Motifs, y = prop, fill = value)) +
	    geom_bar(stat = "identity", col = "black", position = "stack") +
	    geom_text(aes(label = paste0(label," ",test), y = ypos), size = 3, color = "black",fontface="bold") +
	    geom_text(data=subset(mm,cum==tot),aes(label=tot,y=cprop),vjust=-1)+
	    theme_minimal() +
	    scale_fill_manual(values = colors)+
	    scale_x_discrete(labels=xlabels)+
	    geom_point(inherit.aes=FALSE,data=mm2,aes(x=Motifs,y=y,color=Reference,shape=Classification_method),size=5) +
	    xlab("Motif finding method") +
	    ylab("Experiments fraction") +
	    theme_minimal() +
	    theme(panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10))+
	    guides(color= guide_legend(override.aes = list(shape = 15)))
	    
	  if(ncol(data)==2){
	  jaccard=get_jaccard(data,classes_sel)
	  p2=ggplot(jaccard,aes(x=Var1,y=Var2,fill=Jaccard))+
	  scale_y_discrete(limits=rev)+
	  geom_point(col="black",shape=21,size=9)+
	  coord_equal()+
	  geom_text(aes(label=Freq))+scale_fill_gradient(low = "white",high = "forestgreen")+theme_bw()+
	  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1))+
	  xlab(colnames(data)[1])+ylab(colnames(data)[2])
	  
	  return(list(p1,p2))
	  }
	  
	  if(ncol(data)==3){
	  jaccard_plots=list()
	  for(i in 1:3){
	  	data2=data[,-i]
	  	jaccard=get_jaccard(data2,classes_sel)
	  	p2=ggplot(jaccard,aes(x=Var1,y=Var2,fill=Jaccard))+
		  scale_y_discrete(limits=rev)+
		  geom_point(col="black",shape=21,size=9)+
		  coord_equal()+
		  geom_text(aes(label=Freq))+scale_fill_gradient(low = "white",high = "forestgreen")+theme_bw()+
		  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1))+
		  xlab(colnames(data2)[1])+ylab(colnames(data2)[2])
		jaccard_plots[[i]]=p2
	  }
	  return(c(list(p1), jaccard_plots))}
	  
}

## TF ##

setwd("/adat/Progetti/KDM/MATERIALE/TF/CLASSIFICAZIONE")

load("class_all_fb_hocomoco_kdmtom.Rdata")
load("class_all_fb_hocomoco_tomtom.Rdata")

tt1<-data.frame(Factorbook_HOCOMOCOv13_Tomtom=class_all_fb_hocomoco_kdmtom[,3],Factorbook_HOCOMOCOv13_KDMtom=class_all_fb_hocomoco_tomtom[,3])


load("class_all_meme_hocomoco_tomtom.Rdata")
load("class_all_kdm_hocomoco_kdmtom.Rdata")

tt2<-data.frame(
	Meme_HOCOMOCOv13_Tomtom=class_all_meme_hocomoco_tomtom[,3],
	KDM_HOCOMOCOv13_KDMtom=class_all_kdm_hocomoco_kdmtom[,3]
)

#load("class_all_kdmr_hocomoco_kdmtom.Rdata")
load("class_all_kdmr_hocomoco_kdmtom_nosel.Rdata")

tt3<-data.frame(
	Meme_HOCOMOCOv13_Tomtom=class_all_meme_hocomoco_tomtom[,3],
	KDMM_HOCOMOCOv13_KDMtom=class_all_kdmr_hocomoco_kdmtom_nosel[,3]
)


load("class_sel_meme_memeall_tomtom.Rdata")
load("class_sel_kdm_wall_kdmtom.Rdata")
load("class_sel_fb_fball_tomtom.Rdata")
tt4<-data.frame(
	Factorbook_FBALL_Tomtom=class_sel_fb_fball_tomtom[,3],
	Meme_MEMEALL_Tomtom=class_sel_meme_memeall_tomtom[,3],
	KDM_WALL_KDMtom=class_sel_kdm_wall_kdmtom[,3]
)

p1<-plotClassification(tt1,"KDMtom",5)
ggarrange(plotlist=p1,nrow=1)
p2<-plotClassification(tt2,"KDM",6)
ggarrange(plotlist=p2,nrow=1)
p3<-plotClassification(tt3,"KDMM",6)
ggarrange(plotlist=p3,nrow=1)
p4<-plotClassification(tt4,"KDM",7)
ggarrange(plotlist=p4)



## RBP ##
setwd("/adat/Progetti/KDM/MATERIALE/RBP/CLASSIFICAZIONE")
tt1=readRDS("class_all_cisbp_mcross.rds")
tt2=readRDS("class_sel_ENCODE_mCross.rds")[,-2]
tt3=readRDS("class_sel_ENCODE_mCross.rds")[,-1]
tt4=readRDS("class_sel_ENCODE_reference.rds")

p1<-plotClassification(tt1,"KDMtom",5)
ggarrange(plotlist=p1,nrow=1)
p2<-plotClassification(tt2,"KDMtom",6)
ggarrange(plotlist=p2,nrow=1)
p3<-plotClassification(tt3,"KDMtom",6)
ggarrange(plotlist=p3,nrow=1)
p4<-plotClassification(tt4,"KDM",7)
ggarrange(plotlist=p4,nrow=1)





