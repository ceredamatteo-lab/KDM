library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)


if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}

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


get_jaccard=function(data,classes_sel){
  all_combos <- expand.grid(
    Var1 = classes_sel,
    Var2 = classes_sel
  )
  
  tmp <- data %>%
    count(Var1 = factor(.[[1]], levels = classes_sel),
          Var2 = factor(.[[2]], levels = classes_sel)) %>%
    right_join(all_combos, by = c("Var1","Var2")) %>%
    mutate(Freq = ifelse(is.na(n), 0, n))
  
  tmp <- tmp %>%
    group_by(Var1) %>% mutate(tot.Var1 = sum(Freq)) %>% ungroup() %>%
    group_by(Var2) %>% mutate(tot.Var2 = sum(Freq)) %>% ungroup() %>%
    mutate(Jaccard = Freq / (tot.Var1 + tot.Var2 - Freq))
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
    print(jaccard)
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



# RBP no cell line ----

tt1=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_all_cisbp_mcross.rds")
tt2=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_sel_ENCODE_mCross.rds")[,c(3,1)]
tt3=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_sel_ENCODE_mCross.rds")[,c(3,2)]
tt4=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_sel_ENCODE_reference.rds")[,c(2,1)]

p1<-plotClassification(tt1,"KDMtom",5,adjust = T)
p2<-plotClassification(tt2,"KDMtom",6,adjust = T)
p3<-plotClassification(tt3,"KDMtom",6,adjust = T)
p4<-plotClassification(tt4,"KDM",7,adjust = T)


pdf("Figure/ENCODE_Dataset/New_figure/01_kdmtom_tomtom.pdf",height = 7,width = 10)
ggarrange(plotlist=p1,nrow=1)
dev.off()

pdf("Figure/ENCODE_Dataset/New_figure/02_kdm_streme_mcross.pdf",height = 7,width = 10)
ggarrange(plotlist=p2,nrow=1)
dev.off()

pdf("Figure/ENCODE_Dataset/New_figure/03_kdmm_streme_mcross.pdf",height = 7,width = 10)
ggarrange(plotlist=p3,nrow=1)
dev.off()

pdf("Figure/ENCODE_Dataset/New_figure/04_kdm_streme_wall.pdf",height = 7,width = 10)
ggarrange(plotlist=p4,nrow=1)
dev.off()

# RBP cell line ----

tt2=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_sel_ENCODE_mCross_cellLine.rds")[,-2]
tt3=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_sel_ENCODE_mCross_cellLine.rds")[,-1]
tt4=readRDS("Rdata/ENCODE_Dataset/CLASSIFICAZIONE/class_sel_ENCODE_reference_cellLine.rds")

p2<-plotClassification(tt2,"KDMtom",6,adjust = T)
p3<-plotClassification(tt3,"KDMtom",6,adjust = T)
p4<-plotClassification(tt4,"KDM",7,adjust = T)

pdf("Figure/ENCODE_Dataset/New_figure/cell_line_02_kdm_streme_mcross.pdf",height = 7,width = 10)
ggarrange(plotlist=p2,nrow=1)
dev.off()

pdf("Figure/ENCODE_Dataset/New_figure/cell_line_03_kdmm_streme_mcross.pdf",height = 7,width = 10)
ggarrange(plotlist=p3,nrow=1)
dev.off()

pdf("Figure/ENCODE_Dataset/New_figure/cell_line_04_kdm_streme_wall.pdf",height = 7,width = 10)
ggarrange(plotlist=p4,nrow=1)
dev.off()

## classificazione filtrata per cell.line ----
tt3=cbind(tt3,exp_info)
tt4=cbind(tt4,exp_info)
p3<-plotClassification(subset(tt3,cell=="K562")[,1:2],"KDMtom",6,adjust = T)
ggarrange(plotlist=p3,nrow=1)
p4<-plotClassification(subset(tt4,cell=="K562")[,1:2],"KDM",7,adjust = T)
ggarrange(plotlist=p4,nrow=1)
