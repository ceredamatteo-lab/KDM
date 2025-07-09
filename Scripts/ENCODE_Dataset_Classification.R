library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}

# A: NO cell line info ----

load("Rdata/ENCODE_Dataset/classification.Rdata")

## 1 KDM vs STREME ----
class_1=do.call(rbind,list(classification$kdm_mcross%>%mutate(method="KDM"),
                           classification$kdmr_mcross%>%mutate(method="KDM_reduced"),
                           classification$streme_mcross%>%mutate(method="STREME")))%>%filter(class!="Not in reference")%>%
  mutate(class=factor(class,levels=rev(c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif"))))

class_2=do.call(rbind,list(classification$kdm_wall%>%mutate(method="KDM"),
                           classification$streme_pwmall%>%mutate(method="STREME")))%>%
  mutate(class=factor(class,levels=rev(c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference"))))


colors=c("grey90","grey80",brewer.pal(5,"Reds"))
names(colors)=levels(class_2$class)

class_1_summary=class_1 %>%count(method, class)%>%group_by(method)%>%mutate(prop = n / sum(n))%>%
  arrange(method, desc(class))%>%mutate(cprop=cumsum(prop),label = paste0(n,"(",round(cprop * 100), "%)"),ypos = cumsum(prop) - prop / 2)
class_2_summary=class_2 %>%count(method, class)%>%group_by(method)%>%mutate(prop = n / sum(n))%>%
  arrange(method, desc(class))%>%mutate(cprop=cumsum(prop),label = paste0(n,"(",round(cprop * 100), "%)"),ypos = cumsum(prop) - prop / 2)
pdf("Figure/ENCODE_Dataset/Classification_01_Classes.pdf",height = 7,width = 10)
ggarrange(
  ggplot(class_1_summary, aes(x = method, y = prop, fill = class)) +
    geom_bar(stat = "identity", col = "black", position = "stack") +
    geom_text(aes(label = label, y = ypos), size = 3, color = "black") +
    theme_minimal() +
    theme(axis.title = element_blank(),panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Classification against mCross"),
  ggplot(class_2_summary, aes(x = method, y = prop, fill = class)) +
    geom_bar(stat = "identity", col = "black", position = "stack") +
    geom_text(aes(label = label, y = ypos), size = 3, color = "black") +
    theme_minimal() +
    theme(axis.title = element_blank(),panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Classification train vs test"),
  align = "hv",nrow = 1)
dev.off()


## 2 Classificationc vs nPeaks ----
exp_info=readRDS("Rdata/ENCODE_Dataset/exp_info.rds")

class_2$nTrainPeaks=exp_info$nTrainPeaks

pdf("Figure/ENCODE_Dataset/Classification_02_Classes_vs_nPeaks.pdf",height = 7,width = 10)
ggplot(class_2,aes(x=class,y=nTrainPeaks,fill=class))+geom_boxplot(notch = T)+facet_grid(~method)+scale_fill_manual(values = colors)+
  scale_y_log10()+theme_bw() +
  theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_label(data=class_2%>%count(class,method),aes(label=n),y=log10(max(class_2$nTrainPeaks)),size=3,vjust=-1,fill="white")
dev.off()

## 3 Classificationc vs AUCs ----
aucs=readRDS("Rdata/ENCODE_Dataset/AUC_KDM.rds")%>%select(-Train_AUC,-Train_N.Motifs,-Test_AUC,-Test_N.Motifs)%>%
  rename(Train_AUC=Train_AUC_Windows,Train_N.Motifs=Train_N.Motifs_Windows,Test_AUC=Test_AUC_Windows,Test_N.Motifs=Test_N.Motifs_Windows)%>%
  mutate(ratio=log10(Train_AUC/Test_AUC))
res=merge(class_2%>%filter(method=="KDM")%>%mutate(Experiment=exp_info$Experiment),aucs%>%select(Experiment,Train_AUC,Test_AUC,ratio),by="Experiment")

method="KDM"
pdf("Figure/ENCODE_Dataset/Classification_03_Classes_vs_AUCs.pdf",height = 5,width = 10)
ggarrange(
  ggplot(res,aes(x=class,y=Train_AUC,fill=class))+geom_boxplot(notch = T)+scale_fill_manual(values = colors)+
    theme_bw() +theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
                      axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    ggtitle(paste0("AUC TRAIN on TEST - ",method)),
  
  ggplot(res,aes(x=class,y=Test_AUC,fill=class))+geom_boxplot(notch = T)+scale_fill_manual(values = colors)+
    theme_bw() +theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
                      axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    ggtitle(paste0("AUC TEST on TRAIN - ",method)),
  
  ggplot(res,aes(x=class,y=ratio,fill=class))+geom_boxplot(notch = T)+scale_fill_manual(values = colors)+
    theme_bw() +theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
                      axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    ggtitle(paste0("TRAIN/TEST AUC ratio - ",method)),
  align = "hv",nrow = 1,common.legend = T,legend = "right")

dev.off()



# B: Cell line info ----

classification=readRDS("Rdata/ENCODE_Dataset/classification_new.rds")
  
## 1 KDM vs STREME ----
class_1=do.call(rbind,list(classification$kdm_mcross%>%mutate(method="KDM"),
                           classification$kdmr_mcross%>%mutate(method="KDM_reduced"),
                           classification$streme_mcross%>%mutate(method="STREME")))%>%filter(class!="Not in reference")%>%
  mutate(class=factor(class,levels=rev(c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif"))))

class_2=do.call(rbind,list(classification$kdm_wall%>%mutate(method="KDM"),
                           classification$streme_pwmall%>%mutate(method="STREME")))%>%
  mutate(class=factor(class,levels=rev(c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference"))))


colors=c("grey90","grey80",brewer.pal(5,"Reds"))
names(colors)=levels(class_2$class)

class_1_summary=class_1 %>%count(method, class)%>%group_by(method)%>%mutate(prop = n / sum(n))%>%
  arrange(method, desc(class))%>%mutate(cprop=cumsum(prop),label = paste0(n,"(",round(cprop * 100), "%)"),ypos = cumsum(prop) - prop / 2)
class_2_summary=class_2 %>%count(method, class)%>%group_by(method)%>%mutate(prop = n / sum(n))%>%
  arrange(method, desc(class))%>%mutate(cprop=cumsum(prop),label = paste0(n,"(",round(cprop * 100), "%)"),ypos = cumsum(prop) - prop / 2)
pdf("Figure/ENCODE_Dataset/Classification_01_Classes_CellLine.pdf",height = 7,width = 10)
ggarrange(
  ggplot(class_1_summary, aes(x = method, y = prop, fill = class)) +
    geom_bar(stat = "identity", col = "black", position = "stack") +
    geom_text(aes(label = label, y = ypos), size = 3, color = "black") +
    theme_minimal() +
    theme(axis.title = element_blank(),panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Classification against mCross"),
  ggplot(class_2_summary, aes(x = method, y = prop, fill = class)) +
    geom_bar(stat = "identity", col = "black", position = "stack") +
    geom_text(aes(label = label, y = ypos), size = 3, color = "black") +
    theme_minimal() +
    theme(axis.title = element_blank(),panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggtitle("Classification train vs test"),
  align = "hv",nrow = 1)
dev.off()


## 2 Classificationc vs nPeaks ----
exp_info=readRDS("Rdata/ENCODE_Dataset/exp_info.rds")

class_2$nTrainPeaks=exp_info$nTrainPeaks

pdf("Figure/ENCODE_Dataset/Classification_02_Classes_vs_nPeaks_CellLine.pdf",height = 7,width = 10)
ggplot(class_2,aes(x=class,y=nTrainPeaks,fill=class))+geom_boxplot(notch = T)+facet_grid(~method)+scale_fill_manual(values = colors)+
  scale_y_log10()+theme_bw() +
  theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_label(data=class_2%>%count(class,method),aes(label=n),y=log10(max(class_2$nTrainPeaks)),size=3,vjust=-1,fill="white")
dev.off()

## 3 Classificationc vs AUCs ----
aucs=readRDS("Rdata/ENCODE_Dataset/AUC_KDM.rds")%>%select(-Train_AUC,-Train_N.Motifs,-Test_AUC,-Test_N.Motifs)%>%
  rename(Train_AUC=Train_AUC_Windows,Train_N.Motifs=Train_N.Motifs_Windows,Test_AUC=Test_AUC_Windows,Test_N.Motifs=Test_N.Motifs_Windows)%>%
  mutate(ratio=log10(Train_AUC/Test_AUC))
res=merge(class_2%>%filter(method=="KDM")%>%mutate(Experiment=exp_info$Experiment),aucs%>%select(Experiment,Train_AUC,Test_AUC,ratio),by="Experiment")

method="KDM"
pdf("Figure/ENCODE_Dataset/Classification_03_Classes_vs_AUCs_CellLine.pdf",height = 5,width = 10)
ggarrange(
  ggplot(res,aes(x=class,y=Train_AUC,fill=class))+geom_boxplot(notch = T)+scale_fill_manual(values = colors)+
    theme_bw() +theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
                      axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    ggtitle(paste0("AUC TRAIN on TEST - ",method)),
  
  ggplot(res,aes(x=class,y=Test_AUC,fill=class))+geom_boxplot(notch = T)+scale_fill_manual(values = colors)+
    theme_bw() +theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
                      axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    ggtitle(paste0("AUC TEST on TRAIN - ",method)),
  
  ggplot(res,aes(x=class,y=ratio,fill=class))+geom_boxplot(notch = T)+scale_fill_manual(values = colors)+
    theme_bw() +theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),
                      axis.text.x = element_blank(),axis.ticks.x = element_blank())+
    ggtitle(paste0("TRAIN/TEST AUC ratio - ",method)),
  align = "hv",nrow = 1,common.legend = T,legend = "right")

dev.off()


## 4 Jaccard Index ----
return_jacccard=function(a,b){
  l1=length(a)
  l2=length(b)
  l3=sum(a%in%b)
  return(l3/(l1+l2-l3))
}

result_jaccard=rbind()
for(c in unique(class_1$class)){
  sel=subset(class_1,class==c)%>%mutate(ID=paste0(target,"_",cell))
  kdm=unique(subset(sel,method=="KDM")$ID)
  kdmr=unique(subset(sel,method=="KDM_reduced")$ID)
  str=unique(subset(sel,method=="STREME")$ID)
  result_jaccard=rbind(result_jaccard,c(c,"KDM-KDMR",return_jacccard(kdm,kdmr)))
  result_jaccard=rbind(result_jaccard,c(c,"KDM-STREME",return_jacccard(kdm,str)))
  result_jaccard=rbind(result_jaccard,c(c,"KDMR-STREME",return_jacccard(kdmr,str)))
}
result_jaccard=data.frame(result_jaccard)
colnames(result_jaccard)=c("class","pair","jaccard")
result_jaccard=result_jaccard%>%mutate(jaccard=as.numeric(jaccard),class=factor(class,levels=c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference")))

ggplot(result_jaccard,aes(x=class,y=jaccard,fill=class))+geom_bar(stat="identity",col="black")+
  facet_grid(~pair)+theme_bw() +
  theme(axis.title = element_blank(),panel.grid = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values = colors)+geom_text(aes(label=round(jaccard,2)),vjust=-1)


result_jaccard=rbind()
for(c in unique(class_2$class)){
  sel=subset(class_2,class==c)%>%mutate(ID=paste0(target,"_",cell))
  kdm=unique(subset(sel,method=="KDM")$ID)
  str=unique(subset(sel,method=="STREME")$ID)
  result_jaccard=rbind(result_jaccard,c(c,"KDM-STREME",return_jacccard(kdm,str)))
}
result_jaccard=data.frame(result_jaccard)
colnames(result_jaccard)=c("class","pair","jaccard")
result_jaccard=result_jaccard%>%mutate(jaccard=as.numeric(jaccard),class=factor(class,levels=c("(0,1]","(1,3]","(3,10]","(10,Inf]","No match","No central motif","Not in reference")))

ggplot(result_jaccard,aes(x=class,y=jaccard,fill=class))+geom_bar(stat="identity",col="black")+
  facet_grid(~pair)+theme_bw() +
  theme(axis.title = element_blank(),panel.grid = element_blank(),plot.title = element_text(size = 10)) +
  scale_fill_manual(values = colors)+geom_text(aes(label=round(jaccard,2)),vjust=-1)

tmp=subset(class_1,method!="KDM")%>%mutate(ID=paste0(target,"_",cell))


