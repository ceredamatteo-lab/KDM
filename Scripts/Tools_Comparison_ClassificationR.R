library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}


load("Rdata/classification.Rdata")

# 1 KDM vs STREME ----
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
pdf("Figure/Classification_01_Classes.pdf",height = 7,width = 10)
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


# 2 Classificationc vs nPeaks ----
exp_info=readRDS("Rdata/exp_info.rds")

class_2$nTrainPeaks=exp_info$nTrainPeaks

pdf("Figure/Classification_02_Classes_vs_nPeaks.pdf",height = 7,width = 10)
ggplot(class_2,aes(x=class,y=nTrainPeaks,fill=class))+geom_boxplot(notch = T)+facet_grid(~method)+scale_fill_manual(values = colors)+
  scale_y_log10()+theme_bw() +
  theme(axis.title.x = element_blank(),plot.title = element_text(size = 10),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  geom_label(data=class_2%>%count(class,method),aes(label=n),y=log10(max(class_2$nTrainPeaks)),size=3,vjust=-1,fill="white")
dev.off()

# 3 Classificationc vs AUCs ----
aucs=readRDS("Rdata/AUC_KDM.rds")%>%select(-Train_AUC,-Train_N.Motifs,-Test_AUC,-Test_N.Motifs)%>%
  rename(Train_AUC=Train_AUC_Windows,Train_N.Motifs=Train_N.Motifs_Windows,Test_AUC=Test_AUC_Windows,Test_N.Motifs=Test_N.Motifs_Windows)%>%
  mutate(ratio=log10(Train_AUC/Test_AUC))
res=merge(class_2%>%filter(method=="KDM")%>%mutate(Experiment=exp_info$Experiment),aucs%>%select(Experiment,Train_AUC,Test_AUC,ratio),by="Experiment")

method="KDM"
pdf("Figure/Classification_03_Classes_vs_AUCs.pdf",height = 5,width = 10)
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
