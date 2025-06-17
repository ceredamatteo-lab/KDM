library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)

if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}


# 1. AUCs with KDM ----

aucs=readRDS("Rdata/ENCODE_Dataset/AUC_KDM.rds")%>%mutate(ID=paste0(target,"_",cell))%>%
  select(-Train_AUC,-Train_N.Motifs,-Test_AUC,-Test_N.Motifs)%>%
  rename(Train_AUC=Train_AUC_Windows,Train_N.Motifs=Train_N.Motifs_Windows,Test_AUC=Test_AUC_Windows,Test_N.Motifs=Test_N.Motifs_Windows)%>%
  mutate(ratio=log10(Train_AUC/Test_AUC))

## 1.1 TRAIN and TEST ----

pdf("Figure/ENCODE_Dataset/AUC_01_TRAIN_and_TEST.pdf",width = 10,height = 5)
ggarrange(ggplot(aucs,aes(x=Train_AUC,y=Test_AUC))+geom_point()+stat_cor()+ggtitle("TRAIN AUC vs TEST AUC")+theme_bw(),
          ggplot(aucs,aes(x=nTrainPeaks,y=abs(ratio)))+
            geom_point()+scale_x_log10()+
            stat_cor()+ggtitle("TRAIN/TEST RATIO vs N.PEAKS")+theme_bw()+
            ylab(expression(abs(log[10](frac(TRAIN, TEST)))))+theme(axis.title.y=element_text(angle=0,vjust=0.5))+xlab("Input peaks"),
          align = "hv",nrow = 1)
dev.off()

## 1.2 TRAIN_AUC SUMMARY -----

aucs=aucs%>%rename(AUC=Train_AUC,N.Motifs=Train_N.Motifs)
ggtexttable(
  t(data.frame(MIN=min(aucs$AUC),MAX=max(aucs$AUC),MEDIAN=median(aucs$AUC),MAD=mad(aucs$AUC),AVERAGE=mean(aucs$AUC),SD=sd(aucs$AUC))%>%
      mutate(across(where(is.numeric), ~ round(.x, 4))))%>%data.frame()%>%dplyr::rename(Val=1),
  theme = ttheme("light", base_size = 14))

pdf("Figure/ENCODE_Dataset/AUC_02_AUC_SUMMARY.pdf",width = 15,height = 8)
ggarrange(
  ggplot(aucs,aes(x=cell,y=AUC,fill=cell))+geom_boxplot(notch = T,width=0.3)+stat_compare_means(label.x = 1.5,hjust=0.5,label = "p.format")+
    scale_fill_manual(values = colors_cell)+theme_bw()+theme(panel.grid = element_blank())+
    geom_text(data=aucs%>%group_by(cell)%>%summarise(n=n()),aes(label=n),y=1)+ggtitle("AUCS ~ cell.line"),
  
  ggplot(aucs,aes(x=cell,y=nTrainPeaks,fill=cell))+geom_boxplot(notch = T,width=0.3)+stat_compare_means(label.x = 1.5,hjust=0.5,label = "p.format")+
    scale_fill_manual(values = colors_cell)+theme_bw()+theme(panel.grid = element_blank())+
    #geom_text(data=aucs%>%group_by(cell)%>%summarise(n=n()),aes(label=n),y=1)+
    scale_y_log10()+ggtitle("N.Peaks ~ cell.line")+ylab("Input peaks"),
  
  ggplot(aucs,aes(x=nTrainPeaks,y=AUC))+geom_point(aes(fill=cell),col="black",shape=21,size=3)+
    scale_x_log10()+stat_cor()+theme_bw()+scale_fill_manual(values = colors_cell)+ggtitle("AUCS ~ N.Peaks")+xlab("Input peaks"),
  
  ggplot(aucs,aes(x=nTrainPeaks,y=N.Motifs))+geom_point(aes(fill=cell),col="black",shape=21,size=3)+scale_x_log10()+stat_cor()+theme_bw()+scale_fill_manual(values = colors_cell)+ggtitle("N.Motifs ~ N.Peaks")+xlab("Input peaks"),
  ggplot(aucs,aes(x=N.Motifs,y=AUC))+geom_point(aes(fill=cell),col="black",shape=21,size=3)+stat_cor()+theme_bw()+scale_fill_manual(values = colors_cell)+ggtitle("AUCS ~ N.Motifs"),
  
  ggplot(dcast(aucs,target~cell,value.var = "AUC",fun.aggregate = mean),aes(x=HepG2,y=K562))+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")+geom_point(size=3)+stat_cor()+theme_bw()+
    ggtitle(paste0("RBP on both cell.line (n=",nrow(dcast(aucs,target~cell,value.var = "AUC",fun.aggregate = mean)%>%filter(!is.na(HepG2),!is.na(K562))),")")), 
  align ="hv"
)

dev.off()

## 1.3 VAN NOSTRAND CLUSTERS----
colors_vn=c("Intron"="#2DAEE0","5ss"="#7D619A","Non-coding Exon"="#F8AE1A","3'UTR + CDS"="#197C3F","CDS"="#31AF47","CDS + other"="#8FC349")
vn=readRDS("Rdata/ENCODE_Dataset/VN_clusters.rds")
aucs=merge(aucs,vn%>%select(ID,VN_Class),by="ID",all.x = T)

pdf("Figure/ENCODE_Dataset/AUC_03_VanNostrand_Clusters.pdf",width = 15,height = 5)
ggarrange(
  ggplot(aucs%>%filter(!is.na(VN_Class)),aes(x=VN_Class,y=AUC,fill=VN_Class))+geom_boxplot(alpha=0.8,notch = T)+
    stat_compare_means(ref.group = "Intron",label = "p.format")+
    scale_fill_manual(values = colors_vn)+
    theme_bw()+theme(axis.text.x=element_blank(), panel.grid= element_blank(),axis.ticks.x = element_blank())+ggtitle("AUC vs VN.Class"),
  
  ggplot(aucs%>%filter(!is.na(VN_Class)),aes(x=VN_Class,y=nTrainPeaks,fill=VN_Class))+geom_boxplot(alpha=0.8,notch = T)+
    stat_compare_means(ref.group = "Intron",label = "p.format")+
    scale_fill_manual(values = colors_vn)+
    theme_bw()+theme(axis.text.x=element_blank(),panel.grid= element_blank(),axis.ticks.x = element_blank())+scale_y_log10()+
    ggtitle("N.Peaks vs VN.Class")+ylab("Input peaks"),
  common.legend = T,legend = "right",
  align = "hv",nrow = 1
)
dev.off()


# 2. COMPARISON with other methods ----

all_auc=readRDS("Rdata/ENCODE_Dataset/AUC_Comparison.rds")
exp_info=readRDS("Rdata/ENCODE_Dataset/exp_info.rds")
all_auc$nTrainPeaks=exp_info$nTrainPeaks
all_auc$target=exp_info$target
all_auc$cell=exp_info$cell
all_auc=merge(all_auc%>%mutate(ID=paste0(target,"_",cell)),vn%>%select(ID,VN_Class),by="ID",all.x = T)
all_auc$N.Class=cut(all_auc$nTrainPeaks, breaks = quantile(all_auc$nTrainPeaks, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), include.lowest = TRUE, labels =  c("LOW", "MEDIUM", "HIGH", "VERY.HIGH"))

all_auc_m=melt(all_auc,id.vars = c("Experiment","ID","nTrainPeaks","target","cell","VN_Class","N.Class"))
all_auc_m$value[all_auc_m$value==0]<-NA

order=as.character(all_auc_m%>%group_by(variable)%>%summarise(median=median(value,na.rm=T))%>%arrange(median)%>%pull(variable))
my_comparisons <- list( c("STREME", "KDM"), c("GraphProt", "KDM"), c("RNAProt", "KDM") )

## 2.1 AUC Distributions ----

pdf("Figure/ENCODE_Dataset/AUC_04_AUCs_Boxplot.pdf",height = 5,width = 5)
ggplot(all_auc_m%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")
dev.off()

## 2.2 AUC Distributions ~ CLUSTERS----
med<-all_auc_m%>%filter(!is.na(VN_Class))%>%group_by(variable, VN_Class)%>%summarise(median = median(value, na.rm = TRUE), .groups = "drop") %>%group_by(VN_Class) %>%filter(median == max(median)) %>%ungroup()
pdf("Figure/ENCODE_Dataset/AUC_05_AUCs_Boxplot_Clusters.pdf",height = 10,width = 10)
ggplot(all_auc_m%>%filter(!is.na(VN_Class))%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_hline(data = med, aes(yintercept = median),linetype = "dashed")+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+facet_wrap(~VN_Class,scales = "free")+ggtitle("AUC distributions ~ VN_Classes")
dev.off()

## 2.3 AUC Distributions ~ N.PEAKS----
med<-all_auc_m%>%group_by(variable, N.Class)%>%summarise(median = median(value, na.rm = TRUE), .groups = "drop") %>%group_by(N.Class) %>%filter(median == max(median)) %>%ungroup()
pdf("Figure/ENCODE_Dataset/AUC_06_AUCs_Boxplot_nPeaks.pdf",height = 10,width = 10)
ggplot(all_auc_m%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_hline(data = med, aes(yintercept = median),linetype = "dashed")+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+facet_wrap(~N.Class,scales = "free")+ggtitle("AUC distributions ~ N.Peaks (1st,2nd,3rd,4th quartile)")
dev.off()

## 2.4 AUC Distributions ~ Cell----
med<-all_auc_m%>%group_by(variable, cell)%>%summarise(median = median(value, na.rm = TRUE), .groups = "drop") %>%group_by(cell) %>%filter(median == max(median)) %>%ungroup()
pdf("Figure/ENCODE_Dataset/AUC_07_AUCs_Boxplot_Cell.pdf",height = 5,width = 7)
ggplot(all_auc_m%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_hline(data = med, aes(yintercept = median),linetype = "dashed")+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+facet_wrap(~cell,scales = "free")+ggtitle("AUC distributions ~ Cell.line")
dev.off()
