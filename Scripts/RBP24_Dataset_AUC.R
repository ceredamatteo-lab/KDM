library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)

if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}

# 0. RBP_23.optimization = 23 experiments, max 500 peaks ----
## Ottimizzazione dei parametri Rank, MRM Length e MRM Gaps

auc=readRDS("Rdata/RBP_24/Optimization_Parameters_AUC.rds")%>%
  mutate(Combination=paste0(Length,"-",Gaps,"-",Rank),MRM=paste0(Length,"-",Gaps))

d=auc%>%group_by(Length,Gaps,Rank)%>%summarise(Val=median(AUC),n=n())%>%
  mutate(Combination=paste0(Length,"-",Gaps,"-",Rank),MRM=paste0(Length,"-",Gaps))

pdf("Figure/RBP_24/01_Best_parameters.pdf",width = 7,height = 7)
ggplot(d%>%mutate(Rank=as.numeric(Rank)),aes(x=Rank,y=Val))+geom_hline(yintercept = subset(d,Length==6&Gaps==2&Rank==70)$Val,lty="dashed")+
  geom_line(aes(col=MRM),lwd=1)+geom_point(aes(fill=MRM),size=3,col="black",shape=21)+theme_minimal()+
  ylab("MEDIAN AUC")+scale_x_continuous(breaks = unique(as.numeric(d$Rank)))+
  scale_fill_brewer(type = "qual",palette = "Dark2")+scale_color_brewer(type = "qual",palette = "Dark2")+
  guides(fill = guide_legend( override.aes = list(size = 3)))+theme(axis.title.y = element_text(angle=0,vjust=0.5))+
  ggtitle("Median AUC across 23 optimiazion datasets\nfor each Length~Gaps combination at each Rank")
dev.off()


# 1. RBP_24_RNAPRot: 23 dataset con max 5000 picchi per ogni esperimento ----
auc=readRDS("Rdata/RBP_24/RBP24_RNAProt_AUC.rds")

auc_mean=auc[,c(1,grep("AUC",colnames(auc)))]
auc_mean=melt(auc_mean,id.vars = "Experiment")
auc_mean$variable=gsub("\\.AUC","",auc_mean$variable)

order=as.character(auc_mean%>%group_by(variable)%>%summarise(median=median(value,na.rm=T))%>%arrange(median)%>%pull(variable))
my_comparisons <- list(c("RNAProt", "KDM"), c("DeepCLIP", "KDM"), c("GraphProt", "KDM") )

pdf("Figure/RBP_24/02_RBP24_RNAProt_AUC.pdf",width = 7,height = 7)
ggplot(auc_mean%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+scale_fill_brewer(palette = "Dark2")+
  geom_text(data=auc_mean%>%group_by(variable)%>%summarise(median=median(value,na.rm=T),mean=mean(value,na.rm = T)),aes(y=median,label=round(median,3)),show.legend = F,hjust=-1)+
  ggtitle("AUCs on RBP_24 with RNAProt datsets")
dev.off()

## 1.1 AUC Distribtions vs Clusters (Annovar) ----
info=read.delim("Tables/RBP_24_Info.tsv")%>%arrange(Experiment)
auc_mean=merge(auc_mean,info%>%select(Experiment,CLUSTER),by="Experiment")

pdf("Figure/RBP_24/03_RBP24_RNAProt_AUC_Annovar_Clusters.pdf",width = 10,height = 10)
ggplot(auc_mean%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+scale_fill_brewer(palette = "Dark2")+
  ggtitle("AUCs on RBP_24 with RNAProt datsets ~ Annovar Clusters")+facet_wrap(~CLUSTER,scales = "free")
dev.off()

# 2. RBP_24_GraphProt: 23 dataset con tutti i picchi original di Graphprot ----
info=read.delim("Tables/RBP_24_Info.tsv")%>%arrange(Experiment)
auc=readRDS("Rdata/RBP_24/RBP24_GraphProt_AUC.rds")%>%mutate(N.Peaks=info$N.Peaks,Cluster=info$CLUSTER)
auc_mean=melt(auc,id.vars = c("RNA.binding.protein","N.Peaks","Cluster"))

order=as.character(auc_mean%>%group_by(variable)%>%summarise(median=median(value,na.rm=T))%>%arrange(median)%>%pull(variable))
my_comparisons <- list(c("GraphProt", "KDM"), c("deepnet.rbp", "KDM"), c("DeepBind", "KDM"),c("iDeepV", "KDM"),c("iDeepE", "KDM") 
                       ,c("CNN.LSTM.E", "KDM") ,c("ResNet.E", "KDM") ,c("Pse.SVM", "KDM"))

pdf("Figure/RBP_24/04_RBP24_GraphProt_AUC.pdf",width = 10,height = 7)
ggplot(auc_mean%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+#scale_fill_brewer(palette = "Dark2")+
  geom_text(data=auc_mean%>%group_by(variable)%>%summarise(median=median(value,na.rm=T),mean=mean(value,na.rm = T)),
            aes(y=median,label=round(median,2)),show.legend = F,hjust=-0.7)+
  ggtitle("AUCs on RBP_24 with GraphProt datsets")
dev.off()

## 2.1 Correlazine AUC ~ N.Peaks ----
pdf("Figure/RBP_24/05_RBP24_GraphProt_AUC_Correlation.pdf",width = 11,height = 11)
ggarrange(
ggplot(auc_mean,aes(x=N.Peaks,y=value,col=variable,fill=variable))+stat_cor(aes(label=paste(..r.label..,..p.label.., unique(auc_mean$variable)[..group..], sep = " - ")),
                                                                            show.legend = F,label.y = seq(1,1.1,length.out=length(unique(auc_mean$variable))))+
  geom_smooth(method='lm',alpha=0.1)+theme_bw()+xlab("Number of peaks")+ylab("AUC")+ggtitle("AUC ~ N.Peaks"),

ggplot(auc_mean,aes(x=N.Peaks,y=value,col=variable,fill=variable))+stat_cor(aes(label=paste(..r.label..,..p.label.., unique(auc_mean$variable)[..group..], sep = " - ")),
                                                                            show.legend = F,label.y = seq(1,1.1,length.out=length(unique(auc_mean$variable))))+
  geom_smooth(method='lm',alpha=0.1)+theme_bw()+xlab("Number of peaks")+ylab("AUC")+ggtitle("AUC ~ log10(N.Peaks)")+scale_x_log10(),
ncol = 1,align = "hv",common.legend = T,legend = "right")
dev.off()

## 2.2 AUC Distribtions vs Clusters (Annovar) ---
pdf("Figure/RBP_24/06_RBP24_GraphProt_AUC_Annovar_Clusters.pdf",width = 20,height = 10)
ggplot(auc_mean%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+#scale_fill_brewer(palette = "Dark2")+
  ggtitle("AUCs on RBP_24 with GraphProt datsets ~ Annovar Clusters")+facet_wrap(~Cluster,scales = "free")
dev.off()
