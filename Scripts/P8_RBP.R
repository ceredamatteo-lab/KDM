library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
#library(ggstatsplot)


if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
}


auc=readRDS(paste0(data_folder,"p8.rds"))

auc_mean=auc[,c(1,grep("AUC",colnames(auc)))]
auc_mean=melt(auc_mean,id.vars = "Experiment")
auc_mean$variable=gsub("\\.AUC","",auc_mean$variable)

order=as.character(auc_mean%>%group_by(variable)%>%summarise(median=median(value,na.rm=T))%>%arrange(median)%>%pull(variable))
my_comparisons <- list(c("RNAProt", "KDM"), c("DeepCLIP", "KDM"), c("GraphProt", "KDM") )

colors_fill=c("gold",brewer.pal(length(unique(auc_mean$variable))-1,"Dark2"))
names(colors_fill)[1]="KDM"
names(colors_fill)[-1]=unique(auc_mean$variable)[which(unique(auc_mean$variable)!="KDM")]

pdf(paste0(plot_folder,"P8_RBP.pdf"),width = 7,height = 7)
ggplot(auc_mean%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+scale_fill_manual(values = colors_fill)+#scale_fill_brewer(palette = "Dark2")+
  geom_text(data=auc_mean%>%group_by(variable)%>%summarise(median=median(value,na.rm=T),mean=mean(value,na.rm = T)),aes(y=median,label=round(median,3)),show.legend = F,hjust=-1)+
  ggtitle("AUCs on RBP_24 with RNAProt datsets")
dev.off()
