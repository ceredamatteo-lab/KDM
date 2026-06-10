library(ggplot2)
library(ggpubr)
library(dplyr)

if(Sys.info()['nodename'] == "Matteos-MacBook-Air.local"){
  setwd("/Users/tbecchi/Desktop/repository/KDM/")
  data_folder="Rdata/Paper_Figure/"
  plot_folder="Figure/Paper_Figure/"
  data.info="Rdata/Paper_Figure/ENCODE_eCLIP_DATASET.Rdata"
} else{
  setwd("/adat/Progetti/KDM/MATERIALE/")
  data_folder="RBP/Paper_Figure/"
  plot_folder="PARTI_FIGURE/"
  data.info="RBP/ENCODE_eCLIP_DATASET.Rdata"
}

ic=readRDS(paste0(data_folder,"IC_DNA_RNA.rds"))%>%mutate(reference=ifelse(type=="DNA","HOCOMOCO","mCrossBase"))


p1=ggplot(ic,aes(x=reference,y=average.IC))+geom_boxplot(notch = T,width=0.5,outlier.size = 0.5)+theme_classic()+
  stat_compare_means(label.x = 1.5,hjust=0.5,label = "p.format")+xlab("TYPE")+ylab("Average PWM IC per position")
p2=ggplot(ic,aes(x=reference,y=KDM_IC))+geom_boxplot(notch = T,width=0.5,outlier.size = 0.5)+theme_classic()+
  stat_compare_means(label.x = 1.5,hjust=0.5,label = "p.format")+xlab("TYPE")+ylab("KDM IC")

pdf(paste0(plot_folder,"IC_boxplot.pdf"),width = 5,height = 5)
ggarrange(p1,p2,nrow = 1,align = "hv")
dev.off()


p1=ggplot(subset(ic,type=="DNA"),aes(x=average.IC,y=KDM_IC))+geom_point(size=0.5)+theme_classic()+
  stat_cor()+xlab("HOCOMOCO PWM average IC per position")+ylab("HOCOMOCO KDM IC")#+geom_smooth(method = "lm")
p2=ggplot(subset(ic,type=="RNA"),aes(x=average.IC,y=KDM_IC))+geom_point(size=0.5)+theme_classic()+
  stat_cor()+xlab("mCrossBase PWM average IC per position")+ylab("mCrossBase KDM IC")#+geom_smooth(method = "lm")

pdf(paste0(plot_folder,"IC_scatterplot.pdf"),width = 7,height = 5)
ggarrange(p1,p2,nrow = 1,align = "hv")
dev.off()

