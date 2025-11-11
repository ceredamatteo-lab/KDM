library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
#library(ggstatsplot)


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

kdm_color="gold"
load(data.info)

# 1. COMPARISON with other methods ----

all_auc=readRDS(paste0(data_folder,"p6.rds"))
all_auc$nTrainPeaks=exp_info$nTrainPeaks
all_auc$target=exp_info$target
all_auc$cell=exp_info$cell
#all_auc=merge(all_auc%>%mutate(ID=paste0(target,"_",cell)),vn%>%select(ID,VN_Class),by="ID",all.x = T)
all_auc$N.Class=cut(all_auc$nTrainPeaks, breaks = quantile(all_auc$nTrainPeaks, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), include.lowest = TRUE, labels =  c("LOW", "MEDIUM", "HIGH", "VERY.HIGH"))

all_auc_m=melt(all_auc,id.vars = c("Experiment","nTrainPeaks","target","cell","N.Class"))
all_auc_m$value[all_auc_m$value==0]<-NA

order=as.character(all_auc_m%>%group_by(variable)%>%summarise(median=median(value,na.rm=T))%>%arrange(median)%>%pull(variable))
my_comparisons <- list( c("STREME", "KDM"), c("GraphProt", "KDM"), c("RNAProt", "KDM") )

colors_fill=c("gold",hue_pal()(length(levels(all_auc_m$variable))-1))
names(colors_fill)[1]="KDM"
names(colors_fill)[-1]=levels(all_auc_m$variable)[which(levels(all_auc_m$variable)!="KDM")]

## 2.1 AUC Distributions ----

pdf(paste0(plot_folder,"P6A_RBP.pdf"),height = 5,width = 5)
ggplot(all_auc_m%>%mutate(variable=factor(variable,levels=order)),aes(x=variable,y=value,fill=variable))+
  geom_boxplot(aes(fill=variable),width=0.5,notch=T,show.legend=F,alpha=1)+
  stat_compare_means(comparisons = my_comparisons,method="wilcox",paired=TRUE)+
  theme_bw()+theme(axis.title.x=element_blank(),axis.title.y=element_text(angle=0,vjust=0.5),panel.grid=element_blank())+
  ylab("AUC")+scale_fill_manual(values = colors_fill)
dev.off()


## 2.5 Correlation AUC~N.Peaks ----

coefs<-all_auc_m %>%group_by(variable) %>%
  do({
    #model <- lm(value ~ log10(nTrainPeaks), data = .)
    model <- lm(value ~ log10(nTrainPeaks), data = .)
    data.frame(slope = coef(model)[2], intercept = coef(model)[1])
  })

kdm<-coefs %>% filter(variable == "KDM")
graphprot<-coefs %>% filter(variable == "GraphProt")
rnaprot<-coefs %>% filter(variable == "RNAProt")

intersect_x<-function(a1, b1, a2, b2){ 
  x0 <- (b2 - b1) / (a1 - a2)
  return(x0)}

x1 = intersect_x(kdm$slope, kdm$intercept, graphprot$slope, graphprot$intercept)
x2 = intersect_x(kdm$slope, kdm$intercept, rnaprot$slope, rnaprot$intercept)

data_freq=data.frame(nTrainPeaks=all_auc%>%arrange(nTrainPeaks)%>%pull(nTrainPeaks))%>%
mutate(Percentile=seq(1,250)/250)%>%mutate(xA=abs(nTrainPeaks-10^x1),xB=abs(nTrainPeaks-10^x2))
y1=data_freq$Percentile[which(data_freq$xA==min(data_freq$xA))]
y2=data_freq$Percentile[which(data_freq$xB==min(data_freq$xB))]

pdf(paste0(plot_folder,"P6B_RBP.pdf"),height = 6,width = 8)
ggarrange(
  ggplot(all_auc_m%>%mutate(variable=factor(variable,levels=order)),aes(x=nTrainPeaks,y=value,col=variable))+stat_cor(show.legend = F)+
    geom_smooth(method='lm',alpha=0.2)+theme_bw()+theme(axis.title.x = element_blank())+
    geom_segment(x = x1, xend = x1, y = -Inf, yend = kdm$slope * x1 + kdm$intercept, linetype = "dashed", color = "black",lwd=0.3)+
    geom_segment(x = x2, xend = x2, y = -Inf, yend = kdm$slope * x2 + kdm$intercept, linetype = "dashed", color = "black",lwd=0.3)+
    ylab("AUC")+scale_color_manual(values = colors_fill)+scale_x_log10()
  ,
  ggplot(data_freq,aes(x=nTrainPeaks,y=Percentile))+
    annotate(geom = "rect",xmin = 10^x1,xmax = 10^x2,ymin=0,ymax=y2,fill="gold",alpha=0.5,lwd=0)+
    annotate(geom = "rect",xmin = 0,xmax = 10^x1,ymin=y1,ymax=y2,fill="gold",alpha=0.5,lwd=0)+
    geom_line(col="black",lwd=0.5)+theme_bw()+scale_x_log10(breaks=round(c(10,100,1000,10000,10^x1,10^x2),0))+
    scale_y_continuous(breaks=round(c(seq(0,1,length.out=5),y1,y2),2))+
    annotate("segment",y=y1,yend=y1,x=0,xend=10^x1,lty="dashed",lwd=0.3)+
    annotate("segment",y=y2,yend=y2,x=0,xend=10^x2,lty="dashed",lwd=0.3)+
    annotate("segment",y=0,yend=y1,x=10^x1,xend=10^x1,lty="dashed",lwd=0.3)+
    annotate("segment",y=0,yend=y2,x=10^x2,xend=10^x2,lty="dashed",lwd=0.3)+
    theme(panel.grid=element_blank())

  ,
  align = "hv",ncol=1,heights = c(6,4))
dev.off()

