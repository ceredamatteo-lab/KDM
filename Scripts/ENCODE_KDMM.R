library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggstatsplot)

if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}


c25 <- c("dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A",
  "#FF7F00",
  "gold1",
  "skyblue2", "#FB9A99",
  "palegreen2",
  "#CAB2D6",
  "#FDBF6F",
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown", "black" 
)


# RBP ----

data=readRDS("Rdata/KDMM_plot.rds")

cscale<-c25[1:nrow(data$info2)]
p1<-ggplot(data=data$info2,aes(x=name,y=`Cumulated probability`,fill=name)) + geom_bar(stat="identity",col="black") + scale_fill_manual(values=cscale)
p2<-ggplot(data=data$info3,aes(x=x,y=prob,fill=name)) + geom_bar(position="stack",stat="identity") + scale_fill_manual(values=cscale) + facet_grid(plt~.,scale="free_y")

#p2=ggplot(data=data$info3,aes(x=x,y=prob,color=name, fill=name))+geom_area(aes(group=name), alpha=0.8, position="stack",col="black")+facet_grid(plt~.,scale="free_y")+theme_minimal()+ scale_fill_manual(values=cscale)

mm=merge(data$mm,data.frame(data$anno_freq)%>%rename(Annovar=Var1),by="Annovar")%>%mutate(Group=paste0(Annovar,"\n(",Freq,")"))

p3=grouped_ggbetweenstats(data = mm,x = Group, y = value,grouping.var = variable,
  notch = TRUE,results.subtitle = FALSE, subtitle = NULL,centrality.plotting=F,ggtheme = ggplot2::theme_minimal(),palette = "Set2")
p4<-ggplot(data=data$um,aes(x=X1,y=X2,color=status)) + geom_point(size=0.7)+scale_color_brewer(palette="Set2")

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_01.pdf",width = 7,height = 7)
p1
dev.off()

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_02.pdf",width = 10,height = 5)
p2
dev.off()

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_03.pdf",width = 15,height = 5)
p3
dev.off()

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_04.pdf",width = 7,height = 7)
p4
dev.off()


#ggarrange(ggarrange(p1,p2,ncol=2,widths=c(1,3)),ggarrange(p3,p4,ncol=2,widths=c(3,1)),nrow=2,heights = c(4,6))

# TF ----
data=readRDS("Rdata/KDMM_plot_TF.rds")

cscale<-c25[1:nrow(data$info2)]
p1<-ggplot(data=data$info2,aes(x=name,y=`Cumulated probability`,fill=name)) + geom_bar(stat="identity",col="black") + scale_fill_manual(values=cscale)
p2<-ggplot(data=data$info3,aes(x=x,y=prob,fill=name)) + geom_bar(position="stack",stat="identity") + scale_fill_manual(values=cscale) + facet_grid(plt~.,scale="free_y")

#p2=ggplot(data=data$info3,aes(x=x,y=prob,color=name, fill=name))+geom_area(aes(group=name), alpha=0.8, position="stack",col="black")+facet_grid(plt~.,scale="free_y")+theme_minimal()+ scale_fill_manual(values=cscale)

p3=grouped_ggbetweenstats(data = data$mm,x = status, y = value,grouping.var = Var2,
                          notch = TRUE,results.subtitle = FALSE, subtitle = NULL,centrality.plotting=F,ggtheme = ggplot2::theme_minimal(),palette = "Set2")
p4<-ggplot(data=data$um,aes(x=X1,y=X2,color=status)) + geom_point(size=0.7)+scale_color_brewer(palette="Set2")


ggarrange(ggarrange(p1,p2,ncol=2,widths=c(1,3)),ggarrange(p3,p4,ncol=2,widths=c(3,1)),nrow=2,heights = c(4,6))

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_01.pdf",width = 7,height = 7)
p1
dev.off()

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_02.pdf",width = 10,height = 5)
p2
dev.off()

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_03.pdf",width = 15,height = 5)
p3
dev.off()

pdf("Figure/ENCODE_Dataset/KDMM/KDMM_04.pdf",width = 7,height = 7)
p4
dev.off()




