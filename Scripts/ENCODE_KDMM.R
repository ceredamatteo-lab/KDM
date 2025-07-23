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

ggplot(data$mm,aes(x=variable,y=value))+geom_boxplot(notch = T)+facet_wrap(~Annovar,scales = "free")

mm=merge(data$mm,data.frame(data$anno_freq)%>%rename(Annovar=Var1),by="Annovar")%>%mutate(Group=paste0(Annovar,"\n(",Freq,")"))

p3=grouped_ggbetweenstats(data = mm,x = Group, y = value,grouping.var = variable,
  notch = TRUE,results.subtitle = FALSE, subtitle = NULL,centrality.plotting=F,ggtheme = ggplot2::theme_minimal(),palette = "Set2")
p4<-ggplot(data=data$um,aes(x=X1,y=X2,color=status)) + geom_point(size=0.7)+scale_color_brewer(palette="Set1")

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


grouped_ggbetweenstats(data = data$mm,x = variable, y = value,grouping.var = Annovar,
                       notch = TRUE,results.subtitle = FALSE, subtitle = NULL,ggtheme = ggplot2::theme_minimal(),palette = "Set2")

um=data$um%>%mutate(status=factor(status,levels=c("exonic","splicing","intronic","UTR5","ncRNA_exonic")))
pA=ggarrange(ggbetweenstats(data=um,x=status,y=X1,package = "RColorBrewer",palette = "Set1",results.subtitle = FALSE, subtitle = NULL,boxplot.args = list(notch=TRUE,width=0.2),point.args = list(alpha=0)),
          ggbetweenstats(data=um,x=status,y=X2,package = "RColorBrewer",palette = "Set1",results.subtitle = FALSE, subtitle = NULL,boxplot.args = list(notch=TRUE,width=0.2),point.args = list(alpha=0)),ncol =1)

pB=ggplot(um,aes(x=X1,y=X2))+geom_density_2d_filled(show.legend = F,col="black")+facet_wrap(~status, nrow=1)+theme_minimal()+theme(panel.grid = element_blank())

pC=ggplot(data=um,aes(x=X1,y=X2,color=status)) + geom_point(size=0.7,show.legend = F)+scale_color_brewer(palette="Set1")+facet_grid(~status)+theme_minimal()
ggarrange(pA,pB,pC,ncol = 1,align = "hv",heights = c(4,2,4))

tmp=subset(data$um,status%in%c("intronic","exonic"))%>%mutate(X=ifelse(X1>0,"right","left"),Y=ifelse(X2>0,"top","bottom"),Quadrante=paste0(X,"-",Y))
dcast(tmp%>%group_by(status,Quadrante)%>%summarise(n=n()),status~Quadrante)
#ggarrange(ggarrange(p1,p2,ncol=2,widths=c(1,3)),ggarrange(p3,p4,ncol=2,widths=c(3,1)),nrow=2,heights = c(4,6))

# TF ----
data=readRDS("Rdata/KDMM_plot_TF.rds")

cscale<-c25[1:nrow(data$info2)]
p1<-ggplot(data=data$info2,aes(x=name,y=`Cumulated probability`,fill=name)) + geom_bar(stat="identity",col="black") + scale_fill_manual(values=cscale)
p2<-ggplot(data=data$info3,aes(x=x,y=prob,fill=name)) + geom_bar(position="stack",stat="identity") + scale_fill_manual(values=cscale) + facet_grid(plt~.,scale="free_y")

#p2=ggplot(data=data$info3,aes(x=x,y=prob,color=name, fill=name))+geom_area(aes(group=name), alpha=0.8, position="stack",col="black")+facet_grid(plt~.,scale="free_y")+theme_minimal()+ scale_fill_manual(values=cscale)

p3=grouped_ggbetweenstats(data = data$mm,x = status, y = value,grouping.var = Var2,
                          notch = TRUE,results.subtitle = FALSE, subtitle = NULL,centrality.plotting=F,ggtheme = ggplot2::theme_minimal(),palette = "Set2")
p4<-ggplot(data=data$um,aes(x=X1,y=X2,color=status)) + geom_point(size=0.7)+scale_color_brewer(palette="Set1")


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




