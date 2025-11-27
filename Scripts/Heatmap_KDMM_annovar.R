library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
th=0.05

if(Sys.info()["user"]=="tbecchi"){
  load(file = "/Users/tbecchi/Desktop/repository/KDM/Rdata/Paper_Figure/ENCODE_eCLIP_DATASET.Rdata")
  res=readRDS("/Users/tbecchi/Desktop/CLIP_maps/pvalues_annovar2.rds")%>%filter(Annovar!=".")%>%mutate(Sign=P.val<=th,Score=-log10(P.val))
  plot_folder="/Users/tbecchi/Desktop/CLIP_maps/"
  ic=readRDS("/Users/tbecchi/Desktop/CLIP_maps/ICscore.rds")
} else {
  load(file = "/adat/Progetti/KDM/MATERIALE/RBP/ENCODE_eCLIP_DATASET.Rdata")
  res=readRDS("/adat/Progetti/KDM/MATERIALE/RBP/pvalues_annovar2.rds")%>%filter(Annovar!=".")%>%mutate(Sign=P.val<=th,Score=-log10(P.val))
  plot_folder="/adat/Progetti/KDM/MATERIALE/RBP/"
  ic=readRDS("/adat/Progetti/KDM/MATERIALE/RBP/CLIP_maps/ICscore.rds")
}  
ic=ic%>%dplyr::group_by(rbp,cell)%>%dplyr::summarise(IC.mean=mean(icscore),IC.sd=sd(icscore))
exp_info=merge(exp_info%>%dplyr::mutate(ID=paste0(target,"_",cell)),
               ic%>%dplyr::mutate(ID=paste0(rbp,"_",cell))%>%dplyr::select(ID,IC.mean),
               by="ID",all.x = TRUE)
info_cl=res%>%select(Exp,Cluster,N.motifs)%>%unique()

res_d=dcast(res,Exp+Cluster~Annovar,value.var = "P.val")

print(length(unique(res$Exp)))
print(length(unique(res$Cluster)))
print(mean(info_cl$N.motifs))

n <- length(unique(res$Annovar))
mycols <- colorRampPalette(brewer.pal(9, "Set1"))(n)

ggplot(res,aes(x=Score,fill=Annovar))+geom_density(alpha=1)+facet_wrap(~Annovar,ncol = 1,scales = "free_y")+
  scale_x_log10()+geom_vline(xintercept = -log10(th))+theme_classic()+
  scale_fill_manual(values = mycols)

res2=dcast( res,formula = Exp+Cluster~Annovar,value.var = "P.val")
p_fish=c()
for(i in 1:nrow(res2)){
  x <- res2[i, -c(1,2)][ !is.na(res2[i, -c(1,2)]) ]
  if(0 %in% x){x[x==0]=.Machine$double.xmin}
  X2 <- -2 * sum(log(x))
  p_fisher <- pchisq(X2, df = 2*length(x),lower.tail = FALSE)
  p_fish=c(p_fish,p_fisher)
}
res2$Fisher.Method=p_fish
res2=res2%>%filter(Fisher.Method<=th)


print(length(unique(res2$Exp)))
print(length(unique(res2$Cluster)))
print(mean(subset(info_cl,Cluster%in%res2$Cluster)$N.motifs))

P_method=rbind()
for( exp in unique(res2$Exp)){
  p_fish=c()
  tmp=subset(res2,Exp==exp)
  for(k in 3:15){
    x=tmp[,k]
    x=x[!is.na(x)]
    if (length(x)==0){
      p_fisher=NA
    } else{
      if(0 %in% x){x[x==0]=.Machine$double.xmin}
      X2 <- -2 * sum(log(x))
      p_fisher <- pchisq(X2, df = 2*length(x),lower.tail = FALSE)
    }
    p_fish=c(p_fish,p_fisher)
  }
  P_method=rbind(P_method,p_fish)
}
colnames(P_method)=colnames(res2)[3:15]
rownames(P_method)=unique(res2$Exp)
P_method[P_method>th]=NA
rem=which(colSums(is.na(P_method))==nrow(P_method))
P_method=P_method[,-rem]
P_method[!is.na(P_method) & P_method == 0]=.Machine$double.xmin
P_method=-log10(P_method)
o=colSums(!is.na(P_method))
o=o[order(o)]
o2=rowSums(!is.na(P_method))
o2=o2[order(o2)]

o2=info%>%arrange(nTrainPeaks)%>%pull(Experiment)

tmp=info_cl%>%filter(Cluster%in%res2$Cluster)%>%group_by(Exp)%>%summarise(N.Clusters=n(),Average.Motifs=mean(N.motifs))%>%rename(Experiment=1)
info=merge(exp_info%>%filter(Experiment %in% res2$Exp),tmp,by="Experiment")
rownames(info)=info$Experiment

col_fun <- colorRamp2(
  breaks = c(min(P_method, na.rm = TRUE), max(P_method, na.rm = TRUE)),
  colors = c("lightblue", "blue")
)

cell_vec <- info[rev(names(o2)), "cell"]
clus_vec <- info[rev(names(o2)), "N.Clusters"]
mot_vec <- info[rev(names(o2)), "Average.Motifs"]
seq_vec <- log10(info[rev(names(o2)), "nTrainPeaks"])

cell_types <- unique(cell_vec)
cell_colors <- structure(c("orchid","forestgreen"),names = cell_types)

P_method=P_method[rev(names(o2)),rev(names(o))]
row_names_heatmap <- info[rev(names(o2)), "target"]

row_ha <- rowAnnotation(
  Cell = cell_vec,
  N.Clusters = clus_vec,
  Avg.Motifs=mot_vec,
  `log10\n(N.Regions)`=seq_vec,
  col = list(
    Cell = cell_colors,
    N.Clusters = colorRamp2(c(min(clus_vec, na.rm = TRUE), max(clus_vec, na.rm = TRUE)),c("moccasin", "navajowhite4") ),
    Avg.Motifs = colorRamp2(c(min(mot_vec, na.rm = TRUE), max(mot_vec, na.rm = TRUE)),c("pink", "firebrick") ),
    `log10\n(N.Regions)` = colorRamp2(c(min(seq_vec, na.rm = TRUE), max(seq_vec, na.rm = TRUE)),c("mistyrose1", "mistyrose4") )
  ),
  annotation_width = unit(c(4, 4), "mm"),annotation_name_gp = gpar(fontsize = 8)
)

ht=Heatmap(
  P_method,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_labels = row_names_heatmap,
  row_names_gp = gpar(fontsize = 8),  
  column_names_gp = gpar(fontsize = 8),
  col = col_fun,
  na_col = "white",
  name = "-log10\n(Fisher's\nmethod)",
  
  left_annotation = row_ha, 
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "white", fill = fill))
  }
)
ht

pdf(paste0(plot_folder,"annovar_fisher_method.pdf"),width = 6,height = 20)
ht
dev.off()

subset(exp_info,cell=="HepG2")$Experiment


info2=res2
info2$n_sign <- rowSums(info2[, 3:15] < th, na.rm = TRUE)
info2=info2%>%dplyr::select(Exp,Cluster,Fisher.Method,n_sign)
info2=merge(info2,info_cl[,2:3],by="Cluster")
info2=merge(info2%>%dplyr::rename(Experiment=2),info[,c("Experiment","nTrainPeaks","IC.mean")],by="Experiment",all.x = TRUE)
info2$Fisher.Method[info2$Fisher.Method==0]=.Machine$double.xmin
info2$log_fisher=-log10(info2$Fisher.Method)
info2=info2%>%dplyr::rename(`Significative classes`=n_sign,`Number of motifs`=N.motifs,`Average IC`=IC.mean,`Number of regions`=nTrainPeaks,`-log10(Fisher)`=log_fisher)

corrs=data.frame()
for(k in 4:7){
  for (w in (k+1):8){
    a1=info2[,k]
    a2=info2[,w]
    t=cor.test(a1,a2)
    corrs=rbind(corrs,data.frame(v1=colnames(info2)[k],v2=colnames(info2)[w],corr=t$estimate,pval=t$p.value))
  }
}
corrs=corrs%>%mutate(Sign= pval<=th,v1=factor(v1,levels=colnames(info2)[4:8]),v2=factor(v2,levels=colnames(info2)[4:8]))
p1=ggplot(corrs,aes(x=v1,y=v2,fill=corr,alpha=Sign))+
  geom_tile(col="grey60")+scale_fill_gradient2(low = "navy",high = "firebrick",mid = "white",midpoint = 0)+
  scale_alpha_manual(values = c(0,1), guide = "none")+coord_equal()+theme_minimal()+theme(axis.title = element_blank(),panel.grid = element_blank())+
  geom_text(aes(label=round(corr,2)))+ggtitle("")+ggtitle("Cluster correlations")

info3=info2%>%dplyr::group_by(Experiment)%>%dplyr::summarise(`Average\nSignificative classes`=mean(`Significative classes`),`Average\nNumber of motifs`=mean(`Number of motifs`),
                                                             `Number of regions`=unique(`Number of regions`),`Average IC`=unique(`Average IC`),
                                                             `Average\n-log10(Fisher)`=mean(`-log10(Fisher)`),`Number of clusters`=n())
corrs=data.frame()
for(k in 2:6){
  for (w in (k+1):7){
    a1=data.frame(info3)[,k]
    a2=data.frame(info3)[,w]
    t=cor.test(a1,a2)
    corrs=rbind(corrs,data.frame(v1=colnames(info3)[k],v2=colnames(info3)[w],corr=t$estimate,pval=t$p.value))
  }
}
corrs=corrs%>%mutate(Sign= pval<=th,v1=factor(v1,levels=colnames(info3)[2:7]),v2=factor(v2,levels=colnames(info3)[2:7]))
p2=ggplot(corrs,aes(x=v1,y=v2,fill=corr,alpha=Sign))+
  geom_tile(col="grey60")+scale_fill_gradient2(low = "navy",high = "firebrick",mid = "white",midpoint = 0)+
  scale_alpha_manual(values = c(0,1), guide = "none")+coord_equal()+theme_minimal()+theme(axis.title = element_blank(),panel.grid = element_blank())+
  geom_text(aes(label=round(corr,2)))+ggtitle("Experiment correlations")

ggarrange(p1,p2,nrow = 1,align = "h")

length(unique(res$Exp))
n=res%>%count(Exp,Cluster)%>%count(Exp)
ggplot(n,aes(x=n))+geom_histogram(col="black")+
  theme_classic()+labs(x="Number of clusters")


tmp=res%>%group_by(Exp,Cluster)%>%summarise(Tot.Sign=sum(Sign))%>%group_by(Exp)%>%summarise(n=sum(Tot.Sign>0))
ggplot(tmp,aes(x=n))+geom_histogram(col="black")+
  theme_classic()+labs(x="Number of significant clusters")



o=res %>% filter(Sign)%>%count(Annovar)%>%arrange(n)%>%pull(Annovar)
ggplot(res %>% filter(Sign)%>%mutate(Annovar=factor(Annovar,levels=o)),aes(y = Annovar, fill = Annovar)) +
  geom_bar(width = 0.8) +
  geom_text(stat = "count",aes(label = ..count.., x = ..count..),  size = 3,hjust=-0.5) +
  scale_fill_manual(values = mycols) +
  theme_classic()+scale_x_continuous(expand = c(0,0,0,100))+labs(x="Count")

tmp=res%>%group_by(Exp,Annovar)%>%summarise(Tot.Sign=sum(Sign))
o=tmp %>% filter(Tot.Sign>0)%>%group_by(Annovar)%>%summarise(n=n())%>%arrange(n)%>%pull(Annovar)
ggplot(tmp %>% filter(Tot.Sign>0)%>%mutate(Annovar=factor(Annovar,levels=o)),aes(y = Annovar, fill = Annovar)) +
  geom_bar(width = 0.8) +
  geom_text(stat = "count",aes(label = ..count.., x = ..count..),  size = 3,hjust=-0.5) +
  scale_fill_manual(values = mycols) +
  theme_classic()+scale_x_continuous(expand = c(0,0,0,100))+labs(x="Count")




for 
