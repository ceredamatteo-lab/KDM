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
} else {
  load(file = "/adat/Progetti/KDM/MATERIALE/RBP/ENCODE_eCLIP_DATASET.Rdata")
  res=readRDS("/adat/Progetti/KDM/MATERIALE/RBP/pvalues_annovar2.rds")%>%filter(Annovar!=".")%>%mutate(Sign=P.val<=th,Score=-log10(P.val))
}  
#res=readRDS("/Users/tbecchi/Desktop/CLIP_maps/pvalues_annovar.rds")%>%filter(Annovar!=".")%>%mutate(Sign=P.val<=th,Score=-log10(P.val))
info_cl=rr%>%select(Exp,Cluster,N.motifs)%>%unique()

print(length(unique(res$Exp)))
print(length(unique(res$Cluster)))
print(mean(info_cl$N.motifs))

n <- length(unique(res$Annovar))
mycols <- colorRampPalette(brewer.pal(9, "Set1"))(n)

ggplot(res,aes(x=Score,fill=Annovar))+geom_density(alpha=1)+facet_wrap(~Annovar,ncol = 1,scales = "free_y")+
  scale_x_log10()+geom_vline(xintercept = -log10(th))+theme_classic()+
  scale_fill_manual(values = mycols)

res2=dcast( rr,formula = Exp+Cluster~Annovar,value.var = "P.val")
p_fish=c()
for(i in 1:nrow(res2)){
  x <- res2[i, -c(1,2)][ !is.na(res2[i, -c(1,2)]) ]
  if(0 %in% x){
    #nx=min(x[x!=0])
    x[x==0]=.Machine$double.xmin
  }
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
      if(0 %in% x){
        x[x==0]=.Machine$double.xmin
      }
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
info=merge(exp_info%>%filter(Experiment %in% res2$Exp)%>%select(Experiment,target,cell),
           res2%>%count(Exp)%>%rename(Experiment=1,n.clusters=2),by="Experiment")
rownames(info)=info$Experiment

col_fun <- colorRamp2(
  breaks = c(min(P_method, na.rm = TRUE), max(P_method, na.rm = TRUE)),
  colors = c("lightblue", "blue")
)

cell_vec <- info[rev(names(o2)), "cell"]
N_vec <- info[rev(names(o2)), "n.clusters"]
cell_types <- unique(cell_vec)
cell_colors <- structure(c("orchid","forestgreen"),names = cell_types)

row_ha <- rowAnnotation(
  Cell = cell_vec,
  N.clusters = N_vec,
  col = list(
    Cell = cell_colors,
    N.clusters = colorRamp2(c(min(N_vec, na.rm = TRUE), max(N_vec, na.rm = TRUE)),c("moccasin", "navajowhite4") )
  ),
  annotation_width = unit(c(4, 4), "mm")
)
P_method=P_method[rev(names(o2)),rev(names(o))]
row_names_heatmap <- info[rev(names(o2)), "target"]
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

pdf("/Users/tbecchi/Desktop/CLIP_maps/annovar_fisher_method.pdf",width = 6,height = 20)
ht
dev.off()

subset(exp_info,cell=="HepG2")$Experiment











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