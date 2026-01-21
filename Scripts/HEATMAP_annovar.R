library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)

th=0.05
dist_pairwise=function(x) {
  n <- nrow(x)
  D <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      good <- !is.na(x[i, ]) & !is.na(x[j, ])
      
      if (sum(good) == 0) {
        d <- 1
      } else {
        d <- sqrt(sum((x[i, good] - x[j, good])^2))
      }
      
      D[i, j] <- D[j, i] <- d
    }
  }
  return(as.dist(D))
}

if(Sys.info()["user"]=="tbecchi"){
  load(file = "/Users/tbecchi/Desktop/repository/KDM/Rdata/Paper_Figure/ENCODE_eCLIP_DATASET.Rdata")
  res=readRDS("/Users/tbecchi/Desktop/CLIP_maps/pvalues_annovar2.rds")%>%filter(Annovar!=".")%>%mutate(Sign=P.val<=th,Score=-log10(P.val))
  plot_folder="/Users/tbecchi/Desktop/CLIP_maps/"
  ic=readRDS("/Users/tbecchi/Desktop/CLIP_maps/ICscore.rds")
  gene_categories=readRDS("/Users/tbecchi/Desktop/repository/Pipelines/GeneCategories/Rdata/gene_categories.20251120.rds")
} else {
  load(file = "/adat/Progetti/KDM/MATERIALE/RBP/ENCODE_eCLIP_DATASET.Rdata")
  res=readRDS("/adat/Progetti/KDM/MATERIALE/RBP/pvalues_annovar2.rds")%>%filter(Annovar!=".")%>%mutate(Sign=P.val<=th,Score=-log10(P.val))
  plot_folder="/adat/Progetti/KDM/MATERIALE/RBP/"
  ic=readRDS("/adat/Progetti/KDM/MATERIALE/RBP/CLIP_maps/ICscore.rds")
} 


gc=dcast(gene_categories%>%mutate(value=TRUE),gene_name~category,value.var = "value",fun.aggregate = function(x) TRUE, fill = FALSE)%>%
  rename(target=gene_name)
gc=cbind(gc,
         dcast(gene_categories%>%mutate(value=TRUE),gene_name~subtype,value.var = "value",fun.aggregate = function(x) TRUE, fill = FALSE)%>%
           select(m6a_eraser,m6a_reader,m6a_regulator,m6a_writer,polyA,splicing))
ic=ic%>%dplyr::group_by(rbp,cell)%>%dplyr::summarise(IC.mean=mean(icscore),IC.sd=sd(icscore))%>%ungroup()
exp_info=merge(exp_info%>%dplyr::mutate(ID=paste0(target,"_",cell)),
               ic%>%dplyr::mutate(ID=paste0(rbp,"_",cell))%>%dplyr::select(ID,IC.mean),
               by="ID",all.x = TRUE) %>% group_by(target, cell) %>% arrange(desc(nTrainPeaks)) %>%
  mutate(TOP = row_number() == 1) %>% ungroup()
info_cl=res%>%select(Exp,Cluster,N.motifs)%>%unique()

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

P_all=P_method
tmp=res2[,c(3:15)]
freq=apply(tmp, 2, function(x) sum(is.na(x)))/nrow(tmp)
rem=names(freq[freq>0.5])

tmp=info_cl%>%filter(Cluster%in%res2$Cluster)%>%group_by(Exp)%>%summarise(N.Clusters=n(),Average.Motifs=mean(N.motifs))%>%rename(Experiment=1)
info=merge(exp_info%>%filter(Experiment %in% res2$Exp),tmp,by="Experiment")
info=merge(info,gc,by="target",all.x = TRUE)
rownames(info)=info$Experiment

for( cell_line in c("HepG2","K562")){
  sel=subset(info,cell==cell_line & !is.na(IC.mean) & TOP & nTrainPeaks>=500)$Experiment
  #sel=subset(info,cell==cell_line & TOP)$Experiment
  
  P_method=P_all[rownames(P_all)%in%sel,!colnames(P_all)%in%rem]
  
  #info_sel=info[info$cell==cell_line,]
  #gc_sel=gc[rownames(P_method),]
  info_sel=info[rownames(P_method),]
  idx=make.unique(info_sel$target)
  rownames(P_method)=idx
  rownames(info_sel)=idx
  
  d_row <- dist_pairwise(P_method)
  d_col <- dist_pairwise(t(P_method))
  clust_row <- hclust(d_row)
  clust_col <- hclust(d_col)
  
  col_fun <- colorRamp2(
    breaks = c(min(P_method, na.rm = TRUE), max(P_method, na.rm = TRUE)),
    #colors = c("lightblue", "blue")
    colors = viridis(2)
  )
  clus_vec=info_sel[, "N.Clusters"]
  mot_vec=info_sel[, "Average.Motifs"]
  seq_vec=log10(info_sel[, "nTrainPeaks"])
  ic_vec=info_sel[, "IC.mean"]
  
  row_ha <- rowAnnotation(
    #Cell = cell_vec,
    EM=info_sel[,"EM"],
    ETM=info_sel[,"ETM"],m6a_eraser=info_sel[,"m6a_eraser"],m6a_reader=info_sel[,"m6a_reader"],m6a_regulator=info_sel[,"m6a_regulator"],m6a_writer=info_sel[,"m6a_writer"],
    PT=info_sel[,"PT"],polyA=info_sel[,"polyA"],
    TF=info_sel[,"TF"],
    Splicing=info_sel[,"splicing"],
    N.Clusters = clus_vec,
    Avg.Motifs=mot_vec,
    `log10\n(N.Regions)`=seq_vec,
    Avg.IC=ic_vec,
    col = list(
      Cell = cell_colors,
      N.Clusters = colorRamp2(c(min(clus_vec, na.rm = TRUE), max(clus_vec, na.rm = TRUE)),c("moccasin", "navajowhite4") ),
      Avg.Motifs = colorRamp2(c(min(mot_vec, na.rm = TRUE), max(mot_vec, na.rm = TRUE)),c("pink", "firebrick") ),
      `log10\n(N.Regions)` = colorRamp2(c(min(seq_vec, na.rm = TRUE), max(seq_vec, na.rm = TRUE)),c("mistyrose1", "mistyrose4") ),
      Avg.IC = colorRamp2(c(min(ic_vec, na.rm = TRUE), max(ic_vec, na.rm = TRUE)),c("lightgreen", "darkgreen") ),
      EM=c("TRUE"="black","FALSE"="white"),
      ETM=c("TRUE"="black","FALSE"="white"),
      PT=c("TRUE"="black","FALSE"="white"),
      TF=c("TRUE"="black","FALSE"="white"),Splicing=c("TRUE"="black","FALSE"="white"),
      m6a_eraser=c("TRUE"="black","FALSE"="white"),m6a_reader=c("TRUE"="black","FALSE"="white"),
      m6a_regulator=c("TRUE"="black","FALSE"="white"),m6a_writer=c("TRUE"="black","FALSE"="white"),polyA=c("TRUE"="black","FALSE"="white")
    ),
    annotation_width = unit(c(1, 1), "mm"),annotation_name_gp = gpar(fontsize = 8),
    gp = gpar(col = "black"),
    annotation_legend_param = list(border = TRUE,labels_gp = gpar(fontsize = 7),title_gp = gpar(fontsize = 7, fontface = "bold")),
    gap = unit(c(3,1,1,1,1, 3, 1, 3, 3, 3, 1, 1, 1, 1), "mm"),na_col = "white"
  )
  set.seed(30580)
  ht=Heatmap(
    P_method,
    cluster_rows = clust_row,
    cluster_columns = clust_col,
    row_names_side = "left",
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
  pdf(paste0(plot_folder,"annovar_fisher_method_subset_",cell_line,".pdf"),width = 9,height = 9)
  #pdf(paste0(plot_folder,"annovar_fisher_method_ALL_",cell_line,".pdf"),width = 9,height = 12)
  print(ht)
  dev.off()
}





info2=res2
info2$n_sign <- rowSums(info2[, 3:15] < th, na.rm = TRUE)
info2=info2%>%dplyr::select(Exp,Cluster,Fisher.Method,n_sign)
info2=merge(info2,info_cl[,2:3],by="Cluster")
info2=merge(info2%>%dplyr::rename(Experiment=2),info[,c("Experiment","cell","nTrainPeaks","IC.mean")],by="Experiment",all.x = TRUE)
info2$Fisher.Method[info2$Fisher.Method==0]=.Machine$double.xmin
info2$log_fisher=-log10(info2$Fisher.Method)
info2=info2%>%dplyr::rename(`Significative classes`=n_sign,`Number of motifs`=N.motifs,`Average IC`=IC.mean,`Number of regions`=nTrainPeaks,`-log10(Fisher)`=log_fisher)%>%
  relocate(cell, .after = last_col())
bkp=info2

for( cell_line in c("HepG2","K562")){
  info2=subset(bkp,cell==cell_line)
  corrs=data.frame()
  for(k in 4:7){
    for (w in (k+1):8){
      a1=info2[,k]
      a2=info2[,w]
      t=cor.test(a1,a2)
      corrs=rbind(corrs,data.frame(v1=colnames(info2)[k],v2=colnames(info2)[w],corr=t$estimate,pval=t$p.value))
    }
  }
  corrs=corrs%>%mutate(Sign= pval<=th,v1=factor(v1,levels=colnames(info2)[4:8]),v2=factor(v2,levels=colnames(info2)[4:8]))%>%
    mutate(
      p_stars = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE         ~ ""
      )
    )
  
  pdf(paste0(plot_folder,"annovar_fisher_method_correlation_",cell_line,"_cluster.pdf"),width = 9,height = 9)
  print(create_heatmap(corrs))
  dev.off()
  
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
  corrs=corrs%>%mutate(Sign= pval<=th,v1=factor(v1,levels=colnames(info3)[2:7]),v2=factor(v2,levels=colnames(info3)[2:7]))%>%
    mutate(
      p_stars = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE         ~ ""
      )
    )
  
  pdf(paste0(plot_folder,"annovar_fisher_method_correlation_",cell_line,"_experiments.pdf"),width = 9,height = 9)
  print(create_heatmap(corrs))
  dev.off()
  
}


p1=ggplot(corrs,aes(x=v1,y=v2,fill=corr,alpha=Sign))+
  geom_tile(col="black")+scale_fill_gradient2(low = "navy",high = "firebrick",mid = "white",midpoint = 0)+
  scale_alpha_manual(values = c(0,1), guide = "none")+coord_equal()+theme_minimal()+theme(axis.title = element_blank(),panel.grid = element_blank())+
  #geom_text(aes(label=round(corr,2)))
  geom_text(aes(label = paste0(round(corr, 2), "\n", p_stars)),size = 4)+ggtitle("")+ggtitle("Cluster correlations")


corrs=data.frame()
for(k in 2:6){
  for (w in (k+1):7){
    a1=data.frame(info3)[,k]
    a2=data.frame(info3)[,w]
    t=cor.test(a1,a2)
    corrs=rbind(corrs,data.frame(v1=colnames(info3)[k],v2=colnames(info3)[w],corr=t$estimate,pval=t$p.value))
  }
}
corrs=corrs%>%mutate(Sign= pval<=th,v1=factor(v1,levels=colnames(info3)[2:7]),v2=factor(v2,levels=colnames(info3)[2:7]))%>%
  mutate(
    p_stars = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  )
p2=ggplot(corrs,aes(x=v1,y=v2,fill=corr,alpha=Sign))+
  geom_tile(col="black")+scale_fill_gradient2(low = "navy",high = "firebrick",mid = "white",midpoint = 0)+
  scale_alpha_manual(values = c(0,1), guide = "none")+coord_equal()+theme_minimal()+theme(axis.title = element_blank(),panel.grid = element_blank())+
  #geom_text(aes(label=round(corr,2)))
  geom_text(aes(label = paste0(round(corr, 2), "\n", p_stars)),size = 4)+ggtitle("Experiment correlations")

ggarrange(p1,p2,nrow = 1,align = "h")


create_heatmap=function(corr){
  corrs=corrs%>%mutate(val=ifelse(pval<=th,corr,NA),lab=paste0(round(corr, 2), "\n", p_stars))
  tmp=dcast(corrs,v2~v1,value.var = "val")
  rownames(tmp)=tmp[,1]
  tmp2=dcast(corrs,v2~v1,value.var = "lab")
  tmp2_mat <- as.matrix(tmp2[,-1])
  rownames(tmp2_mat) <- tmp2[,1]
  tmp2_mat <- tmp2_mat[rownames(tmp), colnames(tmp[,-1])]
  if(min(tmp[,-1],na.rm = T)<0){
    col_fun <- colorRamp2(
      c(min(tmp[,-1], na.rm = TRUE), 0, max(tmp[,-1], na.rm = TRUE)),
      c("blue", "white", "red"))
  } else{
    col_fun <- colorRamp2(
      c( 0, max(tmp[,-1], na.rm = TRUE)),
      c("white", "red"))
  }
    
  return(Heatmap(
    tmp[,-1],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = rev(rownames(tmp)),
    col = col_fun,na_col = "white",row_names_side = "left",
    #rect_gp = gpar(col = "white", lwd = 0.5),
    name = "correlation",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (j <=i) {
        grid.rect(
          x, y, w, h,
          gp = gpar(fill = NA, col = "black", lwd = 0.5)
        )
      }
      if (!is.na(tmp[i, j + 1])) {
        grid.text(
          tmp2_mat[i, j],
          x, y,
          gp = gpar(fontsize = 12)
        )
      }
    }
  )
  )
}




