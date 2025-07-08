library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

if(Sys.info()['user'] == "tbecchi"){setwd("/Users/tbecchi/Desktop/repository/KDM/")}

annovar=readRDS("Rdata/RBP_24/RBP24_annovar_data.rds")
annovar=annovar[,-c(grep("nc",colnames(annovar)),grep("intergenic",colnames(annovar)))]
scaled <- annovar / rowSums(annovar)

info=read.delim("Tables/RBP_24_Info.tsv")
vn=readRDS("Rdata/ENCODE_Dataset/VN_clusters.rds")
info=merge(info,vn%>%filter(cell_line=="HepG2")%>%rename(PROTEIN=gene_name)%>%dplyr::select(PROTEIN,VN_Class),by="PROTEIN",all.x=T)%>%rename(ENCODE_HepG2=VN_Class)
info=merge(info,vn%>%filter(cell_line=="K562")%>%rename(PROTEIN=gene_name)%>%dplyr::select(PROTEIN,VN_Class),by="PROTEIN",all.x=T)%>%rename(ENCODE_K562=VN_Class)

info$ENCODE_K562[grep("IGF2BP",info$PROTEIN)]=unique(subset(vn[grep("IGF2BP",vn$gene_name),],cell_line=="K562")$VN_Class)
info$ENCODE_HepG2[grep("IGF2BP",info$PROTEIN)]=unique(subset(vn[grep("IGF2BP",vn$gene_name),],cell_line=="HepG2")$VN_Class)


info <- info[match(rownames(scaled), info$Experiment), ]
row_anno <- rowAnnotation(
  ENCODE.HepG2=info$ENCODE_HepG2,
  ENCODE.K562=info$ENCODE_K562,
  Method = info$METHOD,
  Cell.Line = info$CELL_LINE,
  gp = gpar(col = "black"),na_col="white",
  col=list(ENCODE.HepG2=c("Intron"="#2DAEE0","5ss"="#7D619A","Non-coding Exon"="#F8AE1A","3'UTR + CDS"="#197C3F","CDS"="#31AF47","CDS + other"="#8FC349"),
           ENCODE.K562=c("Intron"="#2DAEE0","5ss"="#7D619A","Non-coding Exon"="#F8AE1A","3'UTR + CDS"="#197C3F","CDS"="#31AF47","CDS + other"="#8FC349"),
           Method= setNames(brewer.pal(4,"Set1"),unique(info$METHOD)),
           Cell.Line= setNames(brewer.pal(3,"Set2"),unique(info$CELL_LINE))),
  gap = unit(c(1,3,1), "mm")
)

set.seed(30580)
km <- kmeans(scaled, centers = 5)
cluster_labels <- factor(km$cluster, labels = c("3'UTR+Intronic", "Intronic","3'UTR","3'UTR+CDS","CDS"))

Heatmap(
  scaled,
  name = "Scaled Proportion",
  col = colorRamp2(c(0, 1), c("white", "red")),
  rect_gp = gpar(col = "black", lwd = 0.5),
  row_split = cluster_labels,
  row_gap = unit(5, "mm"),
  right_annotation = row_anno,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", scaled[i, j]),  # formato con 2 decimali
      x = x, y = y,
      gp = gpar(fontsize = 8, col = "black")
    )
  }
)

