options(warn = -1)
suppressMessages(library(rkdMotifs))
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))

library(dplyr)
library(rkdMotifs)
library(stringr)
library(glmnet)
library(parallel)
library(optparse)


setwd("/adat/Progetti/RBP/RNA_exons")
source("/adat/Progetti/KDM//FACTORBOOK2/RTX/funzioni_importanti.rtx")
source("Scripts/Splicing_analysis_functions.R")

wd<-getwd()
ref<-paste0(wd,"/Reference")




option_list = list(
  make_option(c("-r", "--rbp"), type="character", default=NULL,
              help="RBP name", metavar="character"),
  make_option(c("-c", "--cell_line"), type="character", default=NULL,
              help="Cell line", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Genome (hg19 or hg38)", metavar="character"),
  make_option(c("-L", "--window"), type="double", default=NULL,
              help="Intronic window", metavar="double"),
  make_option(c("-a", "--alpha"), type="double", default=1,
              help="Alpha of regression", metavar="double")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

rbp<-opt$rbp
message(paste0("RBP: ",rbp))
cell_line<-opt$cell_line
message(paste0("Cell line: ",cell_line))
genome<-opt$genome
message(paste0("Genome: ",genome))
L<-opt$window
message(paste0("Intronic window: ",L))
alpha<-opt$alpha
message(paste0("Alpha for regression: ",alpha))


length<-6 #MRMs length
gaps<-2 #Number of gaps inside MRMs
mot<-70 #Number of motifs
hwin<-15 
flanks<-(L/2)+hwin+1
sign_th<-0.01 # Significance threshold to reduce the number of motifs in the final visualization

# 1. EXONS ANNOTATION + SUMMARY PLOT WITH EXONS INFO

fold<-list.files(paste0(wd,"/RNAseq/",cell_line,"/"))[grep(paste0(rbp,"-"),list.files(paste0(wd,"/RNAseq/",cell_line,"/")))]
exons<-read.table(paste0(wd,"/RNAseq/",cell_line,"/",fold,"/MATS_output/SE.MATS.JunctionCountOnly.txt"),header=TRUE)%>%
mutate(upstreamI=exonStart_0base-upstreamEE,downstreamI=downstreamES-exonEnd,exonLength=exonEnd-exonStart_0base)
exons$IncLevelDifference= -exons$IncLevelDifference

print("Differentially included: FDR <=0.1")
exons<-exons%>%mutate(Type=ifelse(FDR>0.1 | abs(IncLevelDifference)<= 0.1,"Control",ifelse(IncLevelDifference>0,"Enhanced","Silenced")))
E<-sum(exons$Type=="Enhanced")
S<-sum(exons$Type=="Silenced")
print(paste0("Enhanced: ",E," Silenced: ",S))
if(E+S>=100){plot_summary(exons,"FDR")}

if(E+S<100){
    print("Differentially included: pValue <=0.05")
    exons<-exons%>%mutate(Type=ifelse(PValue>0.05 | abs(IncLevelDifference)<= 0.1,"Control",ifelse(IncLevelDifference>0,"Enhanced","Silenced")))
    E<-sum(exons$Type=="Enhanced")
    S<-sum(exons$Type=="Silenced")
    print(paste0("Enhanced: ",E," Silenced: ",S))
    plot_summary(exons,"pValue")
}


up<-exons%>%filter(Type=="Enhanced")%>%filter(upstreamI>100,downstreamI>100)
down<-exons%>%filter(Type=="Silenced")%>%filter(upstreamI>100,downstreamI>100)
bkg<-exons%>%filter(Type=="Control")%>%filter(upstreamI>100,downstreamI>100)

# 2. CREATE INTRONICS REGIONS FROM EXONS COORDINATES

U<-do.call(rbind, unlist(apply(up,1,from_exons_to_bed), recursive = FALSE))%>%as.data.frame()
D<-do.call(rbind, unlist(apply(down,1,from_exons_to_bed), recursive = FALSE))%>%as.data.frame()
B<-do.call(rbind, unlist(apply(bkg,1,from_exons_to_bed), recursive = FALSE))%>%as.data.frame()

U_upstream_bed<-U%>%filter(region=="upstream")%>%select(-region)
U_downstream_bed<-U%>%filter(region=="downstream")%>%select(-region)
D_upstream_bed<-D%>%filter(region=="upstream")%>%select(-region)
D_downstream_bed<-D%>%filter(region=="downstream")%>%select(-region)

bkg_upstream_bed<-B%>%filter(region=="upstream")%>%select(-region)%>%arrange(abs(as.numeric(score)))
bkg_downstream_bed<-B%>%filter(region=="downstream")%>%select(-region)%>%arrange(abs(as.numeric(score)))

# 3. FREQUENCY OF EACH MRM INSIDE INTRONIC REGIONS


U_upstream_counts<-from_bed_to_counts(U_upstream_bed)
U_downstream_counts<-from_bed_to_counts(U_downstream_bed)
D_upstream_counts<-from_bed_to_counts(D_upstream_bed)
D_downstream_counts<-from_bed_to_counts(D_downstream_bed)
counts<-new("kdistr",do.call(cbind,list(U_upstream_counts,U_downstream_counts,D_upstream_counts,D_downstream_counts)),l=as.integer(length),k=as.integer(gaps),doubleStrand=FALSE,estimated=FALSE)

# 4. FACTORIZATION

W<-kdmMotifs(counts,initialization="cluster",np=0,rank=mot)

# 5. GLM MODELS with 10-FOLD CV + COEFFICIENTS ANALYSIS

coeff_U<-calculate_coeff(U_upstream_counts,U_downstream_counts,W)
coeff_D<-calculate_coeff(D_upstream_counts,D_downstream_counts,W)

# 6. OVERREPRESENTATION ANALYSYS
    
anno_U_upstream<-annotate_motif(U_upstream_bed,bkg_upstream_bed)    
anno_U_downstream<-annotate_motif(U_downstream_bed,bkg_downstream_bed)    
anno_D_upstream<-annotate_motif(D_upstream_bed,bkg_upstream_bed)    
anno_D_downstream<-annotate_motif(D_downstream_bed,bkg_downstream_bed)    

RES<-data.frame()
for(m in 1:mot){
RES<-rbind(RES,c(m,"UP","upstream",anno_U_upstream$enrichment[[m]]%>%filter(win==(L+1))%>%select(enr,enr_q.value)%>%unlist()))
RES<-rbind(RES,c(m,"UP","downstream",anno_U_downstream$enrichment[[m]]%>%filter(win==(L+1))%>%select(enr,enr_q.value)%>%unlist()))
RES<-rbind(RES,c(m,"DOWN","upstream",anno_D_upstream$enrichment[[m]]%>%filter(win==(L+1))%>%select(enr,enr_q.value)%>%unlist()))
RES<-rbind(RES,c(m,"DOWN","downstream",anno_D_downstream$enrichment[[m]]%>%filter(win==(L+1))%>%select(enr,enr_q.value)%>%unlist()))
}
colnames(RES)<-c("model","Alt","window","ENR","q.value")
RES<-RES%>%mutate(ENR=as.numeric(ENR),q.value=as.numeric(q.value))%>%mutate(Value=ifelse(q.value>sign_th,0,log10(ENR)))

# 7. MRMDs WITH A SIGNFICANT OVERREPRESNTATION IN AT LEAST ONE SUBREGION

sign<-RES%>%group_by(model)%>%summarise(n=sum(Value==0))%>%filter(n!=4)

# 8. MRMDs SCORE ACROSS REGIONS

scores<-do.call(rbind,list(
calculate_score(anno_U_upstream)%>%mutate(Alt="UP",window="upstream")%>%filter(Motif%in%sign$model),
calculate_score(anno_U_downstream)%>%mutate(Alt="UP",window="downstream")%>%filter(Motif%in%sign$model),
calculate_score(anno_D_upstream)%>%mutate(Alt="DOWN",window="upstream")%>%filter(Motif%in%sign$model),
calculate_score(anno_D_downstream)%>%mutate(Alt="DOWN",window="downstream")%>%filter(Motif%in%sign$model)
))

# 9. MRMDs to PWMs CONVERTION

pwms_converted<-kdmTopwm(W,nflank=3,cutoff=1,um=FALSE,fraction=0.99,estimate=FALSE)
  
# 10. mCROSS PWMs to MRMDs CONVERTION + COMPARISON WITH EACH MRMDS

mcross<-readRDS(paste0(ref,"/reference_",cell_line,"_ss_met.rds"))
pwms<-kdmLoadPWMSet(paste0(ref,"/mCrossBase_",cell_line,"_pwms.meme"))
colnames(mcross)<-c(names(pwms),"5ss","3ss","DRACH")
kdmtom<-kdmTom(W,mcross,threshold=0.8,best.hit=FALSE)
kdmtom$mcross<-do.call(rbind,str_split(kdmtom$target,"\\."))[,2]
kdmtom<-kdmtom[,c("query","rank","mcross")]%>%filter(rank<=5)%>%group_by(query)%>%summarise(mCross=paste0(unique(mcross),collapse="-"))

# 11. FINAL DATA PREPARATION

total_coeff<-rbind(coeff_U%>%mutate(Alt="UP"),coeff_D%>%mutate(Alt="DOWN"))%>%mutate(Type=paste0(Alt,"_",window))
total_coeff<-total_coeff%>%mutate(window=factor(window,levels=c("upstream","downstream")))
RES<-RES%>%mutate(window=factor(window,levels=c("upstream","downstream")))
scores<-scores%>%mutate(window=factor(window,levels=c("upstream","downstream")),Alt=factor(Alt,levels=c("UP","DOWN")))
kdmtom$label=gsub("-","\n",kdmtom$mCross)
kdmtom$class="1"
kdmtom$class[grep(rbp,kdmtom$mCross)]<-"2"
kdmtom$class[grep("5ss",kdmtom$mCross)]<-"3"
kdmtom$class[grep("3ss",kdmtom$mCross)]<-"4"

final_data<-list(profile=scores,test=RES,coeff=total_coeff,pwms=pwms_converted,kdmtom=kdmtom,sel=sign)

tt<-"lasso"
if(alpha==0){tt<-"ridge"}
saveRDS(final_data,paste0(wd,"/Results/",rbp,"_",cell_line,"_",genome,"_",tt,".rds"))

# 12. FINAL PLOT

plot_motif(final_data)
