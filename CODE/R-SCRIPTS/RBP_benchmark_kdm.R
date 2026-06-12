library(rkdMotifs)
#library(stringr)
library(glmnet)
library(ROCR)

setwd("/adat/Progetti/KDM/MATERIALE/ALL/")

# KDM ----

s<-"lambda.min"
set=1
experiment=do.call(rbind,strsplit(list.files(paste0("DATASETS/RBP_BENCHMARK/set",set,"_hg38_bed_fasta/"),pattern="positives.fa"),"\\."))[,1]
ov_win=10
length=6
gaps=2
rank=100

setwd("/adat/Progetti/KDM/MATERIALE/ALL/DATASETS/RBP_BENCHMARK/")
result=data.frame(experiment=experiment,AUROC=NA)
for(i in 1:length(experiment)){
    print(i)
    exp=experiment[i]
    pos_seq<-scan(paste0("set",set,"_hg38_bed_fasta/",exp,".positives.fa"),what="char")
    pos_ids=gsub(">","",pos_seq[seq(1,length(pos_seq)-1,2)])
    pos_seq=pos_seq[seq(2,length(pos_seq),2)]
    pos_seq=gsub("U", "T", pos_seq)

    neg_seq<-scan(paste0("set",set,"_hg38_bed_fasta/",exp,".negatives.fa"),what="char")
    neg_ids=gsub(">","",neg_seq[seq(1,length(neg_seq)-1,2)])
    neg_seq=neg_seq[seq(2,length(neg_seq),2)]
    neg_seq=gsub("U", "T", neg_seq)
    
    total_win=unique(nchar(c(pos_seq,neg_seq)))
    win_kmer= (total_win+2*ov_win)/3
    
    test<-scan(paste0("set",set,"_hold_out_comparison_test_ids/",exp,".test_ids"),what="char")
    test_pos_x=which(pos_ids %in% test)
    test_neg_x=which(neg_ids %in% test)
    
    pos_train=pos_seq[-test_pos_x]
    neg_train=neg_seq[-test_neg_x]
    train_seq=c(pos_train,neg_train)
    
    upstream_train<-substr(train_seq, 1 , win_kmer+1)
    centre_train<-substr(train_seq, ceiling(nchar(train_seq)/2)-(win_kmer/2), floor(nchar(train_seq)/2)+(win_kmer/2)+1)
    downstream_train<-substr(train_seq, nchar(train_seq)-win_kmer+1 , nchar(train_seq))
    
    pos_test=pos_seq[test_pos_x]
    neg_test=neg_seq[test_neg_x]
    test_seq=c(pos_test,neg_test)
    
    upstream_test<-substr(test_seq, 1 , win_kmer+1)
    centre_test<-substr(test_seq, ceiling(nchar(test_seq)/2)-(win_kmer/2), floor(nchar(test_seq)/2)+(win_kmer/2)+1)
    downstream_test<-substr(test_seq, nchar(test_seq)-win_kmer+1 , nchar(test_seq))
    
    label=c(rep(1,length(pos_train)),rep(0,length(neg_train)))
    test_label=c(rep(1,length(pos_test)),rep(0,length(neg_test)))

    freq_W<-kdmDistr(train_seq,length,gaps,double.strand=FALSE,strict=FALSE)
    W<-kdmMotifs(freq_W,initialization="cluster",np=0,rank=rank)
    saveRDS(W,paste0("W_KDM/",exp,".rds"))

    train_features=do.call(cbind,list(kdmFeatures(W,upstream_train),kdmFeatures(W,centre_train),kdmFeatures(W,downstream_train)))
    
    set.seed(123)
    mm<-cv.glmnet(y=label,x=train_features, family="binomial", type.measure="auc",alpha=1,nfold=10)
    
    test_features=do.call(cbind,list(kdmFeatures(W,upstream_test),kdmFeatures(W,centre_test),kdmFeatures(W,downstream_test)))
    
    auc=performance(prediction(labels=test_label,predictions=predict(mm,test_features,type="link",s=s)),"auc")@y.values[[1]]
    
    result$AUROC[i]=auc 
}

library(readxl)
setwd("/adat/Progetti/KDM/MATERIALE/ALL/")
other=read_xlsx("RBP_PREDICTION/RBP_BENCHMARK/Supplementary_tables_S1-S7.xlsx",sheet="Table S5",skip=4)[-1,]
other=other[-24,]
colnames(other)[1]="Experiment"
result[c(13,14),]=result[c(14,13),]
other$KDM=result$AUROC
other$DeepRAM=as.numeric(other$DeepRAM)
other$RNAProt=as.numeric(other$RNAProt)

# Graphrot ----
other$GraphProt=NA
for(i in 1:nrow(other)){
    exp=other$Experiment[i]
    if(exp=="PARCLIP_MOV10"){exp="PARCLIP_MOV10_Sievers"}
    test<-scan(paste0("DATASETS/RBP_BENCHMARK/set",set,"_hold_out_comparison_test_ids/",exp,".test_ids"),what="char")
    pos=sum(grepl("pos",test))
    neg=sum(grepl("neg",test))

    pred=read.table(paste0("RBP_PREDICTION/RBP_BENCHMARK/graphprot/PREDICTION/",exp,".predictions"))
    auc=performance(prediction(labels=c(rep(1,pos),rep(0,neg)),predictions=pred[,3]),"auc")@y.values[[1]]
    other$GraphProt[i]=auc
}

saveRDS("RBP_PREDICTION/RBP_BENCHMARK/AUC_RNAProt23_RBP.rds")


