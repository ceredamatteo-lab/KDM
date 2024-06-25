start<-Sys.time()

library(optparse)

option_list = list(
  make_option(c("-p", "--project"), type="character", default=NULL,
              help="Project name", metavar="character"),
  make_option(c("-l", "--length"), type="double", default=NULL,
              help="Kmer length", metavar="double"),
  make_option(c("-g", "--gaps"), type="double", default=0,
              help="Number of gaps (default: 0)", metavar="double"),
  make_option(c("-i", "--initialization"), type="character", default="cluster",
              help="Inizialization method (default: cluster)", metavar="character"),
  make_option(c("-a", "--alpha"), type="double", default=1,
              help="Elastic net mixing parameter (default: 1)", metavar="double"),
  make_option(c("-C", "--combine"), type="logical", default=FALSE,
              help="Combine negative and positive counts as input (default: FALSE)", metavar="logical"),
  make_option(c("-m", "--motifs"), type="double", default=NULL,
              help="Number of motifs", metavar="double")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

project<-opt$project
message(paste0("Project: ",project))
length<-opt$length
message(paste0("Kmer length: ",length))
gaps<-opt$gaps
message(paste0("Gaps: ",gaps))
init<-opt$initialization
message(paste0("Initialization: ",init))
alpha<-opt$alpha
message(paste0("Elastic net mixing parameter: ",alpha))
combine<-opt$combine
message(paste0("Combine negative and positive counts as input: ",combine))
mot<-opt$motifs
message(paste0("Number of motifs: ",mot))


suppressMessages(library(glmnet))
suppressMessages(library(rkdMotifs))
library(rkdMotifs)
library(glmnet)
library(ROCR)
library(stringr)

if( Sys.info()['nodename']%in%c('aulonocara.bp.lnf.it')){

  DATA_DIR="/home/tommaso.becchi/Documents/GraphProt/"
  #CGB_DIR=''
  #GIT_LOCAL_DIR=''
  #CGB_SHARED='~/Dropbox (HuGeF)/SPLICING/CGB_shared/'

  }

setwd(DATA_DIR)
exp<-readRDS(paste0("Rdata/",project,"_cv.rds"))

pos<-exp$pos_seq
neg<-exp$neg_seq
pos_ids<-exp$pos_ids
neg_ids<-exp$neg_ids
lista_ids<-exp$cv

create_windows<-function(input){
    left<-substr(input, floor(nchar(input)/2)-46 , floor(nchar(input)/2)-16)
    centre<-substr(input, floor(nchar(input)/2)-15, floor(nchar(input)/2)+15)
    right<-substr(input, floor(nchar(input)/2)+16 , floor(nchar(input)/2)+46)
    
    C1<-kdmDistr(left,length,gaps,double.strand=FALSE,strict=TRUE)
    CC<-kdmDistr(centre,length,gaps,double.strand=FALSE,strict=TRUE)
    C3<-kdmDistr(right,length,gaps,double.strand=FALSE,strict=TRUE)
    
    return(list(left=C1,centre=CC,right=C3))
    #return(new("kdistr",cbind(C1,CC,C3),l=as.integer(length),k=as.integer(gaps),doubleStrand=FALSE,estimated=FALSE))
    
}

pos_counts<-create_windows(pos)
neg_counts<-create_windows(neg)


pos_CT_init<-new("kdistr",cbind(pos_counts$left[,1:1000] ,pos_counts$centre[,1:1000],pos_counts$right[,1:1000]),l=as.integer(length),k=as.integer(gaps),doubleStrand=FALSE,estimated=FALSE)


pos_CT<-new("kdistr",cbind(pos_counts$left,pos_counts$centre,pos_counts$right),l=as.integer(length),k=as.integer(gaps),doubleStrand=FALSE,estimated=FALSE)
neg_CT<-new("kdistr",cbind(neg_counts$left,neg_counts$centre,neg_counts$right),l=as.integer(length),k=as.integer(gaps),doubleStrand=FALSE,estimated=FALSE)

if(combine==FALSE){
W0<-kdmMotifs(pos_CT_init,initialization="cluster",rank=mot,np=0,init.only=TRUE)
W<-kdmMotifs(pos_CT,initialization="custom",W0=W0@.Data,np=0,rank=mot)
}

saveRDS(W,paste0("Models/W_",project,"_",length,"_",gaps,"_",mot,".rds"))

pos_features<-cbind(kdmFeatures(W,pos_counts$left),kdmFeatures(W,pos_counts$centre),kdmFeatures(W,pos_counts$right))
neg_features<-cbind(kdmFeatures(W,neg_counts$left),kdmFeatures(W,neg_counts$centre),kdmFeatures(W,neg_counts$right))

coefficients_all<-data.frame()
result<-data.frame()
for (i in 1:10){
train_features<-rbind(pos_features[pos_ids!=i,],neg_features[neg_ids!=i,])
test_features<-rbind(pos_features[pos_ids==i,],neg_features[neg_ids==i,])
train_label<-c(rep(1,sum(pos_ids!=i)),rep(0,sum(neg_ids!=i)))
test_label<-c(rep(1,sum(pos_ids==i)),rep(0,sum(neg_ids==i)))


mm<-cv.glmnet(y=train_label,x=train_features, family="binomial", type.measure="auc",alpha=alpha,foldid=lista_ids[[i]])
result<-rbind(result,data.frame(i,test_label,predict(mm,test_features,type="link",s="lambda.min")))
coefficients_all<-rbind(coefficients_all,data.frame(model=c("intercept",rep(seq(1,mot),3)),window=c("intercept",rep("left",mot),rep("centre",mot),rep("right",mot)),cv=i,coeff=coef(mm,s="lambda.min")[,1]))
}

write.table(coefficients_all,paste0("Models/",project,"_",length,"_",gaps,"_",mot,".txt"),col.names=TRUE,row.names=FALSE,quote=FALSE)

data<-read.table("GraphProt_CLIP_performances/CV_new_performances.txt",header=TRUE)
pf<-performance(prediction(labels=result$test_label,predictions=result$lambda.min),"auc")@y.values[[1]]
colnames(result)[1]<-"index"
perfs<-c()
for(i in 1:10){
tmp<-subset(result,index==i)
perfs<-c(perfs,performance(prediction(labels=tmp$test_label,predictions=tmp$lambda.min),"auc")@y.values[[1]])
}
pf_mean<-mean(perfs)
pf_sd<-sd(perfs)

#message(c(pf,pf_mean,pf_sd))
data<-rbind(data,c(project,length,gaps,mot,NA,init,alpha,combine,pf,pf_mean,pf_sd))

colnames(data)<-c("PROJECT","LENGTH","GAPS","MOTIFS","MOTIFS_USED","INITIALIZATION","ALPHA","POSITIVE_and_NEGATIVE","PERFORMANCE","mean_PERFORMANCE","sd_PERFORMANCE")
write.table(data,"GraphProt_CLIP_performances/CV_new_performances.txt",col.names=TRUE,row.names=FALSE,quote=FALSE)

end<-Sys.time()
time<-end-start
print(time)


