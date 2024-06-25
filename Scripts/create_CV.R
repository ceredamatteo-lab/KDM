if( Sys.info()['nodename']%in%c('aulonocara.bp.lnf.it')){

  DATA_DIR="/home/tommaso.becchi/Documents/GraphProt/"
  #CGB_DIR=''
  #GIT_LOCAL_DIR=''
  #CGB_SHARED='~/Dropbox (HuGeF)/SPLICING/CGB_shared/'

  }

setwd(DATA_DIR)


library(stringr)
library(stringi)

return_150margins_seq<-function(input){
seq<-str_replace_all(str_replace_all(input[seq(2,length(input),2)],"u","t"),"U","T")
data<-data.frame(index=seq(1:length(seq)),marg_left=regexpr("[A-Z]", seq)-1,marg_right=regexpr("[A-Z]", stri_reverse(seq))-1)
seq<-seq[subset(data,marg_right==150 & marg_left ==150)$index]
return(seq)
}


projects<-read.table("GraphProt_experiments.txt")$V1

for (project in projects){

print(project)

pos<-return_150margins_seq(scan(paste0("GraphProt_CLIP_sequences/",project,".train.positives.fa"),what="char"))
neg<-return_150margins_seq(scan(paste0("GraphProt_CLIP_sequences/",project,".train.negatives.fa"),what="char"))


pos_ids<-sample(seq(1:10),length(pos),replace=TRUE,prob=rep(1/10,10))
neg_ids<-sample(seq(1:10),length(neg),replace=TRUE,prob=rep(1/10,10))

lista_ids<-list()
for(i in 1:10){
pos_res<-length(pos)-sum(pos_ids==i)
neg_res<-length(neg)-sum(neg_ids==i)
res<-pos_res+neg_res
cv_ids<-sample(seq(1:10),res,replace=TRUE)
lista_ids[[i]]<-cv_ids
}

final<-list()
final$pos_seq<-pos
final$neg_seq<-neg
final$pos_ids<-pos_ids
final$neg_ids<-neg_ids
final$cv<-lista_ids

saveRDS(final,paste0("Rdata/",project,"_cv.rds"))
}
