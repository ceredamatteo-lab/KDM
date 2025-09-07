
top_centrimo=readRDS("profiles/centrimo_top5.rds")


result=data.frame(Experiment=exp_info$Experiment,target=exp_info$target,
  sel_centrimo=0,sel_kdm=0,centrimo.target=FALSE,kdm.target=FALSE,Intersect=0)

for(u in 1:nrow(exp_info)){
  print(u)
  exp=exp_info$Experiment[u]
  target=exp_info$target[u]

  pn=readRDS(paste0("profiles/",exp,".rds"))
  ne<-kdmCentrimo(pn,symmetric=FALSE,tail="upper")
  res<-c()
  for(k in 1:ncol(mcross_W)){
    E<-ne[[k]]
    cpos<-pn$x[2:(length(pn$x)-1)][floor((E$Start+E$End)/2)]
    cond<-which(cpos<=0 & (E$Centrality_q.value * ncol(mcross_W)) < 0.05 & (E$Enrichment_q.value  * ncol(mcross_W))<0.05 & E$Centrality>1 & E$Enrichment>1)
    if(length(cond)>0){
      ss<-E[cond[order(rank(E$Centrality[cond])+rank(E$Enrichment[cond]),E$Centrality[cond],decreasing=c(TRUE,TRUE))],]
      res<-rbind(res,data.frame(motif=k,ss[1,]))
        
    }
  }

  selected_kdm=c()
  if(!is.null(res)){
    res<-res[order(rank(res$Centrality)+rank(res$Enrichment),res$Centrality,decreasing=c(TRUE,TRUE)),]
    selected_kdm<-colnames(mcross_W)[res$motif[1:min(5,nrow(res))]]
    result$sel_kdm[u]=length(selected_kdm)
    result$kdm.target[u]=target%in%do.call(rbind,strsplit(selected_kdm,"\\."))[,2]
  }

  t=which(names(top_centrimo)==exp)
  sel_centrimo=top_centrimo[[t]] 
  if(length(sel_centrimo)>0){
    result$sel_centrimo[u]=length(sel_centrimo)
    result$centrimo.target[u]=target%in%do.call(rbind,strsplit(sel_centrimo,"\\."))[,2]
  }

  if(length(selected_kdm)>0&length(sel_centrimo)>0){result$Intersect[u]=sum(selected_kdm%in%sel_centrimo)}
}
