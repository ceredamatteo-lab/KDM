library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(irr)

if(Sys.info()["user"]=="tbecchi"){load("/Users/tbecchi/Desktop/CLIP_maps/jaccard_RBP_5cases_nopos_binom.Rdata")
  } else {load("/adat/Progetti/KDM/MATERIALE/kdmCentrimo/jaccard_RBP_5cases_nopos_binom.Rdata")}

# BARPLOT ----


get_info=function(data){
    data=data%>%mutate(ID=paste0(target,"_",cell))
    data2=melt(data%>%select(ID,rank_centrimo,rank_kdm),id.vard="ID")%>%
        mutate(Class=ifelse(value<=5,"<=5",ifelse(value!=Inf,">5","Not.found")))%>%
        count(variable, Class)%>%group_by(variable)%>%mutate(tot=sum(n),prop = n / tot)%>%
        arrange(variable, Class)%>%
        mutate(cum=cumsum(n),cprop=cumsum(prop),label = paste0(n,"(",round(prop * 100), "%)"),ypos = cumsum(prop) - prop / 2)%>%mutate(Class=factor(Class,levels=c("Not.found",">5","<=5")))
    test <- binom.test(sum(data$rank_kdm <= 5), nrow(data), sum(data$rank_centrimo <= 5) / nrow(data))
    data2 <- data2 %>%
    mutate(test = if_else(variable == "rank_kdm" & Class == "<=5",
        case_when(test$p.value > 0.05 ~ "=",
                    sum(data$rank_centrimo <= 5) > sum(data$rank_kdm <= 5) ~ "-",
                    TRUE ~ "+"
                    ),""))

    return(data2)
}

anno=readRDS("/adat/Progetti/KDM/MATERIALE/kdmCentrimo/annotation_ranking_methods.rds")
p0=ggplot(anno,aes(x=Method,y=Measure))+geom_tile(fill="white",col="black")+geom_text(aes(label=Value))+facet_grid(~tmp)+theme_minimal()+scale_y_discrete(limits=rev)+theme(panel.grid=element_blank(),axis.title=element_blank(),strip.text=element_blank())

plots=list()
for(i in c(3,2,1,5,4)){
    mm=get_info(jaccard[[i]])
    plots[[length(plots)+1]]=ggplot(mm, aes(x = variable, y = prop, fill = Class)) +
                            geom_bar(stat = "identity", col = "black", position = "stack",width=0.8) +
                            geom_text(aes(label = paste0(label," ", test), y = ypos), size = 3, color = "black",fontface="bold") +
                            geom_text(data=subset(mm,cum==tot),aes(label=tot,y=cprop),vjust=-1)+
                            theme_minimal() +
                            scale_fill_manual(values =c("grey90",brewer.pal(3,"Blues")[c(1,3)]))+
                            theme(panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10))+
                            xlab("Method")+ylab("Percentage")
}


# JACCARD ----


breaks=c(0,5,10000,Inf)
labels=c("1-5",">5","Not selected")
plots2=list()
for(i in c(3,2,1,5,4)){
  d=jaccard[[i]]%>%
    mutate(class.centrimo=cut(rank_centrimo,breaks = breaks,labels=labels),
           class.kdm=cut(rank_kdm,breaks = breaks,labels=labels))
  all_combos <- expand.grid(class.centrimo = unique(d$class.centrimo),class.kdm = unique(d$class.kdm))
  tmp <- d %>%
    count(class.centrimo,class.kdm) %>%
    right_join(all_combos, by = c("class.centrimo","class.kdm")) %>%
    mutate(Freq = ifelse(is.na(n), 0, n)) %>%
    group_by(class.centrimo) %>% mutate(tot.Var1 = sum(Freq)) %>% ungroup() %>%
    group_by(class.kdm) %>% mutate(tot.Var2 = sum(Freq)) %>% ungroup() %>%
    mutate(Jaccard = Freq / (tot.Var1 + tot.Var2 - Freq))
  kk=kappa2(data.frame(d$class.centrimo, d$class.kdm))$value
  plots2[[length(plots2)+1]]=ggplot(tmp,aes(x=class.kdm,y=class.centrimo,fill=Jaccard))+
    scale_y_discrete(limits=rev)+
    geom_point(col="black",shape=21,size=9)+
    coord_equal()+
    geom_text(aes(label=Freq))+scale_fill_gradient(low = "white",high = "purple")+theme_bw()+
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,hjust=1),plot.title = element_text(size=8))+
    ylab("Centrimo")+xlab("KDM.Centrimo")+ggtitle(paste0("K.Cohen coefficient: ",round(kk,3)))
}

# BARPLOT + JACCARD ----

ggarrange(p0,ggarrange(plotlist=plots,nrow=1,ncol=5,common.legend=T,legend = "bottom"),
          ggarrange(plotlist=plots2,nrow=1,ncol=5,common.legend=T,legend = "bottom"),nrow=3,heights=c(1,4,4))


# JACCARD P.VALUES ----

make_vec <- function(L, x, I, U, seed = 30580) {
  set.seed(seed)
  A = c(rep(1, x), rep(0, L - x))
  
  nB = I + (U - x)
  B = rep(0, L)
  inter_idx = sample(1:x, I)
  B[inter_idx] = 1
  remaining_1 = nB - I
  if (remaining_1 > 0) {
    zero_idx = which(A == 0)
    add_idx = sample(zero_idx, remaining_1)
    B[add_idx] = 1
  }
  
  return(list(A,B))
}

compute_pval_row = function(j, nA, intAB, unionAB, union_size, n_perm) {
  if (is.na(j)) return(c(NA,NA))
  if (j == 0) return(c(1,1))
  V = make_vec(union_size, nA, intAB, unionAB)
  test = jaccard.test(V[[1]], V[[2]], method = "bootstrap", B = n_perm, seed = 30580)
  return(c(test$pvalue, test$expectation))
}

n_perm=10000
union_size=1071
library(jaccard)

new_jaccard=list()
for(k in 1:length(jaccard)){
  print(k)
  tmp=jaccard[[k]]
  tmp$J_pval=NA
  tmp$Exp=NA
  
  for(i in 1:nrow(tmp)){
    tt=compute_pval_row(j = tmp$j[i],nA = tmp$ncentrimo[i],intAB = tmp$int[i],unionAB = tmp$union[i],union_size = union_size,n_perm = n_perm)
    tmp$J_pval[i]=tt[1]
    tmp$Exp[i]=tt[2]
  }
  
  tmp$J_pval_adj = p.adjust(tmp$J_pval, method = "bonferroni")
  new_jaccard[[k]]=tmp
}
#saveRDS(new_jaccard,"jaccard_RBP_5cases_nopos_binom_pval.Rdata")

labels = paste0("M", seq_along(new_jaccard))

jaccard_tot = bind_rows(
  lapply(seq_along(new_jaccard), function(i) {
    cbind(new_jaccard[[i]], Model = labels[i])
  })
)

jaccard_tot <- jaccard_tot %>%group_by(Model) %>%mutate(J_pval_adj_new = p.adjust(J_pval, method = "fdr")) %>%ungroup()

ggplot(jaccard_tot%>%filter(!is.na(j),Model%in%c("M1","M3","M5"))%>%mutate(Model=factor(Model,levels=c("M3","M1","M5"))),aes(x=Model,y=j,fill=Model))+geom_boxplot(notch = T)+
  theme_bw()+stat_compare_means(comparisons = list(c("M1","M3"),c("M1","M5"),c("M3","M5")))
table(jaccard[[4]]$j==jaccard[[5]]$j)

ggplot(jaccard_tot%>%filter(!is.na(j),Model%in%c("M1","M3","M5"))%>%mutate(Model=factor(Model,levels=c("M3","M1","M5")),Is.Significant=J_pval_adj<=0.1&j>=Exp))+
  geom_bar(position = "fill",aes(x=Model,alpha=Is.Significant,fill=Model),col="black")+theme_bw()+scale_alpha_manual(values = c(0.3,1))+
  theme(panel.grid = element_blank())

df <- jaccard_tot %>%
  filter(!is.na(j), Model %in% c("M1","M3","M5")) %>%
  mutate(Model = factor(Model, levels = c("M3","M1","M5")),
         Is.Significant = J_pval_adj_new <= 0.1 & j >= Exp
  ) %>%
  group_by(Model, Is.Significant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Model) %>%
  mutate(prop = n / sum(n),ypos = prop,label = ifelse(Is.Significant, paste0(round(100*prop), "%"), ""))

ggplot(df, aes(x = Model, y = prop, fill = Model, alpha = Is.Significant)) +
  geom_bar(stat = "identity", col = "black", position = "fill") +
  geom_text(aes(y = ypos, label = label), color = "black",vjust=1.5,show.legend = F) +
  scale_alpha_manual(values = c(0.3, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())
