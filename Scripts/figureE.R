library(dplyr)
library(reshape2)
library(ggplot2)

load("/adat/Progetti/KDM/MATERIALE/kdmCentrimo/jaccard_RBP_5cases_nopos.Rdata")

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

mm=get_info(jaccard[[1]])
ggplot(mm, aes(x = variable, y = prop, fill = Class)) +
    geom_bar(stat = "identity", col = "black", position = "stack",width=0.8) +
    geom_text(aes(label = paste0(label," ", test), y = ypos), size = 3, color = "black",fontface="bold") +
    geom_text(data=subset(mm,cum==tot),aes(label=tot,y=cprop),vjust=-1)+
    theme_minimal() +
    scale_fill_manual(values =c("grey90",brewer.pal(3,"Blues")[c(1,3)]))+
    theme(panel.grid = element_blank(),axis.line = element_line(linewidth = 0.5),plot.title = element_text(size = 10))+
    xlab("Method")+ylab("Percentage")
    
    
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

ggarrange(p0,ggarrange(plotlist=plots,nrow=1,ncol=5,common.legend=T),nrow=2,heights=c(2,8))