library(ggplot2)
library(tidyr)
library(data.table)
library(gridExtra)
library(stringr)
library(ggplotify)


# upload the curated version of the drug table and clinical trial table

drugs <- fread("data/kras.drugs.csv",na.strings = "")
kras.trials <- fread("data/kras.clinical.trials.csv",na.strings = "")

# a table 

x <- drugs[CTrials!="",list("NCI\ncode"=code,Name=short.Name,Type,Target,CTrials)]
x <- x[order(Target,Type)]
x[,CTrials:=str_wrap(CTrials %>% gsub("\\|",", ",.),width = 13*4)]
names(x)[5] <- "CRC Clinical Trials"

x[grep("^G",Target),Target:=paste("KRAS",Target)]
median(1:nrow(x))

ggplot(x) + geom_bar(aes(Target,fill=Type),color="grey25",linewidth=.5) +
  scale_fill_manual(values=c(Inhibitor="#99d6ea",
                             Degrader="#90a8c3",
                             Vaccine="#97a97c",
                             TCR="#cfe1b9",
                             Antisense="#e5f993")) +
  ylab("Number of Compounds") +
  xlab("\nTarget K-Ras Mutation or Protein") +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme(text = element_text(size=16),
    panel.background = element_rect(fill="grey98",color="black",linewidth = 1),
        panel.grid = element_line(color="grey75"),
        legend.key = element_rect(color=0)) +
  ggtitle("K-Ras targeting compounds in clinical trials with CRC patients") -> g1
  
tableA <- as.ggplot(tableGrob(x[1:42],rows=NULL))
tableB <- as.ggplot(tableGrob(x[43:85],rows=NULL))

pdf("output/table.drugs.pdf",height=18,width=20)
grid.arrange(tableA,tableB,g1,layout_matrix=matrix(c(1,2,3,NA),ncol=2,byrow=T),heights=c(10,3))
dev.off()

tiff("output/table.drugs.tiff",height=18,width=20,res = 150,units = "in")
grid.arrange(tableA,tableB,g1,layout_matrix=matrix(c(1,2,3,NA),ncol=2,byrow=T),heights=c(10,3))
dev.off()

png("output/table.drugs.png",height=18,width=20, res = 150, units = "in")
grid.arrange(tableA,tableB,g1,layout_matrix=matrix(c(1,2,3,NA),ncol=2,byrow=T),heights=c(10,3))
dev.off()

# the table figure has very small fontype
# exploring other alternatives

library(networkD3)
library(ggalluvial)

library(dplyr)
library(forcats)



links1 <- drugs[,.(value=.N),by=list(source=ifelse(CTrials=="","Not in Clinical Trials","In Clinical Trials"),target=Type)]
links2 <- drugs[,.(value=.N),by=list(source=Type,target=Target)]
#links3 <- drugs[,.(value=.N),by=list(source=ifelse(CTrials=="","Not in CTs","In CT"),target=Target)]
links <- rbind(links1,links2)

nodes <- data.frame(name=unique(c(links$source,links$target)))

links[,.(Total=sum(value)),by=source]

links2 <- links

links2$source <- match(links$source,nodes$name) - 1
links2$target <- match(links$target,nodes$name) - 1

sankeyNetwork(Links = links2, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize = 25, nodeWidth = 75,fontFamily = "Arial",
              height=600,width=1200)



ggplot(links) + aes(axis1 = "source", axis)


df <- drugs[,list(Clinical_Trials=ifelse(CTrials=="","NO","YES"),Type,Target,value=1)]

df <- df %>%
  mutate(
    Clinical_Trials = fct_infreq(Clinical_Trials),
    Type = fct_infreq(Type),
    Target = factor(Target,c("G12C","G12D","G12V","G13N","multi-KRAS","pan-KRAS","pan-RAS","SOS1","SHP2")))


df$Type
table(df$Type)

alluvium.colors <- c(YES="gold2",NO="grey")

foo <- colorRampPalette(c("blue3","lightblue","white"))(5)
names(foo) <- levels(df$Type)

alluvium.colors <- c(alluvium.colors,foo)

foo <- colorRampPalette(c("pink3","yellow3","green3"))(9)
names(foo) <- levels(df$Target)

alluvium.colors <- c(alluvium.colors,foo)
rm(foo)

ggplot(df,
       aes(axis1 = Type, axis2 = Target, axis3 = Clinical_Trials, y = value)) +
  geom_alluvium(aes(fill=Target),color="white",width = 1/6,alpha=.5,show.legend = F,lwd=.5) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6, color = "black",show.legend = F,alpha=.7) +
  geom_point(stat = "stratum",aes(x=after_stat(x),y=after_stat(y))) +
  geom_label_repel(stat = "stratum",
                   point.padding = 0,
                   min.segment.length = 0,
                   aes(label = sprintf("%s  (n=%i)",after_stat(stratum),after_stat(count))), 
                   size = 4,
                   nudge_x = rep(c(-5,0,1.2)*.2,c(5,9,2)),
                   box.padding = .4,
                   fill="#FFFFFFAA") +
  geom_text(aes(x=1:3,y=110,label=c("Drug\nType","\nTarget","CRC\nClinical Trials")),
            data=data.frame(),inherit.aes = F,size=5,vjust=0) +
  theme_void() +
  scale_fill_manual("In Clinical Trials",values=alluvium.colors) +
  coord_cartesian(xlim = c(.61,3.4),ylim=c(-1,120)) -> alluvial.plot

ggsave("output/alluvial.png",alluvial.plot,device="png",width=12,height=5,dpi=450,bg = "white")





# plots of the Clinical Trials -----

# adjust dates for ggplot representation

kras.trials[nchar(`Start Date`) > 8,start:=as.Date(`Start Date`)]
kras.trials[nchar(`Start Date`) == 7,start:=as.Date(paste0(`Start Date`,"-01"))]

kras.trials[nchar(`Primary Completion Date`) > 8,completion:=as.Date(`Primary Completion Date`)]
kras.trials[nchar(`Primary Completion Date`) == 7,completion:=as.Date(paste0(`Primary Completion Date`,"-01"))]
kras.trials[is.na(completion)]

my.palette <- list()
my.palette$study.status <- c(TERMINATED="#bb5256",
                             COMPLETED="#E73D23",
                             UNKNOWN="#CCC",
                             ACTIVE_NOT_RECRUITING="#438BC7",
                             RECRUITING="#3AA691",
                             ENROLLING_BY_INVITATION="#59C5AF",
                             NOT_YET_RECRUITING="#95DACC")


x <- kras.trials
x[completion > `Last Update Posted`,last:=completion]
x[completion <= `Last Update Posted`,last:=`Last Update Posted`]
x[,target := Target %>% split1 %>% unlist %>% sort %>% unique %>% join1,by=`NCT Number`]
x <- x[order(Target,start)]
x[,`NCT Number`:=factor(`NCT Number`,levels=as.character(`NCT Number`))]

x$`Study Status` <- factor(x$`Study Status`,names(my.palette$study.status))

trial.plot <- list(
  aes(start,`NCT Number`),
  geom_vline(xintercept = as.Date("2025-03-15"),lwd=3,color="lightgrey"),
  geom_segment(aes(yend=`NCT Number`,xend=`Last Update Posted`),lwd=0.5,lty=5),
  geom_segment(aes(yend=`NCT Number`,xend=completion,color=`Study Status`),lwd=2,show.legend=T),
  xlab(""),
  ylab(""),
  scale_x_date(limits=as.Date(c("2012-01-01","2045-12-31"))),
  scale_y_discrete(limits=rev),
  geom_point(aes(x=as.Date(`Last Update Posted`)),pch=21,fill="white"),
  geom_text(aes(label=target,x=start),hjust=1,nudge_x=-60,size=3),
  scale_color_manual("Study Status",values=my.palette$study.status,drop=F),
  theme(strip.text.x.top = element_text(size=16,face = "bold"),
        panel.background = element_rect(fill="white",color="black"),
        panel.grid = element_line(color="#EEE"),
        strip.background = element_rect(color="black"),
        legend.key.height = unit(10, "pt"))
  
)

active <- ! x$`Study Status` %in% c("COMPLETED","TERMINATED")

g0 <- ggplot(x[!active]) + trial.plot +
  facet_wrap(~"Completed or terminated") +   
  scale_x_date(limits=as.Date(c("2000-01-01","2035-12-31"))) +
  theme(legend.position = "none") +
  geom_text(aes(label=sprintf("%s",Drug),x=last),hjust=0,nudge_x = 240,size=3.2) 

g1 <- ggplot(x[active & Phases=="PHASE1"]) + trial.plot + facet_wrap(~ "Phase I") + 
  theme(legend.position = "none") +
  geom_text(aes(label=sprintf("%s",Drug),x=last),hjust=0,nudge_x = 60,size=3.2) 


g2 <- ggplot(x[active & grepl("PHASE2",Phases)]) + trial.plot + facet_wrap(~ "Phase I/II or II") + 
  theme(legend.position = "none") +
  geom_text(aes(label=sprintf("%s",Drug),x=last),hjust=0,nudge_x = 60,size=3.2) 


g3 <- ggplot(x[active & Phases=="PHASE3"]) + trial.plot + facet_wrap(~ "Phase III") + 
  geom_text(aes(label=sprintf("%s (%s)",Drug,Acronym),x=last),hjust=0,nudge_x = 60)  

g4 <- grid.arrange(g0,g1,g2,g3,layout_matrix=matrix(c(2,3,1,2,4,4),byrow=T,ncol=3),heights=c(9,2),widths=c(3,3,2.5))

g3 <- ggplot(x[active & Phases=="PHASE3"]) + trial.plot + facet_wrap(~ "Phase III") + 
  geom_text(aes(label=sprintf("%s (%s)",Drug,Acronym),x=last),hjust=0,nudge_x = 60) +
  scale_color_manual("Study Status",values=my.palette$study.status[3:7],drop=F) +
  theme(legend.position = "bottom",legend.key = element_rect(color=0)) +
  guides(colour = guide_legend(direction = "vertical"))

g3b <- ggplot(x[active & Phases=="PHASE3"]) + trial.plot + facet_wrap(~ "Phase III") + 
  geom_text(aes(label=sprintf("%s (%s)",Drug,Acronym),x=last),hjust=0,nudge_x = 60) +
  theme(legend.position = "none") 
  

g5 <- ggplot(x2) + aes(Phases,N) + geom_col(aes(fill=`Study Status`)) + 
  coord_flip() +
  scale_x_discrete(limits=rev) +
  geom_text(aes(label=N,group=`Study Status`),position=position_stack(vjust=0.5),size=3) +
  scale_fill_manual(values=my.palette$study.status) +
  ylab("\nNumber of clinical trials") +
  xlab("") +
  theme(panel.background = element_rect(fill="white",color="black"),
        panel.grid = element_line(color="#EFEFEF"),
        legend.key = element_rect(color=0),
        legend.key.size = unit(.4,'cm'),
        legend.text = element_text(size=10)) 


g6 <- grid.arrange(g1,g2,g3b,g0,g5,layout_matrix=matrix(c(1,2,1,3,1,4,5,5),byrow=T,ncol=2),heights=c(5,1,3.5,1.5))



ggsave(filename = "output/ClinicalTrials.tiff",plot = g4,width = 18,height = 10,dpi=200)
ggsave(filename = "output/ClinicalTrials.png",plot = g4,width = 18,height = 10,dpi=200)
ggsave(filename = "output/ClinicalTrials.pdf",plot = g4,width = 18,height = 10)

ggsave(filename = "output/ClinicalTrials2.tiff",plot = g6,width = 12,height = 16,dpi=200)
ggsave(filename = "output/ClinicalTrials2.png",plot = g6,width = 12,height = 16,dpi=200)
ggsave(filename = "output/ClinicalTrials2.pdf",plot = g6,width = 12,height = 16)



