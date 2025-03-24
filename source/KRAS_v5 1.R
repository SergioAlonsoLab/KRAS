# load required libraries ----

library(data.table)
library(tidyr)
library(forcats)
library(ggrepel)
library(ggplot2)
library(ggtext)
library(ggnewscale)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(corrplot)
library(DescTools)
library(emmeans)

# Set working directory -----

try(setwd("/imppc/labs/cge/salonso/KRAS/"),silent=T)
try(setwd("~/Documents/Work/KRAS/"),silent=T)

# color palettes and ggplot theme ----

palette <- list()

palette$msi <- c(MSS="#BCBDDC",MSI="#6A51A3","MSI/HYPER"="#6A51A3")
palette$mut <- c(WT="#FEE6CE",MUT="#F16913")
palette$study <- c("#C19F08","#369E4B","#357EB9","#FA792F","#ce3206")
names(palette$study) <- c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis")
palette$study2 <- brewer.pal(9,"Greens")[c(3,5:9)]
names(palette$study2) <- c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis")

palette$stage <- brewer.pal(9,"Blues")[c(2,4,6,8)]
names(palette$stage) <- c("I","II","III","IV")
palette$sample.type <- c("Primary"="#8EBEDA","Metastasis"="#2478B1")
palette$SigLevel <- c("<0.0001"="#c32f27",
                             "<0.001"="#d8572a",
                             "<0.01"="#db7c26",
                             "<0.05"="#f7b538",
                             "NS"="grey75")

mytheme <- theme(panel.background = element_rect(fill="grey90"),
                 strip.background = element_rect(fill="grey75",colour="grey30"),
                 panel.grid = element_line(color="white"),
                 panel.border = element_rect(color="grey30",fill=NA,linewidth=rel(1)),
) 


# panel labeller

panelTag <- function(tag="A",x=1,y=0,fsize=20) {
  textGrob(label=tag,x=unit(x,"cm"),y=unit(y,"cm"),gp = gpar(fontsize=fsize))
}



# load the bowel data ----

bowel <- readRDS("data/bowel.rds")
genes <- fread("data/all.genes.length.csv")

# KRAS, NRAS and HRAS mutations in the different datasets ----

# Figure 1: Ras vs primary / metastasis ------

selected.genes <- intersect(genes$hgnc_symbol,colnames(bowel))

x <- melt(bowel,id.vars = c("SEX","GROUP","SAMPLE_TYPE","STUDY","STAGE2"),measure.vars = selected.genes,variable.name = "GENE",value.name = "STATUS")
x[,SAMPLE_TYPE:=factor(SAMPLE_TYPE,c("Primary","Metastasis"))]

goi <- c("APC","TP53","FBXW7","SMAD4","KRAS","NRAS","HRAS","PIK3CA","BRAF")
attr(goi,"definition") <- "Genes of"

gg1a <- x[GENE %in% goi & !is.na(STAGE2),list(MUT=sum(STATUS!="WT"),.N),by=list(GROUP,STAGE2,GENE)][,GENE:=factor(GENE,goi)][,FREQ:=MUT/N*100] %>% 
  ggplot() + aes(FREQ,GENE) + 
  geom_col(aes(fill=STAGE2),position=position_dodge2(reverse=TRUE)) +
  geom_text(aes(group=STAGE2,label=sprintf(ifelse(FREQ==0 | FREQ > 1," %1.0f%% (%i/%i)"," %1.1f%% (%i/%i)"),FREQ,MUT,N)), position=position_dodge2(reverse=TRUE,width=.9),hjust=0,size=3) +
  scale_y_discrete(limits=rev) +
  facet_grid(cols=vars(GROUP)) +
  scale_x_continuous(limits=c(0,100)) +
  ylab("") + xlab("Percentage of mutated cases") +
  myLightTheme +
  theme(text=element_text(size=12),axis.text.y=element_text(face="italic",size=13),strip.text = element_text(size=13)) +
  scale_fill_manual("Stage",values=c("#91c1d8","#0070ac","#003371")) 

lapply(goi, function(gene) {
  
  m1 <- glm(I(STATUS!="WT") ~ SEX + GROUP + STAGE2,x[GENE==gene],family=binomial)
  m1 <- list(emmeans(m1,~SEX),
             emmeans(m1,~GROUP),
             emmeans(m1,~STAGE2))
  m1 <- lapply(m1,pairs,reverse=F)
  m1 <- cbind(lapply(m1,summary,adjust="fdr") %>% rbindlist,
              lapply(m1,function(x) confint(x)[,5:6]) %>% rbindlist,
              Gene=gene,
              Model="All Samples")
  
  m2 <- glm(I(STATUS!="WT") ~ SEX + GROUP + STAGE2 + STUDY,x[GENE==gene & STAGE2!="M"],family=binomial)
  m2 <- list(emmeans(m2,~SEX),
             emmeans(m2,~GROUP),
             emmeans(m2,~STAGE2),
             emmeans(m2,~STUDY))
  m2 <- lapply(m2,pairs,reverse=F)
  m2 <- cbind(lapply(m2,summary,adjust="fdr") %>% rbindlist,
              lapply(m2,function(x) confint(x)[,5:6]) %>% rbindlist,
              Gene=gene,
              Model="Primary Tumors")
  
  rbind(m1,m2)
  
}) %>% rbindlist -> stats1


stats1[,Gene:=factor(Gene,goi)]
stats1$contrast %>% gsub(" Primary","",.) %>% gsub("[\\(\\)]","",.) %>% gsub("I-","I/",.) %>% gsub("-","↔",.) -> stats1$contrast

stats1[,contrast:=factor(contrast,unique(contrast))]
stats1[,SigLevel:=cut(p.value,c(0,1e-4,1e-3,1e-2,5e-2,1))]
levels(stats1$SigLevel) <-c ("<0.0001","<0.001","<0.01","<0.05","NS")
stats1[,SigLevel := fct_rev(SigLevel)]

# for graphical purposes, it's preferable to switch the sign of the estimate and CIs

stats1[,estimate:=-estimate]
stats1[,foo:=-asymp.LCL]
stats1[,asymp.LCL:=-asymp.UCL]
stats1[,asymp.UCL:=foo]
stats1[,foo:=NULL]


gg1bc <- list(aes(estimate,contrast),
  geom_vline(xintercept = 0,lty=2),
  geom_segment(aes(x=asymp.LCL,xend=asymp.UCL,yend=contrast,color=SigLevel),lwd=3,show.legend=T),
  geom_point(),
  scale_y_discrete(limits=rev),
  ylab(""),
  xlab("log Odds Ratio"),
  facet_grid(rows=vars(Model),cols=vars(Gene),scales = "free"),
  scale_color_manual("Sig Level",values=palette$SigLevel,
                     drop=FALSE),
  myLightTheme,
  theme(strip.text.x.top = element_text(face="italic",size=13),
        axis.text.y = element_text(size=12)))


gg1b <- stats1[Gene %in% c("KRAS","NRAS","HRAS")] %>% ggplot() + gg1bc


# figure 1  
grid.arrange(arrangeGrob(gg1a,top=panelTag("A")),
             arrangeGrob(gg1b,top=panelTag("B")),
             heights=c(2,1.1))

 
# supplementary figure 1

grid.arrange(stats1[Gene %in% c("APC","TP53","FBXW7")] %>% ggplot() + gg1bc,
stats1[Gene %in% c("SMAD4","PIK3CA","BRAF")] %>% ggplot() + gg1bc)


# some stats shown in the body of the paper

bowel[,table(KRAS!="WT",SAMPLE_TYPE)]
stats1[Gene=="KRAS",sprintf("\n%s OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.2g",contrast,exp(estimate),exp(asymp.LCL),exp(asymp.UCL),p.value)] %>% cat()
stats1[Gene=="KRAS",sprintf("\n%s OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.2g",contrast,exp(-estimate),exp(-asymp.UCL),exp(-asymp.LCL),p.value)] %>% cat()
stats1[Gene=="HRAS",sprintf("\n%s OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.2g",contrast,exp(estimate),exp(asymp.LCL),exp(asymp.UCL),p.value)] %>% cat()

# identify paired primary - metastases samples -----
bowel[,SAMPLE_TYPE:=factor(SAMPLE_TYPE,c("Primary","Metastasis"))]
x <- bowel[order(patientId,SAMPLE_TYPE)][SAMPLE_COUNT > 1,list("Sample ID"=paste(sampleId,collapse="\n"),
                                                     "Sample Type"=paste(SAMPLE_TYPE,collapse="\n"),
                                                     "MSI status"=paste(GROUP,collapse="\n"),
                                                     KRAS=paste(KRAS,collapse="\n"),
                                                     NRAS=paste(NRAS,collapse="\n"),
                                                     HRAS=paste(HRAS,collapse="\n"),
                                                     BRAF=paste(BRAF,collapse="\n")),by=patientId] 
names(x)[1] <- "Patient ID"

dcast(bowel[SAMPLE_COUNT>1], patientId ~ SAMPLE_TYPE,fun.aggregate = length)[Primary>0 & Metastasis>0] 

png("output/table1.png",750,1300)
grid.table(x,rows=NULL)
dev.off()


# Co-mutation analysis using logistic regression ----

x <- bowel[SAMPLE_TYPE=="Primary",..selected.genes] != "WT"
x <- x[,colSums(x) >= 5]
x <- cbind(bowel[SAMPLE_TYPE=="Primary",list(GROUP)],x)

stats2 <- data.table()

for(y in c(goi,"POLE","KMT2D","RNF43")) {
  for(j in c(selected.genes)) {
    
    m1 <- x[,table(get(y),get(j),GROUP)]
    ft1 <- fisher.test(m1[,,1])[c("estimate","conf.int","p.value")] %>% unlist
    ft2 <- fisher.test(m1[,,2])[c("estimate","conf.int","p.value")] %>% unlist
    
    m1 <- data.table(rbind(c(ft1,m1[,,1]),c(ft2,m1[,,2])))
    m1[,Gene1:=y]
    m1[,Gene2:=j]
    m1[,GROUP:=c("MSS","MSI/HYPER")]
   
    stats2 <- rbind(stats2,m1)
    
  }
}

names(stats2)[5:8] <- c("WW","MW","WM","MM")
names(stats2)[1:3] <- c("OR","CI1","CI2")
stats2 <- stats2[Gene1!=Gene2]

stats2[,YulesQ:=(OR-1)/(OR+1)]
stats2[,YulesQlow:=(CI1-1)/(CI1+1)]
stats2[,YulesQhigh:=(CI2-1)/(CI2+1)]
stats2[is.infinite(OR),YulesQ:=1]
stats2[is.infinite(CI2),YulesQ:=1]

stats2[,Gene1:=factor(Gene1,c(goi,"POLE","KMT2D","RNF43"))]
stats2[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]


stats2[,p.adjusted:=p.adjust(p.value,"fdr"),by=list(GROUP,Gene1)]

stats2[,SigLevel:=cut(p.adjusted,c(0,1e-4,1e-3,1e-2,5e-2,1),c("<0.0001","<0.001","<0.01","<0.05","NS"))]
stats2[,SigLevel:=fct_rev(SigLevel)]

stats2[,MinM:=MM/(MM+MW)]
stats2[,MinW:=WM/(WM+WW)]

layersFig2 <-  list(aes(MinM*100,YulesQ) , 
  geom_hline(yintercept=0,lty=2) ,
  geom_segment(aes(xend=MinM*100,y=YulesQlow,yend=YulesQhigh,color=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),lwd=1,show.legend = T) ,
  geom_point(aes(color=SigLevel,size=SigLevel),show.legend = T) ,
  geom_text_repel(aes(label=Gene2),data=function(x) subset(x,p.adjusted < 0.05),
                  nudge_x = 10,min.segment.length = 0,fontface="italic",size=3) ,
  scale_color_manual("Sig Level",values=c(palette$SigLevel[1:4],NS="grey80"),
                     drop=FALSE) ,
  scale_size_manual("Sig Level",values=c(1,1.5,2,3,4),
                    drop=FALSE) ,
  ylab("Yule's Q") ,
  xlab("Co-mutation index (%)") ,
  facet_grid(cols=vars(Gene1),rows=vars(GROUP)) ,
  myLightTheme ,
  theme(strip.text.x.top = element_text(face="italic",size=13),
        strip.text.y.right = element_text(size=13),
        axis.title = element_text(size=13)))



fig2A <- ggplot(stats2[order(p.adjusted,decreasing = T)][Gene1 %in% c("KRAS","NRAS","HRAS")]) + layersFig2

fig2B <- stats2[p.adjusted < 0.05 & Gene1 %in% c("KRAS","NRAS","HRAS")][order(GROUP,Gene1,p.adjusted)][,list(GROUP,Gene1,Gene2,
                                                                                                    "OR [CI95%]"=sprintf("%1.2f [%1.2f-%1.2f]",OR,CI1,CI2),
                                                                                                    FDR=sprintf("%1.1e",p.adjusted),
                                                                                              "co-Mutated"=sprintf("%1.1f%%",MinM*100),
                                                                                                    "Mutated in WT"=sprintf("%1.1f%%",MinW*100))] 
fig2B <- tableGrob(fig2B,rows=NULL,
            theme=ttheme_default(base_size=10),
            widths=unit(c(2,2,2,3,2,2,2)*1.5,"cm"),
            heights=unit(rep(.55,nrow(fig2B)),"cm"))


grid.arrange(arrangeGrob(fig2A + theme(plot.margin=unit(c(.2,.6,.2,1),"cm")),top=panelTag("A")),
             arrangeGrob(fig2B,top=panelTag("B")),
             heights=c(3,2.2))


# for supplementary 

grid.arrange(
  ggplot(stats2[order(p.adjusted,decreasing = T)][Gene1 %in% c("APC","TP53","FBXW7")]) + layersFig2A,
  ggplot(stats2[order(p.adjusted,decreasing = T)][Gene1 %in% c("SMAD4","PIK3CA","BRAF")]) + layersFig2A
)

# co-mutational analysis by mutated codon ----

stats2B <- data.table()

selected.kras.codons <- c("G12","G13","Q61","K117","A146","Other")

for(i in selected.kras.codons) {
  for(j in c(selected.genes)) {
    tryCatch({
    m1 <- bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT",i),table(KRAS.CODON==i,get(j)!="WT",GROUP)]
    m2 <- bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON != "WT",table(KRAS.CODON==i,get(j)!="WT",GROUP)]
    ft1 <- fisher.test(m1[,,1])[c("estimate","conf.int","p.value")] %>% unlist
    ft2 <- fisher.test(m1[,,2])[c("estimate","conf.int","p.value")] %>% unlist
    
    ft3 <-  fisher.test(m2[,,1])[c("estimate","conf.int","p.value")] %>% unlist
    ft4 <-  fisher.test(m2[,,2])[c("estimate","conf.int","p.value")] %>% unlist
    
    
    m1 <- data.table(rbind(c(ft1,m1[,,1]),c(ft2,m1[,,2])))
    m1[,Mutation:=i]
    m1[,Gene2:=j]
    m1[,Comparison:="vsWT"]
    m1[,GROUP:=c("MSS","MSI/HYPER")]
    
    m2 <- data.table(rbind(c(ft3,m2[,,1]),c(ft4,m2[,,2])))
    m2[,Mutation:=i]
    m2[,Gene2:=j]
    m2[,Comparison:="vsMutants"]
    m2[,GROUP:=c("MSS","MSI/HYPER")]
    
    
    stats2B <- rbind(stats2B,rbind(m1,m2))},

    error=function(e) print(e))
  }
}

names(stats2B)[5:8] <- c("WW","MW","WM","MM")
names(stats2B)[1:3] <- c("OR","CI1","CI2")
stats2B <- stats2B[Gene2!="KRAS"]

stats2B[,YulesQ:=(OR-1)/(OR+1)]
stats2B[,YulesQlow:=(CI1-1)/(CI1+1)]
stats2B[,YulesQhigh:=(CI2-1)/(CI2+1)]
stats2B[is.infinite(OR),YulesQ:=1]
stats2B[is.infinite(CI2),YulesQ:=1]


stats2B[,p.adjusted:=p.adjust(p.value,"fdr"),by=list(GROUP,Mutation,Comparison)]

stats2B[,SigLevel:=cut(p.adjusted,c(0,1e-4,1e-3,1e-2,5e-2,1),c("<0.0001","<0.001","<0.01","<0.05","NS"))]
stats2B[,SigLevel:=fct_rev(SigLevel)]

stats2B[,MinM:=MM/(MM+MW)]
stats2B[,MinW:=WM/(WM+WW)]

stats2B[,Mutation:=factor(Mutation,selected.kras.codons,paste("KRAS",selected.kras.codons))]
stats2B[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]

ggplot(stats2B) + aes(MinM,YulesQ) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point(size=1,color="grey80") +
  geom_point(aes(color=SigLevel),data=stats2B[p.adjusted < 0.05],show.legend = T) +
  geom_text_repel(aes(label=Gene2),data=stats2B[p.adjusted < 0.05],fontface="italic",size=3,nudge_x = .1) +
  
  facet_wrap(~ GROUP + Mutation,ncol=3) +
  scale_color_manual("Sig Level",values=c(palette$SigLevel[1:4],NS="grey80"),
                     drop=FALSE) +
  scale_size_manual("Sig Level",values=c(1,1.5,2,3,4),
                    drop=FALSE) +
  ylab("Yule's Q") +
  xlab("Co-mutation index (%)") +
  myLightTheme

stats2B[,compLabel:=factor(Comparison,c("vsWT","vsMutants"),c("vs\nKRAS WT","vs\nother KRAS mutants"))]

fig6a <- ggplot(stats2B[Mutation %in% paste("KRAS",c("G12","G13"))]) + aes(MinM*100,YulesQ) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point(size=1,color="grey80") +
  geom_segment(aes(y=YulesQlow,yend=YulesQhigh,color=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),lwd=1,show.legend=T) +
  geom_point(aes(color=SigLevel,size=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),show.legend = T) +
  geom_text_repel(aes(label=Gene2),data=function(x) subset(x,p.adjusted < 0.05),fontface="italic",size=3,nudge_x = .1) +
  facet_grid(rows=vars(GROUP),cols=vars(paste(Mutation,compLabel,sep="\n"))) +
  scale_color_manual("Sig Level",values=c(palette$SigLevel[1:4],NS="grey80"),
                     drop=FALSE) +
  scale_size_manual("Sig Level",values=c(1,1.5,2,3,4),
                    drop=FALSE) +
  scale_x_continuous(limits=c(0,100)) +
  ylab("Yule's Q") +
  xlab("Co-mutation index (%)") +
  myLightTheme +
  theme(text=element_text(size=12),strip.text = element_text(size=13))


fig6b <- ggplot(stats2B[Mutation %in% paste("KRAS",c("Q61","K117"))]) + aes(MinM*100,YulesQ) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point(size=1,color="grey80") +
  geom_segment(aes(y=YulesQlow,yend=YulesQhigh,color=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),lwd=1,show.legend=T) +
  geom_point(aes(color=SigLevel,size=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),show.legend = T) +
  geom_text_repel(aes(label=Gene2),data=function(x) subset(x,p.adjusted < 0.05),fontface="italic",size=3,nudge_x = .1) +
  facet_grid(rows=vars(GROUP),cols=vars(paste(Mutation,compLabel,sep="\n"))) +
  scale_color_manual("Sig Level",values=c(palette$SigLevel[1:4],NS="grey80"),
                     drop=FALSE) +
  scale_size_manual("Sig Level",values=c(1,1.5,2,3,4),
                    drop=FALSE) +
  scale_x_continuous(limits=c(0,100)) +
  ylab("Yule's Q") +
  xlab("Co-mutation index (%)") +
  myLightTheme +
  theme(text=element_text(size=12),strip.text = element_text(size=13))


fig6c <- ggplot(stats2B[Mutation %in% paste("KRAS",c("A146"))]) + aes(MinM*100,YulesQ) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point(size=1,color="grey80") +
  geom_segment(aes(y=YulesQlow,yend=YulesQhigh,color=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),lwd=1,show.legend=T) +
  geom_point(aes(color=SigLevel,size=SigLevel),data=function(x) subset(x,p.adjusted < 0.05),show.legend = T) +
  geom_text_repel(aes(label=Gene2),data=function(x) subset(x,p.adjusted < 0.05),fontface="italic",size=3,nudge_x = .1) +
  facet_grid(rows=vars(GROUP),cols=vars(paste(Mutation,compLabel,sep="\n"))) +
  scale_color_manual("Sig Level",values=c(palette$SigLevel[1:4],NS="grey80"),
                     drop=FALSE) +
  scale_size_manual("Sig Level",values=c(1,1.5,2,3,4),
                    drop=FALSE) +
  scale_x_continuous(limits=c(0,100)) +
  ylab("Yule's Q") +
  xlab("Co-mutation index (%)") +
  myLightTheme +
  theme(text=element_text(size=12),strip.text = element_text(size=13))

fig6a

grid.arrange(fig6b,fig6c,layout_matrix=matrix(c(1,2,1,NA),ncol=2),widths=c(1.15,1))



bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","Q61","A146","K117") & GROUP=="MSS",glm(I(APC!="WT") ~ KRAS.CODON ,family=binomial)] -> glm1
summary(glm1)

kruskal2 <- function(sample.type="Primary",group,gene,codons=c("G12","G13","Q61","A146","K117")) {
  print(gene)
  bowel[SAMPLE_TYPE==sample.type & KRAS.CODON %in% codons & GROUP==group,kruskal.test(I(get(gene)!="WT") ~ KRAS.CODON)] 
}

glm2 <- function(sample.type="Primary",group,gene,codons=c("G12","G13","Q61","A146","K117")) {
  print(gene)
  bowel[SAMPLE_TYPE==sample.type & KRAS.CODON %in% codons & GROUP==group,glm(I(get(gene)!="WT") ~ KRAS.CODON,family=binomial)] -> glm1
  print(summary(glm1))
  print(drop1(glm1,test = "Chisq"))
  emmeans(glm1,~KRAS.CODON) %>% pairs() 
}

for(i in c("APC","PIK3CA","FBXW7","TP53","BRAF","NRAS")) {
  sprintf("\n\n\nNow analysing %s -----------------------\n",i) %>% cat
  kruskal2(gene=i,group="MSS") %>% print
  glm2(gene=i,group="MSS") %>% print
}

for(i in c("APC","PIK3CA","FBXW7","TP53","BRAF","NRAS")) {
  kruskal2(gene=i,group="MSI/HYPER") %>% print
  glm2(gene=i,group="MSI/HYPER") %>% print
}


matrix(c(374-301,301,22-13,13),ncol=2) %>% fisher.test()


melt(bowel,id.vars=c("SAMPLE_TYPE","GROUP","KRAS.CODON"),measure.vars = c("APC","PIK3CA","FBXW7","TP53","BRAF","NRAS"),variable.name = "Gene") -> foo
foo[,list(MUT=sum(value!="WT"),.N),by=list(SAMPLE_TYPE,GROUP,KRAS.CODON,Gene)] -> foo
foo[,label:=factor(sprintf("_%s_ mutations",Gene),sprintf("_%s_ mutations",levels(Gene)))]
foo[,KRAS.CODON:=fct_reorder(KRAS.CODON,-N)]




fig6layers <- list( aes(MUT/N*100,KRAS.CODON), 
                    geom_col(aes(fill=KRAS.CODON)),
                    scale_y_discrete(limits=rev,labels=function(x) sprintf("_KRAS_ %s",x)),
                    scale_x_continuous(limits=c(0,100)),
                    geom_text(aes(label=sprintf(ifelse(MUT/N>0.01 | MUT==0," %1.0f%% (%i/%i)"," %1.1f%% (%i/%i)"),MUT/N*100,MUT,N)),hjust=0,data=function(x) subset(x,MUT/N<.5)),
                    geom_text(aes(label=sprintf("%1.0f%% (%i/%i) ",MUT/N*100,MUT,N)),hjust=1,color="black",data=function(x) subset(x,MUT/N>=.5)),
                    myLightTheme,
                    theme(strip.text.x.top = element_markdown(size=13),
                          text=element_text(size=12),
                          axis.text.y = element_markdown(size=12),
                          legend.title = element_markdown()),
                    ylab(""),
                    xlab("Mutated cases (%)"),
                    scale_fill_manual("_KRAS_ Mutation",values=c("#AAAAFF",colorRampPalette(c("#EE6633","#FFAA33"))(5))),
                    facet_wrap(~ label))


fig6A <- foo[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","A146","Q61","K117") & GROUP=="MSS"] %>% ggplot() + fig6layers
fig6B <- foo[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","A146","Q61","K117") & GROUP=="MSI/HYPER"] %>% ggplot() + fig6layers

grid.arrange(arrangeGrob(fig6A + ggtitle("MSS Primary CRCs"),top=panelTag("A")),
             arrangeGrob(fig6B + ggtitle("MSI/hypermutated Primary CRCs"),top=panelTag("B")))



# Figure 7 (co-mutation with NF1 and RASA1)

melt(bowel,id.vars=c("SAMPLE_TYPE","GROUP","KRAS.CODON"),measure.vars = c("NF1","RASA1"),variable.name = "Gene") -> foo
foo[,list(MUT=sum(value!="WT"),.N),by=list(SAMPLE_TYPE,GROUP,KRAS.CODON,Gene)] -> foo
foo[,label:=factor(sprintf("_%s_ mutations",Gene),sprintf("_%s_ mutations",levels(Gene)))]
foo[,KRAS.CODON:=fct_reorder(KRAS.CODON,-N)]

fig7A <- foo[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","A146","Q61","K117") & GROUP=="MSS"] %>% ggplot() + fig6layers
fig7B <- foo[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","A146","Q61","K117") & GROUP=="MSI/HYPER"] %>% ggplot() + fig6layers

grid.arrange(arrangeGrob(fig7A + ggtitle("MSS Primary CRCs"),top=panelTag("A")),
             arrangeGrob(fig7B + ggtitle("MSI/hypermutated Primary CRCs"),top=panelTag("B")))

bowel[,table(GROUP,NF1!="WT",SAMPLE_TYPE)][,,1] %>% fisher.test()

kruskal2(sample.type = "Primary",gene = "NF1",group = "MSS",codons = c("WT","G12","G13","Q61","A146","K117"))
kruskal2(sample.type = "Primary",gene = "NF1",group = "MSS",codons = c("G12","G13","Q61","A146","K117"))
kruskal2(sample.type = "Primary",gene = "NF1",group = "MSI/HYPER",codons = c("WT","G12","G13","Q61","A146","K117"))
kruskal2(sample.type = "Primary",gene = "NF1",group = "MSI/HYPER",codons = c("G12","G13","Q61","A146","K117"))

kruskal2(sample.type = "Primary",gene = "RASA1",group = "MSS",codons = c("WT","G12","G13","Q61","A146","K117"))
kruskal2(sample.type = "Primary",gene = "RASA1",group = "MSS",codons = c("G12","G13","Q61","A146","K117"))



bowel[SAMPLE_TYPE=="Primary"& GROUP=="MSS" & KRAS.CODON %in% c("G12","G13","Q61","A146","K117"),glm(I(NF1!="WT") ~KRAS.CODON,family=binomial)] %>% emmeans(.,~KRAS.CODON) %>% pairs(adjust="fdr")
bowel[SAMPLE_TYPE=="Primary"& GROUP=="MSS" & KRAS.CODON %in% c("G12","G13","Q61","A146","K117"),glm(I(NF1!="WT") ~KRAS.CODON,family=binomial)] %>% summary

bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","Q61","K117","A146"),glm(I(NF1!="WT") ~ STUDY + AGE + STAGE2 + GROUP + KRAS.CODON )] -> glm1 
bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("G12","G13","Q61","K117","A146"),glm(I(NF1!="WT") ~ STUDY + AGE + STAGE2 + GROUP + KRAS.CODON )] -> glm2 

summary(glm1)
summary(glm2)

drop1(glm1,test = "Chisq")
drop1(glm2,test = "Chisq")



bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("WT","G12","G13","Q61","K117","A146"),glm(I(RASA1!="WT") ~ STUDY + AGE + STAGE2 + GROUP + KRAS.CODON )] -> glm1 
bowel[SAMPLE_TYPE=="Primary" & KRAS.CODON %in% c("G12","G13","Q61","K117","A146"),glm(I(RASA1!="WT") ~ STUDY + AGE + STAGE2 + GROUP + KRAS.CODON )] -> glm2 

summary(glm1)
summary(glm2)

drop1(glm1,test="Chisq")
drop1(glm2,test="Chisq")



# Mutation analysis by gene size ------

# a possible supplementary figure

geneLengths <- genes[transcript_is_canonical==1 & hgnc_symbol!="",list(cds_length=max(cds_length)),by=hgnc_symbol]
sandbox <- new.env()
load("data/raw.data.RData",sandbox,verbose=T)
cohorts <- sandbox$mutational_data %>% names
lapply(cohorts[1:3],function(i) {
  
  x <- sandbox$mutational_data[[i]]
  x <- x[,list(mutation=paste(proteinChange,collapse="|")),by=list(sampleId,hugoGeneSymbol)]
  x <- merge(bowel[SAMPLE_TYPE=="Primary",list(sampleId,GROUP)],x)
  bowel[sampleId %in% x$sampleId,list(M=.N),by=list(GROUP)] -> mss.msi
  x <- x[,.N,by=list(GROUP,hugoGeneSymbol)]
  x <- x[hugoGeneSymbol!=""]
  x <- merge(x,geneLengths,by.x="hugoGeneSymbol",by.y="hgnc_symbol") %>% unique
  
  x <- merge(x,mss.msi,by="GROUP")
  x[,Freq:=N/M]
  x[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]
  
  bySize <- ggplot(x) + aes(cds_length,N/M) +
    geom_point(color="grey80",size=1) +
    geom_smooth(method="glm",
                method.args=list(family="binomial"),fill="lightblue",color="#9999FF") +
    geom_point(color="orange",data=x[hugoGeneSymbol %in% goi],size=3) +
    geom_point(color="red",data=x[hugoGeneSymbol %in% c("NF1","RASA1")],size=3) +
    geom_text_repel(aes(label=hugoGeneSymbol),data=x[hugoGeneSymbol %in% c(goi,"NF1","RASA1")],nudge_x = 0.01,size=3,fontface="italic",max.overlaps=20) +
    myLightTheme +
    theme(text=element_text(size=12),strip.text.x.top = element_text(size=13)) +
    scale_x_log10() +
    xlab("CDS length (log10 scale)") +
    ylab("Mutation Frequency") +
    scale_y_continuous(limits=c(0,1)) +
    facet_wrap(~GROUP)
  
  
}) -> bySizePlots

grid.arrange(arrangeGrob(bySizePlots[[1]] + ggtitle("DFCI Primary CRCs"), top=panelTag("A",x = .5)),
             arrangeGrob(bySizePlots[[2]] + ggtitle("TCGA Primary CRCs"), top=panelTag("B",x = .5)),
             arrangeGrob(bySizePlots[[3]] + ggtitle("MSKCC Primary CRCs"), top=panelTag("C",x = .5)))
            








#  Mutational analysis by codon and specific mutation -------

byCodon <- rbind(
  bowel[,list(STUDY="ALL",GROUP="ALL",.N),by=list(KRAS.CODON)],
  bowel[,list(STUDY="ALL",.N),by=list(GROUP,KRAS.CODON)],
  bowel[,.N,by=list(STUDY,GROUP,KRAS.CODON)] 
)
byCodon[,N2:=sum(N),by=list(STUDY,GROUP)]
byCodon[,Freq:=N/N2*100]

byCodon[KRAS.CODON!="WT",N3:=sum(N),by=list(STUDY,GROUP)]
byCodon[,Freq.mutations:=N/N3*100]

byCodon[GROUP=="ALL"][order(Freq,decreasing = T)]

fig3A <- ggplot(byCodon[GROUP!="ALL"]) + 
  aes(Freq,KRAS.CODON) +
  geom_col(aes(fill=STUDY),show.legend = F) + 
  geom_text(aes(label=sprintf(ifelse(Freq>1," %1.0f%%"," %1.1f%%"),Freq),group=STUDY),hjust = 0,size=3) +
  scale_y_discrete(limits=rev) +
  facet_grid(col=vars(STUDY),row=vars(GROUP)) +
  scale_fill_manual("Study",values = palette$study) +
  xlab("Mutated Cases (%)") +
  ylab("_KRAS_ mutated codon") +
  scale_x_continuous(limits=c(0,100)) +
  mytheme +
  theme(axis.title.y = element_markdown(),strip.text=element_text(size=13)) 

# All mutations 

byMutation <- data.table(bowel)
byMutation[grep("\\|",KRAS),KRAS:="Multiple"]
byMutation <- rbind(
  byMutation[,list(STUDY="ALL",GROUP="ALL",N=.N),by=KRAS],
  byMutation[,list(STUDY="ALL",N=.N),by=list(GROUP,KRAS)],
  byMutation[,list(N=.N),by=list(STUDY,KRAS,GROUP)])

byMutation[,N2:=sum(N),by=list(GROUP,STUDY)]
byMutation[,Freq:=N/N2*100]

byMutation[KRAS!="WT",N3:=sum(N),by=list(GROUP,STUDY)]
byMutation[,Freq.mutations := N/N3*100]

byMutation[GROUP=="ALL" & STUDY=="ALL" & KRAS != "WT" & N > 1][order(Freq,decreasing=T),KRAS] -> mutorder



fig4 <- ggplot(byMutation[KRAS!="WT" & KRAS %in% mutorder & GROUP != "ALL"]) + aes(Freq,factor(KRAS,mutorder)) +
  geom_col(aes(fill=STUDY),show.legend = F) +
  scale_y_discrete(limits=rev) +
  geom_text(aes(label=sprintf(ifelse(Freq>1," %1.0f%%"," %1.1f%%"),Freq)),hjust=0,size=3) +
  facet_grid(row=vars(GROUP),col=vars(STUDY)) +
  scale_fill_manual("Study",values=palette$study) +
  scale_x_continuous(limits=c(0,25),breaks=c(0,10,20)) +
  xlab("Mutated Cases (%)") +
  ylab("_KRAS_ mutation") +
  mytheme +
  theme(axis.title.y = element_markdown(),strip.text=element_text(size=11)) 



# statistics of figure 3
stats3 <- data.table()

for(i in levels(bowel$KRAS.CODON)) {
  
  m1 <- bowel[,table(GROUP,KRAS.CODON==i)] 
  test1 <- fisher.test(m1)[c("estimate","conf.int","p.value")] %>% unlist
  stats3 <- rbind(stats3,data.table(CODON=i,STUDY="ALL",t(test1))) 
  
  for(j in levels(bowel$STUDY)) {
    
    m1 <- bowel[STUDY==j,table(GROUP,KRAS.CODON==i)] 
    test1 <- fisher.test(m1)[c("estimate","conf.int","p.value")] %>% unlist
    stats3 <- rbind(stats3,data.table(CODON=i,STUDY=j,t(test1))) 
    
  }
  
}

names(stats3)[3:5] <- c("OR","CI1","CI2")

stats3[,YulesQ:=(OR-1)/(OR+1)]
stats3[is.infinite(OR),YulesQ:=1]
stats3[,YulesQlow:=(CI1-1)/(CI1+1)]
stats3[,YulesQhigh:=(CI2-1)/(CI2+1)]
stats3[is.infinite(CI2),YulesQhigh:=1]

stats3[,p.adjusted:=p.adjust(p.value,"fdr"),by=STUDY]
stats3[,SigLevel:=cut(p.adjusted,c(0,1e-4,1e-3,1e-2,5e-2,1),c("<0.0001","<0.001","<0.01","<0.05","NS")) %>% fct_rev]
stats3[,CODON:=factor(CODON,c("WT","G12","G13","Q61","K117","A146","Other","Multiple"))]
stats3[,STUDY:=factor(STUDY,c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis"))]

fig3B <- ggplot(stats3) + aes(YulesQ,CODON) +
  geom_vline(xintercept = 0,lty=2) +
  geom_segment(aes(x=YulesQlow,xend=YulesQhigh,yend=CODON,color=SigLevel),lwd=3,show.legend=T) +
  geom_point() +
  geom_text(aes(label=label),data=data.table(YulesQ=c(-1,1),CODON=c("WT"),label=c("MSS","MSI")),hjust=rep(c(0,1),5),size=3) +
  xlab("Yule's Q") +
  ylab("_KRAS_ mutated codon") +
  scale_color_manual("Sig Level",values=palette$SigLevel) +
  scale_y_discrete(limits=rev) +
  facet_grid(cols=vars(STUDY)) +
  theme(axis.title.y.left = element_markdown(),strip.text.x.top = element_text(size=12)) +
  myLightTheme

# figure 3

grid.arrange(arrangeGrob(fig3A + theme(plot.margin = unit(c(1,3,.1,1,1),"cm")),
                         top=textGrob("A",x=unit(1,"cm"),y=unit(-.5,"cm"),gp=gpar(fontsize=20))),
             arrangeGrob(fig3B + theme(plot.margin = unit(c(1,1,1,1),"cm")),
                         top=textGrob("B",x=unit(1,"cm"),y=unit(-.5,"cm"),gp=gpar(fontsize=20))),
             heights=3:2)

stats3[p.adjusted < 0.05 & STUDY=="ALL",sprintf("\n%s %s, OR=%1.2f, CI=[%1.1f-%1.1f], q-value=%1.1e",STUDY,CODON,OR,CI1,CI2,p.adjusted)] %>% cat

# explore factors associated to the incidence of G12 mutations

glm1 <- glm(I(KRAS.CODON=="G12") ~ GROUP + SEX + AGE + STAGE2 + STUDY,bowel[SAMPLE_TYPE=="Primary"],family=binomial) 
summary(glm1)
emmeans(glm1,~STUDY) %>% pairs %>% summary %>% data.table -> G12stats 
G12stats <- cbind(G12stats,emmeans(glm1,~STUDY) %>% pairs %>% confint %>% data.table %>% {.[,c(5,6)]}) 

G12stats[,sprintf("\n%s, OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.1e",contrast,exp(estimate),exp(asymp.LCL),exp(asymp.UCL),p.value)] %>% cat




##### POINTER OCT 28 2024 ############


