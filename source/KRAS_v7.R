# Load required libraries ----

library(data.table)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggforce) # for geom_sina
library(ggpubr)  # nicer way of integrating tables
library(grid)    
library(gridExtra)
library(patchwork) # combine ggplots
library(emmeans)

# Note: The original multipanel figures were created with gridExtra.
# We are updating the code to use patchwork, which offers greater flexibility,
# simpler syntax, and improved figure quality by automatically aligning plotting areas.


# list of figures (to be updated before publishing)

# Figure 1: Mutational frequency of Ras genes and other frequently mutated genes in CRC

# Figure 2: Structure of GDP-bound wild type KRAS (mol*)
# Figure 3: KRAS codon-specific mutational frequency in CRC
# Figure 4: KRAS mutational frequency in CRC
# Figure 5: KRAS, NRAS and HRAS co-mutation analysis in primary CRCs
# Figure 6: Co-mutation analysis of KRAS G12 and KRAS G13 mutant primary CRCs
# Figure 7: Co-mutations of KRAS with NF1 and RASA1
# Figure 8: Clinical trials (in a different script)

# Figure S1: Associations of mutations in other genes and clinical parameters
# Figure S2: Table with samples analysed from the same patient
# Figure S3: Co-mutation analysis of other genes
# Figure SX: Co-mutation analysis of KRAS Q61, A146 and K117
# Figure SX: 



# Graphical parameters and functions look & feel ----

palette <- list()

palette$study <- c("#C19F08","#369E4B","#357EB9","#FA792F","#ce3206")
names(palette$study) <- c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis")

palette$SigLevel <- c("<0.0001"="#c32f27",
                      "<0.001"="#d8572a",
                      "<0.01"="#db7c26",
                      "<0.05"="#f7b538",
                      "NS"="grey75")

theme1 <- theme(panel.background = element_rect(fill="grey95"),
                strip.background = element_rect(fill="grey85",colour="grey50"),
                panel.grid = element_line(color="white"),
                panel.border = element_rect(color="grey50",fill=NA,linewidth=rel(1)),
                strip.text = element_text(size=12))

# panel labeller (a function to add a letter in the left corner of a graph)

panelTag <- function(tag="A",x=1,y=0,fsize=20) {
  textGrob(label=tag,x=unit(x,"cm"),y=unit(y,"cm"),gp = gpar(fontsize=fsize))
}

# high resolution pdf, tiff and png figures
# it should be changed to ggsave in a future version

figuresHR <- function(gplot,width=10,height=10,fileBaseName="noname",outputDir="output",res=300,...) {
  pdf(sprintf("%s/%s.pdf",outputDir,fileBaseName),width,height)
  plot(gplot)
  dev.off()
  
  tiff(sprintf("%s/%s.tiff",outputDir,fileBaseName),width,height,units = "in",res=res,compression = "lzw",type="cairo",...)
  plot(gplot)
  dev.off()
  
  png(sprintf("%s/%s.png",outputDir,fileBaseName),width,height,units = "in",res=res,...)
  plot(gplot)
  dev.off()
}



# Load the bowel data ----

bowel <- readRDS("data/bowel.rds")
genes <- fread("data/all.genes.length.csv")
selected.genes <- intersect(genes$hgnc_symbol,colnames(bowel))
bowel[,KRAS.CODON:=factor(KRAS.CODON,c("WT","G12","G13","Q61","K117","A146","Other","Multiple"))]
bowel[,SAMPLE_TYPE:=factor(SAMPLE_TYPE,c("Primary","Metastasis"))]
goi <- c("APC","TP53","FBXW7","SMAD4","KRAS","NRAS","HRAS","PIK3CA","BRAF")

# create an ouput directory, if it does not exists

if(!file.exists("output")) dir.create("output")

# STATS 0: TMB in different cohorts -----

ggplot(bowel) + aes(STUDY,TMB_NONSYNONYMOUS) + 
  geom_violin(aes(fill=STUDY),color=NA,alpha=.75,show.legend=F) + 
  geom_sina(aes(color=GROUP,group=STUDY),size=.5,maxwidth=.4) +
  annotate("segment",x=1:4-.3,xend=1:4+.3,y=c(14,14,25,25),yend=c(14,14,25,25)) +
  annotate("text",x=1:4,y=rep(.13,4),label=paste0("n=",table(bowel$STUDY))) +
  
  scale_fill_manual(values=palette$study) +
  scale_color_manual("GROUP",values=c("grey20","red3")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  xlab(NULL) +
  ggtitle("Tumor Mutational Burden") +
  ylab("TMB (non synonymous mutations / Mb)") +
  scale_y_log10() + theme1 -> fig0


figuresHR(fig0,fileBaseName = "Supplementary_Figure_S1",width=7,height=5)


# STATS 1: Ras vs primary / metastasis ------

STATS1 <- new.env()
source("source/stats1.R",local=STATS1,echo=T)

# STATS 2: Frequency by codon -------

STATS2 <- new.env()
source("source/stats2.R",local=STATS2,echo = T)

# STATS 3: Frequency by mutation -----

STATS3 <- new.env()
source("source/stats3.R",local=STATS3,echo = T)

# STATS 4: Paired primary - metastases samples. -----

STATS4 <- new.env()
source("source/stats4.R",local=STATS4,echo = T)

# STATS 5: Co-mutation analysis by gene ----

STATS5 <- new.env()
source("source/stats5.R",local=STATS5,echo = T)

# STATS 6: co-mutational analysis by KRAS mutated codon ----

STATS6 <- new.env()
source("source/stats6.R",local=STATS6,echo = T)

# STATS 7: Mutation analysis by gene size ------

STATS7 <- new.env()
source("source/stats7.R",local=STATS7,echo = T)


