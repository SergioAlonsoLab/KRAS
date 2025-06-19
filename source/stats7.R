library(patchwork)


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
    theme1 +
    theme(text=element_text(size=12)) +
    scale_x_log10() +
    xlab("CDS length (log10 scale)") +
    ylab("Mutation Frequency") +
    scale_y_continuous(limits=c(0,1)) +
    facet_wrap(~GROUP)
  
  
}) -> bySizePlots


bySizePlots[[1]] + ggtitle("DFCI Primary CRCs") +
  bySizePlots[[2]] + ggtitle("TCGA Primary CRCs") +
  bySizePlots[[3]] + ggtitle("MSKCC Primary CRCs") +
  plot_layout(ncol=1) &
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size=20),
        plot.tag.position = c(.01, .98)) -> figS7
  

figuresHR(figS7,8,12,"FigureS7")

# validation in MS-CHORD

x <- sandbox$clinical_data$msk_chord_2024[CANCER_TYPE=="Colorectal Cancer"]
y <- sandbox$mutational_data$msk_chord_2024[sampleId %in% x$sampleId]



y2 <- dcast(y[hugoGeneSymbol %in% selected.genes],sampleId ~ hugoGeneSymbol,value.var = "proteinChange",fun.aggregate = function(x) paste(sort(x),collapse="|"))
y2[y2==""] <- "WT"
x <- merge(x,y2)
x[,MUTATION_COUNT:=as.numeric(MUTATION_COUNT)]
x[,MSI_SCORE:=as.numeric(MSI_SCORE)]
ggplot(x) + aes(x$MSI_TYPE,x$MSI_SCORE) + geom_boxplot()

x$KRAS.CODON <- character()

for(i in c("WT","G12","G13","Q61","A146","K117")) {
  
  x[grep(i,KRAS),KRAS.CODON:=i]
  
}

x[grep("\\|",KRAS),KRAS.CODON:="Multiple"]
x[is.na(KRAS.CODON),KRAS.CODON:="Other"]

x[,GROUP:="MSS"]
x[MSI_SCORE > 10,GROUP:="MSI/HYPER"]

ggplot(x) + geom_density(aes(x=MUTATION_COUNT %>% log10,fill=GROUP),alpha=.3) + geom_vline(xintercept = log10(25))

bowel[,PanelGenes:=gsub("[A-Z]","",GENE_PANEL) %>% as.numeric]
x[,PanelGenes:=gsub("[A-Z]","",GENE_PANEL) %>% as.numeric]
bowel[,MUTATION_COUNT := as.numeric(MUTATION_COUNT)]
bowel[STUDY %in% c("MSKCC Primary","MSKCC Metastasis"),lm(TMB_NONSYNONYMOUS ~ PanelGenes + MUTATION_COUNT)] -> lm1
summary(lm1)

x[,TMB := predict(lm1,type="response",newdata=.SD)]
ggplot(x) + aes(MUTATION_COUNT,TMB) + geom_point(aes(color=SAMPLE_TYPE))

table(x$TMB >= 25,x$GROUP,useNA="if")

x[MUTATION_COUNT>25,GROUP:="MSI/HYPER"]
x[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]



stats7b <- data.table()

for(i in c("KRAS","NRAS","HRAS")) {
  for(j in selected.genes) {
   m <- x[SAMPLE_TYPE=="Primary",table(factor(get(i)!="WT",c(F,T)),
                 factor(get(j)!="WT",c(F,T)),GROUP)] 
   
   t1 <- fisher.test(m[,,1])
   t2 <- fisher.test(m[,,2])
   
   index1 <- m[2,2,1]/sum(m[2,,1])
   index2 <- m[2,2,2]/sum(m[2,,2])
   
   stats7b <- rbind(stats7b,data.table(Gene1=i,Gene2=j,GROUP="MSS",OR=t1$estimate,CI1=t1$conf.int[1],CI2=t1$conf.int[2],p.value=t1$p.value,index=index1))
   stats7b <- rbind(stats7b,data.table(Gene1=i,Gene2=j,GROUP="MSI/HYPER",OR=t2$estimate,CI1=t2$conf.int[1],CI2=t2$conf.int[2],p.value=t2$p.value,index=index2))
   
  }
}



stats7b <- stats7b[Gene1 != Gene2]
stats7b[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]
stats7b[,Gene1:=factor(Gene1,c("KRAS","NRAS","HRAS"))]
stats7b[,YulesQ:=(OR-1)/(OR+1)]
stats7b[,YulesCI1:=(CI1-1)/(CI1+1)]
stats7b[,YulesCI2:=(CI2-1)/(CI2+1)]

stats7b[,FDR:=p.adjust(p.value,"fdr"),by=list(GROUP,Gene1)]
stats7b[,Sig:=cut(FDR,c(0,0.0001,0.001,0.01,0.05,1),labels=names(palette$SigLevel))]
stats7b[,Sig:=fct_rev(Sig)]

fig7.1 <- ggplot(stats7b) + aes(index*100,YulesQ) +
  geom_hline(yintercept = 0,lty=2) +
  geom_point(size=1,color="grey80") +
  geom_segment(aes(y=YulesCI1,yend=YulesCI2,color=Sig),data=stats7b[FDR < 0.05],show.legend=T,lwd=1) +
  
  geom_point(aes(color=Sig),data=stats7b[FDR < 0.05],show.legend=T,size=2) +
  geom_text_repel(aes(label=Gene2),data=stats7b[FDR < 0.05],fontface="italic",size=2.5,nudge_x = 2,max.overlaps = 20) +
  facet_grid(cols=vars(Gene1),rows=vars(GROUP)) +
  scale_color_manual("Sig Level",values=palette$SigLevel,drop=F) +
  xlab("Co-mutation index (%)") +
  scale_x_continuous(limits=c(0,100)) +
  theme1 +
  theme(strip.text.x.top = element_text(face="italic")) +
  ggtitle("Co-mutation analysis in MSK-CHORD primary CRCs")


stats7b <- data.table()
x$KRAS.CODON %>% table
x2 <- x[SAMPLE_TYPE=="Primary",list(NF1=sum(NF1!="WT"),RASA1=sum(RASA1!="WT"),N=.N),by=list(GROUP,KRAS.CODON)]
x2[,KRAS.CODON := factor(KRAS.CODON,c("WT","G12","G13","A146","Q61","K117"))]
x2 <- melt(x2,id.vars = c("GROUP","N","KRAS.CODON")) %>% na.omit

x2.pvals <- x2[KRAS.CODON!="WT",list(p.value=fisher.test(cbind(N,value))$p.value),by=list(GROUP,variable)]

fig7.2 <- ggplot(x2) + aes(value/N*100,KRAS.CODON) + 
  coord_cartesian(ylim = c(-.5,6)) +
  geom_col(aes(fill=KRAS.CODON),show.legend=F) +
  geom_text(aes(label=sprintf(ifelse(value==0|value/N>0.01,"%1.0f%% (%i/%i)","%1.1f%% (%i/%i"),value/N*100,value,N)),hjust=0,nudge_x = 5) + 
  geom_text(aes(label=sprintf("p-value=%1.2f",p.value)),x=50,y=-.25,data=x2.pvals) +
  facet_grid(cols=vars(variable),rows=vars(GROUP)) +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(limits=c(0,100)) +
  ylab("") +
  xlab("Mutated cases (%)") +
  scale_fill_manual(values=c("#AAAAFF",rep("#FFAA33",5))) +
  theme1 +
  theme(strip.text.x.top = element_text(face="italic")) +
  ggtitle("KRAS vs NF1/RASA1 co-mutation analysis in MSK-CHORD Primary CRCs")

fig7.3 <- fig7.1 + fig7.2 + plot_layout(ncol=1) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size=20),
        plot.tag.position = c(.02, .98))

figuresHR(fig7.3,8,10,"Figure7.1")
             
merge(x,y[,list(MUTS=.N),by=patientId]) %>% ggplot() + aes(MUTATION_COUNT %>% as.numeric,MUTS) + geom_point()

y[,MUTS:=.N,by=patientId]
y[,GROUP:="MSS"]
y[MUTS >= 25,GROUP:="MSI/HYPER"]

y2 <- dcast(y, patientId + GROUP ~ hugoGeneSymbol,value.var = "proteinChange", fun.aggregate = function(x) paste(sort(x),collapse="|"))
y2[y2==""] <- "WT"
y2[,M:=.N,by=GROUP]
y2 <- melt(y2,id.vars = c("patientId","GROUP","M"),variable.name = "Gene",value.name = "Status")
y2 <- y2[,list(N=.N),by=list(GROUP,hgnc_symbol=Gene,M,g=ifelse(Status=="WT","WT","MUT"))]
y2 <- merge(y2,geneLengths)
y2[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]

ggplot(y2[g=="MUT"]) +
  aes(cds_length,N/M) +
  geom_point(color="grey80",size=1) +
  geom_smooth(method="glm",
              method.args=list(family="binomial"),fill="lightblue",color="#9999FF") +
  geom_point(color="orange",data = . %>% filter(hgnc_symbol %in% goi),size=3) +
  geom_point(color="red",data= . %>% filter(hgnc_symbol %in% c("NF1","RASA1")),size=3) +
  geom_text_repel(aes(label=hgnc_symbol),data= . %>% filter(hgnc_symbol %in% c(goi,"NF1","RASA1")),nudge_x = 0.01,size=3,fontface="italic",max.overlaps=20) +
  theme1 +
  theme(text=element_text(size=12)) +
  scale_x_log10() +
  xlab("CDS length (log10 scale)") +
  ylab("Mutation Frequency") +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(~GROUP) +
  ggtitle("Mutation Frequency vs CDS length MSK-CHORD Primary CRCs") -> fig7.4


figuresHR(fig7.4,8,4,"Figure7.2")
  
