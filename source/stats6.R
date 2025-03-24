stats6a <- data.table()

for(i in selected.genes) {
  
  m <- bowel[SAMPLE_TYPE=="Primary",table(KRAS.CODON,factor(get(i)!="WT",c(F,T)),GROUP)]

  p1 <- tryCatch(fisher.test(m[c("G12","G13","Q61","K117","A146"),,1])$p.value,error = function(e) 1)
  p2 <- tryCatch(fisher.test(m[c("G12","G13"),,1])$p.value,error = function(e) 1)
  p3 <- tryCatch(fisher.test(m[c("G12","G13","Q61","K117","A146"),,2])$p.value,error = function(e) 1)
  p4 <- tryCatch(fisher.test(m[c("G12","G13"),,2])$p.value,error = function(e) 1)
  
  stats6a <- rbind(stats6a,data.table(Gene=i,GROUP="MSS",p.all=p1,p.g12g13=p2))
  stats6a <- rbind(stats6a,data.table(Gene=i,GROUP="MSI/HYPER",p.all=p3,p.g12g13=p4))
  
}


melt(bowel,id.vars=c("SAMPLE_TYPE","GROUP","KRAS.CODON"),measure.vars = selected.genes,variable.name = "Gene") -> stats6
stats6[,list(MUT=sum(value!="WT"),.N),by=list(SAMPLE_TYPE,GROUP,KRAS.CODON,Gene)] -> stats6
stats6 <- stats6[KRAS.CODON %in% c("WT","G12","G13","Q61","K117","A146") & SAMPLE_TYPE=="Primary"]
stats6[,KRAS.CODON:=fct_reorder(KRAS.CODON,-N)]

stats6 <- merge(stats6,stats6a,by=c("GROUP","Gene"))

stats6[,p.nice:=sprintf("p-value=%1.2g",p.all)]


layers6 <- list( aes(MUT/N*100,KRAS.CODON), 
                 geom_col(aes(fill=KRAS.CODON),show.legend=F),
                 scale_y_discrete(limits=rev,labels=function(x) sprintf("_KRAS_ %s",x)),
                 coord_cartesian(ylim = c(-.5,6)),
                 scale_x_continuous(limits=c(0,100)),
                 geom_text(aes(label=sprintf(ifelse(MUT/N>0.01 | MUT==0," %1.0f%% (%i/%i)"," %1.1f%% (%i/%i)"),MUT/N*100,MUT,N)),hjust=0,data=function(x) subset(x,MUT/N<.5)),
                 geom_text(aes(label=sprintf("%1.0f%% (%i/%i) ",MUT/N*100,MUT,N)),hjust=1,color="black",data=function(x) subset(x,MUT/N>=.5)),
                 geom_text(aes(x=50,y=-.25,label=p.nice),data=function(x) unique(x[,list(GROUP,Gene,p.nice)])),
                 
                 theme1,
                 theme(strip.text.x.top = element_markdown(size=13),
                       text=element_text(size=12),
                       axis.text.y = element_markdown(size=12),
                       legend.title = element_markdown()),
                 ylab(""),
                 xlab("Mutated cases (%)"),
                 scale_fill_manual("_KRAS_ Mutation",values=c("#AAAAFF",rep("#FFAA33",5))),
                 facet_wrap(~ Gene,labeller = labeller(Gene = function(x) sprintf("_%s_ mutations",x))))


x <- c("APC","PIK3CA","FBXW7","TP53","BRAF","NRAS")
x <- stats6[Gene %in% x][,Gene:=factor(Gene,x)][,p.value:=p.all]

fig6A <- x[GROUP=="MSS"] %>% ggplot() + layers6
fig6B <- x[GROUP=="MSI/HYPER"] %>% ggplot() + layers6
fig6.1 <- grid.arrange(arrangeGrob(fig6A + ggtitle("MSS Primary CRCs"),top=panelTag("A")),
             arrangeGrob(fig6B + ggtitle("MSI/hypermutated Primary CRCs"),top=panelTag("B")))

figuresHR(fig6.1,8,10,"Figure6.1")


x <- stats6a[order(p.all)][p.all < 0.01,Gene] 
x <- stats6[Gene %in% x][,Gene:=factor(Gene,x)][,p.value:=p.all]

fig6A <- x[GROUP=="MSS"] %>% ggplot() + layers6
fig6B <- x[GROUP=="MSI/HYPER"] %>% ggplot() + layers6 

fig6.2 <- grid.arrange(arrangeGrob(fig6A + ggtitle("MSS Primary CRCs"),top=panelTag("A")),
                     arrangeGrob(fig6B + ggtitle("MSI/hypermutated Primary CRCs"),top=panelTag("B")))


figuresHR(fig6B + ggtitle("MSI/hypermutated Primary CRCs"),9,7,"Figure6.2")


x <- c("NF1","RASA1")
x <- stats6[Gene %in% x][,Gene:=factor(Gene,x)][,p.value:=p.all]

fig6A <- x[GROUP=="MSS"] %>% ggplot() + layers6
fig6B <- x[GROUP=="MSI/HYPER"] %>% ggplot() + layers6

fig6.3 <- grid.arrange(arrangeGrob(fig6A + ggtitle("MSS Primary CRCs"),top=panelTag("A")),
                     arrangeGrob(fig6B + ggtitle("MSI/hypermutated Primary CRCs"),top=panelTag("B")))

figuresHR(fig6.3,6,6,"Figure6.3")

