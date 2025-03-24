x <- bowel[,..selected.genes] != "WT"
x <- cbind(bowel[,list(GROUP,SAMPLE_TYPE)],x)
x <- x[SAMPLE_TYPE=="Primary"]

stats5 <- data.table()

for(i in c(goi)) {
  for(j in c(selected.genes)) {
    
    m1 <- x[,table(factor(get(i),c(F,T)),factor(get(j),c(F,T)),GROUP)]
    ft1 <- fisher.test(m1[,,1])[c("estimate","conf.int","p.value")] %>% unlist
    ft2 <- fisher.test(m1[,,2])[c("estimate","conf.int","p.value")] %>% unlist
    
    m1 <- data.table(rbind(c(ft1,m1[,,1]),c(ft2,m1[,,2])))
    m1[,Gene1:=i]
    m1[,Gene2:=j]
    m1[,GROUP:=c("MSS","MSI/HYPER")]
    
    stats5 <- rbind(stats5,m1)
    
  }
}

names(stats5)[5:8] <- c("WW","MW","WM","MM")
names(stats5)[1:3] <- c("OR","CI1","CI2")
stats5 <- stats5[Gene1!=Gene2]

stats5[,YulesQ:=(OR-1)/(OR+1)]
stats5[,YulesQlow:=(CI1-1)/(CI1+1)]
stats5[,YulesQhigh:=(CI2-1)/(CI2+1)]
stats5[is.infinite(OR),YulesQ:=1]
stats5[is.infinite(CI2),YulesQ:=1]

stats5[,Gene1:=factor(Gene1,c(goi,"POLE","KMT2D","RNF43"))]
stats5[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]


stats5[,p.adjusted:=p.adjust(p.value,"fdr"),by=list(GROUP,Gene1)]

stats5[,SigLevel:=cut(p.adjusted,c(0,1e-4,1e-3,1e-2,5e-2,1),c("<0.0001","<0.001","<0.01","<0.05","NS"))]
stats5[,SigLevel:=fct_rev(SigLevel)]

stats5[,MinM:=MM/(MM+MW)]
stats5[,MinW:=WM/(WM+WW)]

layers5 <-  list(aes(MinM*100,YulesQ) , 
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
                    theme1 ,
                    theme(strip.text.x.top = element_text(face="italic",size=13),
                          strip.text.y.right = element_text(size=13),
                          axis.title = element_text(size=13)))



fig5A <- ggplot(stats5[order(p.adjusted,decreasing = T)][Gene1 %in% c("KRAS","NRAS","HRAS")]) + layers5

fig5B <- stats5[p.adjusted < 0.05 & Gene1 %in% c("KRAS","NRAS","HRAS")][order(GROUP,Gene1,p.adjusted)][,list(GROUP,Gene1,Gene2,
                                                                                                             "OR [CI95%]"=sprintf("%1.2f [%1.2f-%1.2f]",OR,CI1,CI2),
                                                                                                             FDR=sprintf("%1.1e",p.adjusted),
                                                                                                             "co-Mutated"=sprintf("%1.1f%%",MinM*100),
                                                                                                             "Mutated in WT"=sprintf("%1.1f%%",MinW*100))] 
fig5B <- tableGrob(fig5B,rows=NULL,
                   theme=ttheme_default(base_size=10),
                   widths=unit(c(2,2,2,3,2,2,2)*1.5,"cm"),
                   heights=unit(rep(.55,nrow(fig5B)),"cm"))

# Figure 5

fig5.1 <- grid.arrange(arrangeGrob(fig5A + theme(plot.margin=unit(c(.2,.6,.2,1),"cm")),top=panelTag("A")),
                     arrangeGrob(fig5B,top=panelTag("B")),
                     heights=c(3,2.2))

figuresHR(fig5.1,10,8,fileBaseName = "Figure5.1")

# Figure 5

fig5.2 <- grid.arrange(
  ggplot(stats5[order(p.adjusted,decreasing = T)][Gene1 %in% c("APC","TP53","FBXW7")]) + layers5,
  ggplot(stats5[order(p.adjusted,decreasing = T)][Gene1 %in% c("SMAD4","PIK3CA","BRAF")]) + layers5
)

figuresHR(fig5.2,9,9,"Figure5.2")