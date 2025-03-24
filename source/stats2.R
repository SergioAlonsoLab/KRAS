

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
fig2A <- ggplot(byCodon[GROUP!="ALL"]) + 
  aes(Freq,KRAS.CODON) +
  geom_col(aes(fill=STUDY),show.legend = F) + 
  geom_text(aes(label=sprintf(ifelse(Freq>1," %1.0f%%"," %1.1f%%"),Freq),group=STUDY),hjust = 0,size=3) +
  scale_y_discrete(limits=rev) +
  facet_grid(col=vars(STUDY),row=vars(GROUP)) +
  scale_fill_manual("Study",values = palette$study) +
  xlab("Mutated Cases (%)") +
  ylab("_KRAS_ mutated codon") +
  scale_x_continuous(limits=c(0,100)) +
  theme1 +
  theme(axis.title.y = element_markdown()) 

# statistics of figure 2
stats2 <- data.table()

for(i in levels(bowel$KRAS.CODON)) {
  
  m1 <- bowel[,table(GROUP,KRAS.CODON==i)] 
  test1 <- fisher.test(m1)[c("estimate","conf.int","p.value")] %>% unlist
  stats2 <- rbind(stats2,data.table(CODON=i,STUDY="ALL",t(test1))) 
  
  for(j in levels(bowel$STUDY)) {
    
    m1 <- bowel[STUDY==j,table(GROUP,KRAS.CODON==i)] 
    test1 <- fisher.test(m1)[c("estimate","conf.int","p.value")] %>% unlist
    stats2 <- rbind(stats2,data.table(CODON=i,STUDY=j,t(test1))) 
    
  }
  
}

names(stats2)[3:5] <- c("OR","CI1","CI2")

stats2[,YulesQ:=(OR-1)/(OR+1)]
stats2[is.infinite(OR),YulesQ:=1]
stats2[,YulesQlow:=(CI1-1)/(CI1+1)]
stats2[,YulesQhigh:=(CI2-1)/(CI2+1)]
stats2[is.infinite(CI2),YulesQhigh:=1]

stats2[,p.adjusted:=p.adjust(p.value,"fdr"),by=STUDY]
stats2[,SigLevel:=cut(p.adjusted,c(0,1e-4,1e-3,1e-2,5e-2,1),c("<0.0001","<0.001","<0.01","<0.05","NS")) %>% fct_rev]
stats2[,CODON:=factor(CODON,c("WT","G12","G13","Q61","K117","A146","Other","Multiple"))]
stats2[,STUDY:=factor(STUDY,c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis"))]

fig2B <- ggplot(stats2) + aes(YulesQ,CODON) +
  geom_vline(xintercept = 0,lty=2) +
  geom_segment(aes(x=YulesQlow,xend=YulesQhigh,yend=CODON,color=SigLevel),lwd=3,show.legend=T) +
  geom_point() +
  geom_text(aes(label=label),data=data.table(YulesQ=c(-1,1),CODON=c("WT"),label=c("MSS","MSI")),hjust=rep(c(0,1),5),size=3) +
  xlab("Yule's Q") +
  ylab("_KRAS_ mutated codon") +
  scale_color_manual("Sig Level",values=palette$SigLevel) +
  scale_y_discrete(limits=rev) +
  facet_grid(cols=vars(STUDY)) +
  theme1 +
  theme(axis.title.y.left = element_markdown()) 
  

# figure 2

figure2 <- grid.arrange(arrangeGrob(fig2A + theme(plot.margin = unit(c(0,3,.1,1,1),"cm"),),
                                    top=textGrob("A",x=unit(1,"cm"),y=unit(-.2,"cm"),gp=gpar(fontsize=20))),
                        arrangeGrob(fig2B + theme(plot.margin = unit(c(0,1,1,1),"cm")),
                                    top=textGrob("B",x=unit(1,"cm"),y=unit(-.2,"cm"),gp=gpar(fontsize=20))),
                        heights=c(3,2))

figuresHR(figure2,11,7,fileBaseName = "Figure2")










