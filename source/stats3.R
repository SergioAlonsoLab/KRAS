x <- bowel[,table(KRAS) %>% sort(decreasing=T) %>% names][1:17]


byMutation <- rbind(
  bowel[,list(STUDY="ALL",GROUP="ALL",.N),by=list(KRAS)],
  bowel[,list(STUDY="ALL",.N),by=list(GROUP,KRAS)],
  bowel[,.N,by=list(STUDY,GROUP,KRAS)] 
)


byMutation[,N2:=sum(N),by=list(STUDY,GROUP)]
byMutation[,Freq:=N/N2*100]

byMutation[KRAS!="WT",N3:=sum(N),by=list(STUDY,GROUP)]
byMutation[,Freq.mutations:=N/N3*100]

byMutation <- byMutation[KRAS %in% x]
byMutation[,KRAS:=factor(KRAS,x)]
byMutation$KRAS



fig3A <- ggplot(byMutation[GROUP!="ALL" & !is.na(KRAS)]) + 
  aes(Freq,KRAS) +
  geom_col(aes(fill=STUDY),show.legend = F) + 
  geom_text(aes(label=sprintf(ifelse(Freq>1," %1.0f%%"," %1.1f%%"),Freq),group=STUDY),
            hjust = 0,size=3) +
  scale_y_discrete(limits=rev) +
  facet_grid(col=vars(STUDY),row=vars(GROUP)) +
  scale_fill_manual("Study",values = palette$study) +
  xlab("Mutated Cases (%)") +
  ylab("_KRAS_ mutation") +
  scale_x_continuous(limits=c(0,100)) +
  theme1 +
  theme(axis.title.y = element_markdown()) 

figuresHR(fig3A,10,6,"Figure3")

# association with MSI

stats3a <- data.table()


for(i in x) {
  for(j in c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis")) {
    
    if(j=="ALL") m1 <- bowel[,table(GROUP,factor(KRAS==i,c(F,T))) %>% fisher.test] 
    else
      m1 <- bowel[STUDY==j,table(GROUP,factor(KRAS==i,c(F,T))) %>% fisher.test]
    
    stats3a <- rbind(stats3a,data.table(Mutation=i,STUDY=j,OR=m1$estimate,
                                        CI1=m1$conf.int[1],CI2=m1$conf.int[2],
                                        p.value=m1$p.value))
    
  }
  
  
}

stats3a[,STUDY:=factor(STUDY,unique(STUDY))]
stats3a[,YulesQ := (OR-1)/(OR+1)]
stats3a[,YulesCI1 := (CI1-1)/(CI1+1)]
stats3a[,YulesCI2 := (CI2-1)/(CI2+1)]
stats3a[is.na(YulesCI2),YulesCI2 := -1]
stats3a[,Sig:=cut(p.value,c(0,.0001,.001,.01,.05,1),labels=names(palette$SigLevel))]
stats3a[,Sig:=factor(Sig,names(palette$SigLevel))]
stats3a[,Sig:=fct_rev(Sig)]
stats3a[,Mutation:=factor(Mutation,x)]
stats3a[CI1==0 & is.infinite(CI2),`:=`(YulesQ=NA,YulesCI1=NA,YulesCI2=NA)]


fig3B <- ggplot(stats3a) + aes(YulesQ,Mutation) +
  geom_vline(xintercept = 0,lty = 2) +
  geom_segment(aes(x=YulesCI1,xend=YulesCI2,color=Sig),lwd=3,show.legend = T) +
  geom_point() +
  scale_y_discrete(limits=rev) +
  scale_color_manual("Sig Level",values=palette$SigLevel,drop=F) +
  coord_cartesian(ylim=c(1,18)) +
  annotate("text",x=-1,y=18,label="MSS",hjust=0,size=3) +
  annotate("text",x=1,y=18,label="MSI",hjust=1,size=3) +
  xlab("Yule's Q") +
  ylab("_KRAS_ mutation") +
  theme1 +
  theme(axis.title.y = element_markdown()) +
  facet_grid(cols=vars(STUDY)) 


fig3 <- grid.arrange(arrangeGrob(fig3A + theme(plot.margin = unit(c(0,2,0,.5),"cm"),strip.text.x.top = element_text(size=12)),top=panelTag("A")),
             arrangeGrob(fig3B + theme(plot.margin = unit(c(0,0,.5,.5),"cm"),strip.text.x.top = element_text(size=12)),top=panelTag("B")),
             heights=c(2,1.25))

figuresHR(fig3,11,9,"Figure3")


# multivariate 

stats3b <- data.table()


for(i in x) {
  print(i)
  m1 <- glm(I(KRAS==i) ~ GROUP + SEX + STAGE2 + STUDY,bowel[SAMPLE_TYPE=="Primary"],family=binomial)
  

  m1 <- list(emmeans(m1,~SEX),
             emmeans(m1,~GROUP),
             emmeans(m1,~STAGE2),
             emmeans(m1,~STUDY))
  
  m1 <- lapply(m1,pairs)
  
  m1 <- cbind(lapply(m1,summary,adjust="fdr") %>% rbindlist,
              lapply(m1,function(x) confint(x)[,5:6]) %>% rbindlist)
  
  
  stats3b <- rbind(stats3b,data.table(Mutation=i,m1))

}

stats3b[,contrast := gsub(" - "," ↔ ",contrast)]
stats3b[,contrast := gsub("-","/",contrast)]
stats3b[,contrast := gsub("\\(|\\)","",contrast)]

# for graphical representation preferable to reverse the log OR

stats3b[,OR := exp(-estimate)]
stats3b[,CI1 := exp(-asymp.UCL)]
stats3b[,CI2 := exp(-asymp.LCL)]

# Transform to YulesQ

stats3b[,YulesQ := (OR-1)/(OR+1)]
stats3b[,YulesCI1 := (CI1-1)/(CI1+1)]
stats3b[,YulesCI2 := (CI2-1)/(CI2+1)]
stats3b[is.infinite(OR),YulesQ:=1]
stats3b[is.infinite(CI1),YulesCI1:=1]
stats3b[is.infinite(CI2),YulesCI2:=1]

stats3b[,contrast := factor(contrast,contrast[1:6])]

stats3b[,Sig:=cut(p.value,c(0,.0001,.001,.01,.05,1),labels=names(palette$SigLevel))]
stats3b[,Sig:=fct_rev(Sig)]

x <- bowel[,table(KRAS) %>% sort(decreasing=T) %>% names][2:11]
x <- stats3b[Mutation %in% x][,Mutation := factor(Mutation,x)]
fig3.1 <- ggplot(x) + aes(YulesQ,contrast) + 
  geom_vline(xintercept=0,lty=2) +
  geom_segment(aes(x=YulesCI1,xend=YulesCI2,color=Sig),lwd=3,show.legend = T) +
  geom_point() +
  xlab("Yule's Q") +
  ylab(NULL) +
  scale_y_discrete(limits=rev) +
  scale_color_manual("Sig Level",values=palette$SigLevel) +
  theme1 +
  theme(strip.text.x.top = element_markdown()) +
  facet_wrap(~ Mutation,labeller = labeller(Mutation = function(x) paste("_KRAS_",x)),ncol=5) 

figuresHR(fig3.1,10,4,"Figure3.1")



