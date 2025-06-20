x <- melt(bowel,id.vars = c("SEX","GROUP","SAMPLE_TYPE","STUDY","STAGE2"),measure.vars = selected.genes,variable.name = "GENE",value.name = "STATUS")
x[,SAMPLE_TYPE:=factor(SAMPLE_TYPE,c("Primary","Metastasis"))]


fig1A <- x[GENE %in% goi & !is.na(STAGE2),list(MUT=sum(STATUS!="WT"),.N),by=list(GROUP,STAGE2,GENE)][,GENE:=factor(GENE,goi)][,FREQ:=MUT/N*100] %>% 
  ggplot() + aes(FREQ,GENE) + 
  geom_col(aes(fill=STAGE2),position=position_dodge2(reverse=TRUE)) +
  geom_text(aes(group=STAGE2,label=sprintf(ifelse(FREQ==0 | FREQ > 1," %1.0f%%"," %1.1f%%"),FREQ)),
            position=position_dodge2(reverse=TRUE,width=.9),hjust=0,size=3.5) +
  scale_y_discrete(limits=rev) +
  facet_grid(cols=vars(GROUP)) +
  scale_x_continuous(limits=c(0,100)) +
  ylab("") + xlab("Percentage of mutated cases") +
  theme1 +
  theme(text=element_text(size=12),axis.text.y=element_text(face="italic",size=13)) +
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
stats1[,SigLevel:=cut(p.value,c(0,1e-4,1e-3,1e-2,5e-2,1),labels=c("<0.0001","<0.001","<0.01","<0.05","NS"))]
stats1[,SigLevel := fct_rev(SigLevel)]

# for graphical purposes, it's preferable to switch the sign of the estimate and CIs

stats1[,estimate:=-estimate]
stats1[,foo:=-asymp.LCL]
stats1[,asymp.LCL:=-asymp.UCL]
stats1[,asymp.UCL:=foo]
stats1[,foo:=NULL]


layers1 <- list(aes(estimate,contrast),
                geom_vline(xintercept = 0,lty=2),
                geom_segment(aes(x=asymp.LCL,xend=asymp.UCL,yend=contrast,color=SigLevel),lwd=3,show.legend=T),
                geom_point(),
                scale_y_discrete(limits=rev),
                ylab(""),
                xlab("log Odds Ratio"),
                facet_grid(rows=vars(Model),cols=vars(Gene),scales = "free"),
                scale_color_manual("Sig Level",values=palette$SigLevel,
                                   drop=FALSE),
                theme1,
                theme(strip.text.x.top = element_text(face="italic",size=13),
                      strip.text.y.right = element_text(size=11),
                      axis.text.y = element_text(size=12)))


fig1B <- stats1[Gene %in% c("KRAS","NRAS","HRAS")] %>% ggplot() + layers1


# Figure stats 1  


fig1 <- fig1A / fig1B + plot_layout(heights=c(2,1.1)) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=20),
        plot.tag.position = c(.05, .98))


figuresHR(fig1,8,8,"Figure1")

# Supplementary Figure 1

figS1 <- grid.arrange(stats1[Gene %in% c("APC","TP53","FBXW7")] %>% ggplot() + layers1,
                      stats1[Gene %in% c("SMAD4","PIK3CA","BRAF")] %>% ggplot() + layers1)

figuresHR(figS1,10,8,"FigureS1")

pdf("output/FigureS1.pdf",10,8)
plot(figS1)
dev.off()

# some stats shown in the body of the paper

bowel[,table(KRAS!="WT",SAMPLE_TYPE)]
stats1[Gene=="KRAS",sprintf("\n%s OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.2g",contrast,exp(estimate),exp(asymp.LCL),exp(asymp.UCL),p.value)] %>% cat()
stats1[Gene=="KRAS",sprintf("\n%s OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.2g",contrast,exp(-estimate),exp(-asymp.UCL),exp(-asymp.LCL),p.value)] %>% cat()
stats1[Gene=="HRAS",sprintf("\n%s OR=%1.1f, CI=[%1.1f-%1.1f], p=%1.2g",contrast,exp(estimate),exp(asymp.LCL),exp(asymp.UCL),p.value)] %>% cat()







