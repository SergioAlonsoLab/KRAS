library(data.table)
library(ggplot2)
library(survival)
library(survminer)

mutations <- fread("~/Downloads/mutations.txt")
clinical <- fread("~/Downloads/combined_study_clinical_data.tsv")

mutations <- mutations[order(SAMPLE_ID)]
clinical <- clinical[order(`Sample ID`)]

all(mutations$SAMPLE_ID == clinical$`Sample ID`)

d0 <- data.table(clinical,mutations)

all(d0$`Sample ID` == d0$SAMPLE_ID)

# remove the MSK study

d0 <- d0[`Study ID` != "crc_apc_impact_2020"]

any(duplicated(d0$`Sample ID`))


# remove samples withut mutational information

d0 <- d0[KRAS != "NS"]
d0$Study <- d0$`Study ID`
d0$Study <- factor(d0$Study) 
levels(d0$Study) <- c("DFCI","TCGA","MSKCC")

d0[Study == "DFCI", `Sample Type` := "Primary"]

# Include stage

d0[is.na(`Stage At Diagnosis`)]

d0[!is.na(`Stage At Diagnosis`),Stage := `Stage At Diagnosis`]
d0[!is.na(`Tumor Stage`),Stage := `Tumor Stage`]
d0[!is.na(`Neoplasm Disease Stage American Joint Committee on Cancer Code`),
   Stage := `Neoplasm Disease Stage American Joint Committee on Cancer Code` %>% 
     gsub("STAGE ","",.) %>% gsub("[ABC ]","",.)]
d0[Study == "MSKCC" & `Sample Type` == "Metastasis",Study := "MSKCC-met"]

x <- table(d0$Study,d0$Stage,useNA="if")
x

(x/rowSums(x) * 100) %>% round(.,0)

d0$STUDY_ID <- NULL
fwrite(d0,"~/Desktop/KRAS.DATA.csv")
fread("~/Desktop/KRAS.DATA.csv")


# sex
x <- table(d0$Study,d0$Sex)
(x/rowSums(x)*100) %>% round(.,0)
x <- table(d0$`Study ID`,d0$Sex)
x/rowSums(x)

ggplot()
# age
table(d0$Study,is.na(d0$`Age at Diagnosis`))
table(d0$Study,is.na(d0$`Diagnosis Age`))

d0[is.na(`Age at Diagnosis`), `Age at Diagnosis` := `Diagnosis Age`]
tapply(d0$`Age at Diagnosis`,d0$`Study ID`,mean,na.rm=T) %>% round(.,0)

TukeyHSD(aov(d0$`Age at Diagnosis` ~ d0$Study))

tapply(d0$`Age at Diagnosis`,d0$Study,mean,na.rm=T)
tapply(d0$`Age at Diagnosis`,d0$`Study ID`,mean,na.rm=T)
tapply(d0$`Age at Diagnosis`,d0$`Study ID`,sd,na.rm=T)

# stage

x <- table(d0$`Study ID`,d0$Stage,useNA="if")
round(x/rowSums(x) * 100,0)

# mark hypermutant samples

table(d0$`MSI Status`)
d0[,MSI := NA]
d0[(`TMB (nonsynonymous)` <= 14 & Study %in% c("TCGA","DFCI")) | 
       (`TMB (nonsynonymous)` <= 25 & Study %in% c("MSKCC","MSKCC-met")),
     MSI := "MSS"]


d0[`MSI Status` %in% c("MSI","MSI-high") | 
     (`TMB (nonsynonymous)` > 14 & Study %in% c("TCGA","DFCI")) | 
     (`TMB (nonsynonymous)` > 25 & Study %in% c("MSKCC","MSKCC-met")),
   MSI := "MSI/HYPER"]


msi <- d0[MSI=="MSI/HYPER"]
mss <- d0[MSI=="MSS"]


table(d0$Study,d0$MSI,useNA="i") -> x
x
round(x/rowSums(x)*100,0)

table(d0$`Study ID`,d0$MSI,useNA="i") -> x
x
round(x/rowSums(x)*100,0)

# plots with the differences between cohorts

cf_sr <- ggplot() + coord_flip() + scale_x_discrete(limits = rev)

# to have all the graphs with the same width, they are plotted WITHOUT legends
# legends can be added later by changing show.legend = T


pdf("~/Desktop/panel1.pdf",6,1.5,pointsize = 12)
cf_sr + geom_bar(aes(Study,fill=Sex),position=position_fill(reverse=T),data=d0,color=1,lwd=.2,show.legend = F) + 
  ylab(NULL) + xlab(NULL) + 
  scale_fill_brewer(palette = "Greens",na.value="grey70") +
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))
dev.off()

pdf("~/Desktop/panel2.pdf",6,1.5,pointsize = 20)
cf_sr + geom_boxplot(aes(Study,`Age at Diagnosis`,fill=Study),data=d0,show.legend = F) + xlab(NULL) +
  scale_fill_brewer(palette = "Blues",na.value="grey70")
dev.off()

pdf("~/Desktop/panel3.pdf",6,1.5,pointsize = 12)
cf_sr + geom_bar(aes(Study,fill=Stage),position=position_fill(reverse=T),data=d0,color=1,lwd=.2,show.legend = F) + 
  ylab(NULL) + xlab(NULL) +
  scale_fill_brewer(palette = "Reds",na.value="grey70") +
  scale_y_continuous(labels=paste0(seq(0,100,25),"%")) 
dev.off()

pdf("~/Desktop/panel4.pdf",6,1.5,pointsize = 12)
cf_sr + geom_bar(aes(Study,fill=MSI),position=position_fill(reverse=T),data=d0,color=1,lwd=.2, show.legend = F) + 
  ylab(NULL) + xlab(NULL) +
  scale_fill_brewer(palette="Blues",na.value="grey20",direction=-1) +
  scale_y_continuous(labels=paste0(seq(0,100,25),"%"))
dev.off()


# Co-mutations with NF1 -----------

table(mss$KRAS != "WT",mss$NF1 != "WT") %>% fisher.test()
table(msi$KRAS != "WT",msi$NF1 != "WT") %>% fisher.test()

cbind(
G12=table(mss[grep("G12",KRAS),NF1!="WT"]),
G13=table(mss[grep("G13",KRAS),NF1!="WT"]),
Q61=table(mss[grep("Q61",KRAS),NF1!="WT"])[c("FALSE","TRUE")],
A146=table(mss[grep("A146",KRAS),NF1!="WT"])
) -> x
x[is.na(x)] <- 0
fisher.test(x)

cbind(
  G12=table(msi[grep("G12",KRAS),NF1!="WT"]),
  G13=table(msi[grep("G13",KRAS),NF1!="WT"]),
  Q61=table(msi[grep("Q61",KRAS),NF1!="WT"])[c("FALSE","TRUE")],
  A146=table(msi[grep("A146",KRAS),NF1!="WT"])
) -> x
fisher.test(x)

round(x[2,] / (x[1,] + x[2,]) * 100,0)


# KRAS / NRAS / HRAS

x <- melt(d0,id.vars = c("Study","Study ID","MSI"),measure.vars = c("KRAS","HRAS","NRAS")) %>% data.table
x[value != "WT",value := "MUT"]
x$MSI <- factor(x$MSI,c("MSS","MSI/HYPER"))
x$variable <- factor(x$variable,c("KRAS","NRAS","HRAS"))

table(x$variable,x$value,x$MSI)


ggplot(x) + aes(variable,fill=value) + 
  geom_bar(position=position_fill(reverse=T),color="black",lwd=.25,show.legend = F) + 
  facet_grid(rows=vars(MSI),cols=vars(Study)) +
  scale_y_continuous(breaks=seq(0,1,.25),labels = sprintf("%i%%",seq(0,100,25))) +
  xlab(NULL) + ylab("Percentage of mutated samples") +
  geom_hline(yintercept = c(.25,.5,.75),lwd=.2,color="grey25") +
  theme(axis.text.x.bottom = element_text(angle=90,vjust = .5,hjust=1)) +
  scale_fill_manual(values=c("darkorange","grey95"))


y <- xtabs(~ variable + value + `Study ID` + MSI,x)
y[,,1,1] %>% cp1


cp1 <- function(x) round(100*x/rowSums(x),1)

table(msi$Study,msi$KRAS != "WT")
table(msi$Study,msi$KRAS != "WT") %>% cp1

table(msi$Study,msi$NRAS != "WT")
table(msi$Study,msi$NRAS != "WT") %>% cp1

table(msi$Study,msi$HRAS != "WT")
table(msi$Study,msi$HRAS != "WT") %>% cp1



xtabs(~ (KRAS != "WT") + Study,d0,subset=MSI=="MSS")
xtabs(~ (KRAS != "WT") + Study,d0,subset=MSI!="MSS")

x <- table(d0$Study,d0$KRAS != "WT",d0$MSI,useNA="i")

x[,,1] / rowSums(x[,,1])



x <- table(d0$Study,d0$NRAS != "WT",d0$MSI,useNA="i")
x
x[,,2]/rowSums(x[,,2]) * 100

x <- table(d0$Study,d0$HRAS != "WT",d0$MSI,useNA="i")
x
x[,,2]/rowSums(x[,,2]) * 100

ggplot(d0) + aes(`TMB (nonsynonymous)`,fill = Study) + 
  geom_histogram(alpha=.8,breaks=c(0:10,12,14,100),closed="left") +
  facet_wrap(vars(Study)) + xlab("TMB (nonsynonymous)") 

ggplot(d0) + aes(`TMB (nonsynonymous)` %>% cut(.,c(seq(0,24,2),1000),c(seq(2,24,2),">24"))) +
  geom_bar(aes(fill=Study)) + facet_wrap(vars(Study)) +
  xlab("TMB (nonsynonymous)") + 
  geom_vline(data = filter(d0$Study == "TCGA"),xintercept = 3)
  
table(cut(d0$`TMB (nonsynonymous)`,c(0:25,1000)),d0$Study)

ggplot(d0) + aes(`TMB (nonsynonymous)` %>% log10,fill = `Study ID`) + 
  geom_density(col="black",alpha=.2) +
 xlab("TMB (log10)") +
  geom_vline(xintercept = log10(c(14,25)),lty=2) 

table(d0$`Study ID`,d0$`TMB (nonsynonymous)` > 10)
table(d0$`Study ID`,d0$`TMB (nonsynonymous)` > 14)


# analyze mutations
# take into account that some samples have more than one mutation

table(d0$KRAS,d0$`Study ID`)



View(d0[,c("Study ID","Sample Type")])


View(mutations)

dups <- mutations[duplicated(SAMPLE_ID),SAMPLE_ID]
mutations[SAMPLE_ID == dups[2]]

table(mutations$KRAS != "WT",mutations$BRAF != "WT")

mutations[KRAS != "WT" & BRAF != "WT"]


table(mutations$STUDY_ID,mutations$KRAS)
table(clinical$`Study ID`)

table(mutations$KRAS) %>% sort

all(mutations$SAMPLE_ID %in% clinical$`Sample ID`)
all(clinical$`Sample ID` %in% mutations$SAMPLE_ID)

x <- merge(clinical,mutations,by.x="Sample ID",by.y="SAMPLE_ID")
dim(x)

any(duplicated(mutations$SAMPLE_ID))




x <- merge(coad_read,kras[,c("Sample ID","Protein Change")],all.x=T)

x[is.na(`Protein Change`),`Protein Change` := "WT"]

factor(x$`Protein Change`) -> x$`Protein Change`
gsub(".$","",x$`Protein Change`,perl = T) %>% 
  factor(.,levels=c("W","G12","G13","Q61","A146")) -> x$KRAS

x[is.na(KRAS),KRAS := "OTHER"]
levels(x$KRAS)[1] <- "WT"
# check survival

x$Stage <- x$`Neoplasm Disease Stage American Joint Committee on Cancer Code`
x$Stage <- gsub("[^IV]","",x$Stage)
x$Stage <- ordered(x$Stage,c("I","II","III","IV"))


s1 <- Surv(x$`Overall Survival (Months)`/12,x$`Overall Survival Status` == "1:DECEASED")

survfit(s1 ~ Stage,x) %>% ggsurvplot()
survfit(s1 ~ Stage,x,subset=x$`TMB (nonsynonymous)` < 10) %>% ggsurvplot()
survfit(s1 ~ Stage,x,subset=x$`TMB (nonsynonymous)` > 10) %>% ggsurvplot()

survfit(s1 ~ `TMB (nonsynonymous)` > 10,x,subset=Stage >= "III") %>% ggsurvplot()

# Differences in survival depending on the mutation?

coxph(s1 ~ KRAS,x)
survfit(s1 ~ KRAS,x) %>% ggsurvplot()



# Differences in age

aov(`Diagnosis Age` ~ KRAS,x,subset=`TMB (nonsynonymous)` < 10) %>% TukeyHSD()
aov(`Diagnosis Age` ~ KRAS,x) %>% anova()

ggplot(x) + aes(KRAS,Stage) + geom_bar()



survdiff(s1 ~ mut2,data = x,subset=tmb < 10 & mut2 %in% c("G12","G13"))

table(x$mut2,x$stage) %>% plot

coxph(s1 ~ I(mut2 == "G12"),x,subset=stage=="I") 
coxph(s1 ~ I(mut2 == "G12"),x,subset=stage=="II") 
coxph(s1 ~ I(mut2 == "G12"),x,subset=stage=="III") 
coxph(s1 ~ I(mut2 == "G12"),x,subset=stage=="IV") 


x <- merge(coad_read,kras[,c("Sample ID","Protein Change")],by="Sample ID",sort=F,all.x=T)
names(x)[62] <- "KRAS"
x <- merge(x,nf1[,c("Sample ID","Protein Change")],by="Sample ID",all.x=T)
names(x)[63] <- "NF1"

x$Hypermutant <- x$`TMB (nonsynonymous)` > 10

table(x$KRAS, !is.na(x$NF1), x$Hypermutant)

table(x$KRAS,x$Hypermutant) %>% as.data.table -> foo
foo[order(N)]

xtabs(~ KRAS + !is.na(NF1) + Hypermutant,x)
xtabs(~ KRAS + !is.na(NF1),x,subset=Hypermutant)

subset(x,!Hypermutant)$`Sample ID` %>% duplicated()

x[duplicated(`Sample ID`),`Sample ID`] %>% unique -> duplicated_samples

x[`Sample ID` %in% duplicated_samples]

table(x$KRAS,!is.na(x$NF1))

plot(a,b)

boxplot(x$`TMB (nonsynonymous)` ~ b > 10)

summary(lm(x$`TMB (nonsynonymous)` ~ a + b))


View(kras)
(table(kras$`Protein Change`)) %>% sort(.,decreasing = T)
(table(kras$`Protein Change`)/526*100) %>% sort(.,decreasing = T)


hist(kras$`TMB (nonsynonymous)` %>% log10)

ggplot(kras) + aes(`Protein Change`,`Diagnosis Age`) + geom_boxplot()
table(kras$`Protein Change`,kras$`TMB (nonsynonymous)` > 10)
