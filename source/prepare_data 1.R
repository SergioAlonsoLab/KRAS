# load required libraries ----

library(data.table)
library(tidyr)
library(ggplot2)

# Set working directory -----

#try(setwd("/imppc/labs/cge/salonso/KRAS/"))
#try(setwd("~/Documents/Work/KRAS/"))

# Load datasets ---

load("data/raw.data.RData",verbose = T)

fwrite(genes[transcript_is_canonical==1],file="data/all.genes.length.csv")

studies <- names(clinical_data)

# select only CRC samples (MSKCC is a pancancer study)

lapply(clinical_data, function(dt) {
  
  dt[CANCER_TYPE=="Colorectal Cancer"]
  
}) -> clinical_data

names(clinical_data) <- studies

# repeat mutation_count

for(i in studies) {
  mcount <- mutational_data[[i]][,list(MUT_COUNT=.N),by=sampleId]
  clinical_data[[i]] <- merge(clinical_data[[i]],mcount,all.x=T)
  clinical_data[[i]][,STUDY:=i]
}


# create a unique clinical table ----

rbindlist(clinical_data[1:3],fill=T) -> clinical

clinical[,SAMPLE_COUNT := as.numeric(SAMPLE_COUNT)]
clinical[,SAMPLE_COUNT2 := .N,by=patientId]
clinical[SAMPLE_COUNT!=SAMPLE_COUNT2]

# remove samples without mutational information

clinical[,table(STUDY,is.na(MUT_COUNT))]
clinical[,table(STUDY,is.na(MUTATION_COUNT))]

clinical <- clinical[!is.na(MUTATION_COUNT)]

# number of samples per patient 

clinical[,SAMPLES_PATIENT:=.N,by=patientId]

clinical[,SAMPLE_COUNT:=as.numeric(SAMPLE_COUNT)]
clinical[SAMPLE_COUNT!=SAMPLES_PATIENT]

# recalculate TMB

clinical[is.na(GENE_PANEL),GENE_PANEL:="WES"]

clinical[,table(GENE_PANEL)]

selected.genes <- impact$IMPACT341

mutations <- rbindlist(mutational_data[1:3],fill=TRUE)[,list(sampleId,hugoGeneSymbol,proteinChange)]

# check if all samples form clinical are in the mutation table

all(clinical$sampleId %in% mutations$sampleId)

# select just genes that have been analysed in all datasets
mutations <- mutations[hugoGeneSymbol %in% selected.genes]
setdiff(selected.genes,mutations$hugoGeneSymbol) # PMAIP1 is never mutated
mutations <- rbind(mutations,data.table(sampleId="foo",hugoGeneSymbol="PMAIP1",proteinChange="WT"))

# for merging purposes, add all samples to the list of mutations

setdiff(clinical$sampleId,mutations$sampleId)
#mutations <- rbind(mutations,data.table(sampleId=clinical$sampleId,hugoGeneSymbol="no_gene",proteinChange="no_change"))

# convert in atable
mutations <- dcast(mutations,sampleId ~ hugoGeneSymbol,value.var = "proteinChange",fun.aggregate = function(x) paste(x,collapse="|"))
mutations[mutations==""] <- "WT"

setdiff(clinical$sampleId,mutations$sampleId) # some samples do not harbour mutations in those genes

# create bowel DT, the table used for most figures an analysis

bowel <- merge(clinical,mutations,all.x=T)

for(i in selected.genes) {
  
  bowel[[i]][is.na(bowel[[i]])] <- "WT"
  
}



# Annotate KRAS, NRAS and HRAS codons ---

bowel$KRAS.CODON <- "Other"
bowel$NRAS.CODON <- "Other"
bowel$HRAS.CODON <- "Other"

for(i in c("WT","G12","G13","Q61","A146","K117","\\|")) {
  j <- i
  if(i=="\\|") j <- "Multiple"
  bowel[grep(i,KRAS),KRAS.CODON:=j]
  bowel[grep(i,NRAS),NRAS.CODON:=j]
  bowel[grep(i,HRAS),HRAS.CODON:=j]
}

rm(i,j)


# Age

bowel[!is.na(AGE_AT_DIAGNOSIS),AGE:=AGE_AT_DIAGNOSIS]
bowel[is.na(AGE)]
bowel[,AGE:=as.numeric(AGE)]

# Annotate sample type ----
# DFCI samples are primary according to the paper, but they are not annotated

bowel[is.na(SAMPLE_TYPE),table(STUDY)]
bowel[is.na(SAMPLE_TYPE),SAMPLE_TYPE:="Primary"]

# Rename STUDY_ID -----

bowel[,STUDY_ID:=STUDY]
bowel[grep("msk",STUDY),STUDY:="MSKCC"]
bowel[grep("tcga",STUDY),STUDY:="TCGA"]
bowel[grep("dfci",STUDY),STUDY:="DFCI"]
bowel[,STUDY:=paste(STUDY,SAMPLE_TYPE)]
bowel[,STUDY:=factor(STUDY,c("DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis"))]

table(bowel$STUDY,bowel$STUDY_ID)

# add information about sex of patients TCGA-5M-AAT5 and TCGA-5M-AATA, previously determined by the 
# methylation profile of chromosome X

bowel[patientId=="TCGA-5M-AAT5",SEX := "Female"]
bowel[patientId=="TCGA-5M-AATA",SEX := "Male"]

# table sex and age ----

bowel[,table(STUDY,SEX)] -> x
x;fisher.test(x)
pairwise.prop.test(x,p.adjust.method = "fdr")
(100*x[,1]/rowSums(x)) %>% round(.,0)

bowel[,sprintf("%1.1f±%1.1f",mean(AGE,na.rm=T),sd(AGE,na.rm=T)),by=STUDY]
aov(AGE ~ STUDY,bowel) %>% summary
aov(AGE ~ STUDY,bowel) %>% TukeyHSD()

# Add and analyze tumour Stage -----

bowel[,table(STAGE_AT_DIAGNOSIS,STUDY)]
bowel[,table(TUMOR_STAGE,STUDY)]
bowel[,table(AJCC_PATHOLOGIC_TUMOR_STAGE,STUDY)]

bowel[!is.na(STAGE_AT_DIAGNOSIS),STAGE:=STAGE_AT_DIAGNOSIS]
bowel[!is.na(TUMOR_STAGE),STAGE:=TUMOR_STAGE]
bowel[!is.na(AJCC_PATHOLOGIC_TUMOR_STAGE),STAGE:=AJCC_PATHOLOGIC_TUMOR_STAGE %>% gsub("[^IV]","",.)]

bowel[,table(STAGE,STUDY,useNA="ifany")]

# Add stage from pathology report from the TCGA ----

casesNoStage <- bowel[is.na(STAGE) & STUDY=="TCGA Primary",patientId]

x <- bowel[is.na(STAGE),list(patientId,AJCC_PATHOLOGIC_TUMOR_STAGE,PATH_T_STAGE,PATH_N_STAGE,PATH_M_STAGE)]
x

openPR <- function(case) {
  sprintf("%s",case) %>% print
  browseURL(sprintf("https://www.cbioportal.org/patient/pathologyReport?studyId=coadread_tcga_pan_can_atlas_2018&caseId=%s",case))
}


# for(i in casesNoStage) {
#   openPR(i)
# }

# T3N2   -> III
# T3N1   -> III
# T2N0   -> I
# T1N0   -> I
# T3N0   -> II
# T3N0   -> II
# T4N2   -> III
# T3N1M1 -> IV
# T3N0   -> II
# T3N0  -> II
# T2N1  -> III
# T3N1  -> III
# T3N0  -> II
# T3N1  -> III

bowel[match(casesNoStage,patientId)]$STAGE <- c("III","III","I","I","II","II","III","IV","II","II","III","III","II","III")

bowel[,table(STAGE,STUDY,useNA="ifany")] -> x
x;chisq.test(x[1:4,])
rm(x)

bowel[,STAGE:=ordered(STAGE)]
bowel[,STAGE2 := STAGE]
bowel[SAMPLE_TYPE=="Metastasis",STAGE2 := "M"]

levels(bowel$STAGE2) <- c("I-II","I-II","III-IV","III-IV","M")


# mutation count ----
bowel[,TMB_NONSYNONYMOUS:=as.numeric(TMB_NONSYNONYMOUS)]
bowel[,HYPER:="ND"]
bowel[STUDY %in% c("DFCI Primary","TCGA Primary"),HYPER := ifelse(TMB_NONSYNONYMOUS>=14,"H","L")]
bowel[STUDY %in% c("MSKCC Primary","MSKCC Metastasis"),HYPER := ifelse(TMB_NONSYNONYMOUS>=25,"H","L")]

bowel[,table(STUDY,HYPER)] %>% {print(.);chisq.test(.)}

# MSI ----

bowel[,MSI:="ND"]
bowel[,MSI_SCORE:=as.numeric(MSI_SCORE)]
bowel[,MSI_SCORE_MANTIS:=as.numeric(MSI_SCORE_MANTIS)]
bowel[,MSI_SENSOR_SCORE:=as.numeric(MSI_SENSOR_SCORE)]

bowel[MSI_SCORE >= 10 | MSI_SENSOR_SCORE > 10,MSI:="MSI"]
bowel[MSI_SCORE < 10 | MSI_SENSOR_SCORE < 10,MSI:="MSS"]

bowel[MSI_SCORE_MANTIS > 0.4,MSI:="MSI"]
bowel[MSI_SCORE_MANTIS < 0.4,MSI:="MSS"]
bowel[,MSI_STATUS]

bowel[MSI_STATUS == "MSS",MSI:="MSS"]
bowel[MSI_STATUS %in% c("MSI","MSI-high"),MSI:="MSI"]

bowel[,table(STUDY,MSI,useNA="ifany")]

# HYPER.MSI ----

bowel[,table(HYPER,MSI,useNA = "ifany")]

bowel[,GROUP:="MSS"]
bowel[MSI=="MSI",GROUP:="MSI/HYPER"]
bowel[HYPER=="H",GROUP:="MSI/HYPER"]
bowel[,GROUP:=factor(GROUP,c("MSS","MSI/HYPER"))]
x <- bowel[,table(GROUP,STUDY)]
x;chisq.test(x)
pairwise.prop.test(t(x),p.adjust.method = "fdr")
rm(x)


fwrite(bowel,"data/bowel.tsv",sep="\t")
saveRDS(bowel,"data/bowel.rds")
