library(cBioPortalData)
library(tidyr)
library(data.table)
library(biomaRt)

#try(setwd("~/Documents/Work/KRAS/"),silent=T)
#try(setwd("/imppc/labs/cge/salonso/KRAS/"),silent=T)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl) %>% as.data.table()
filters <- listFilters(ensembl) %>% as.data.table()

attributes[grep("length",name)]
filters[grep("bioty",name)]

att1 <-  c("ensembl_gene_id",
           "ensembl_transcript_id",
           "entrezgene_id",
           "transcript_is_canonical",
           "hgnc_symbol",
           "gene_biotype")

attributes[match(att1,name)]

genes <- getBM(attributes = att1,
               filters = "biotype",
               values = list(c("protein_coding")),
               mart=ensembl) %>% as.data.table

dim(genes)


gene.length <- getBM(attributes = c("ensembl_transcript_id",
                                    "cds_length"),
                     filters = "biotype",
                     values = list(c("protein_coding")),
                     mart=ensembl) %>% as.data.table
  
dim(gene.length)


genes <- merge(genes,gene.length) # all.transcripts
rm(gene.length)

genes[hgnc_symbol=="BRAF"] # test 

cbio <- cBioPortal()

studies <- getStudies(cbio) %>% data.table
studies[grep("coad",studyId)]

study_ids <- c("coadread_dfci_2016","coadread_tcga_pan_can_atlas_2018","crc_msk_2017","msk_chord_2024")

clinical_data <- lapply(study_ids,function(x) {
  
  clinicalData(cbio,studyId = x) %>% as.data.table
  
})

lapply(clinical_data,dim)

x <- data.table(genePanels(cbio))[grep("^IMPACT[0-9]",genePanelId),genePanelId]
impact <- lapply(x,function(panel) getGenePanel(cbio,panel)$hugoGeneSymbol)
names(impact) <- x                                       

wes <- genes[hgnc_symbol!="",unique(hgnc_symbol)]

mutational_data <- lapply(study_ids,function(x) {
  
  genes.to.download <- wes
  if(x=="msk_chord_2024") genes.to.download <- impact$IMPACT505
  
  getDataByGenes(cbio,
                 studyId = x,
                 genes = genes.to.download,
                 by = "hugoGeneSymbol",
                 molecularProfileIds = paste0(x,"_mutations"))[[1]] %>% as.data.table
  
})

lapply(mutational_data,dim)

names(clinical_data) <- study_ids
names(mutational_data) <- study_ids

save(genes,impact,wes,clinical_data,mutational_data,file="data/raw.data.RData",compress=TRUE)

