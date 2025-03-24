library(data.table)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(gridExtra)


# useful functions to expand and contract terms

join1 <- function(x) paste(x,collapse="|")
split1 <- function(x) strsplit(x,split="\\|")

# ggplot theme

palette <- list()

palette$msi <- c(MSS="#BCBDDC",MSI="#6A51A3","MSI/HYPER"="#6A51A3")
palette$mut <- c(WT="#FEE6CE",MUT="#F16913")
palette$study <- c("#C19F08","#369E4B","#357EB9","#FA792F","#ce3206")
names(palette$study) <- c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis")
palette$study2 <- brewer.pal(9,"Greens")[c(3,5:9)]
names(palette$study2) <- c("ALL","DFCI Primary","TCGA Primary","MSKCC Primary","MSKCC Metastasis")

palette$stage <- brewer.pal(9,"Blues")[c(2,4,6,8)]
names(palette$stage) <- c("I","II","III","IV")
palette$sample.type <- c("Primary"="#8EBEDA","Metastasis"="#2478B1")

mytheme <- theme(panel.background = element_rect(fill="grey95"),
                 strip.background = element_rect(fill="grey90",colour="grey30"),
                 panel.grid = element_line(color="white"),
                 panel.border = element_rect(color="grey30",fill=NA,linewidth=rel(1)),
) 



if(!file.exists("data/Thesaurus.txt")) {
  
  curl::curl_download("https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus.FLAT.zip",destfile = "data/Thesaurus.zip")
  
}


thesaurus <- fread("data/Thesaurus.zip",quote="")
names(thesaurus) <- c("code","concept IRI","parents","synonyms","definition","display name","concept status","semantic type","concept in subset")
thesaurus$synonyms %>% split1 %>% sapply(FUN=function(x) x[1]) -> thesaurus$name # the preferred name
setkey(thesaurus,"code")

# remove retired concepts

thesaurus <- thesaurus[!grepl("Retired",`concept status`)]

# calculate the number of children per term 

x <- thesaurus$parents %>% split1 %>% unlist %>% table %>% data.table()
names(x) <- c("code","children")

thesaurus <- merge(thesaurus,x,by="code",all.x=T)
thesaurus[is.na(children),children:=0]

isDrug <- grep("Pharmacolo",thesaurus$`semantic type`)

# reconstruct the thesaurus tree

library(igraph)

relations <- thesaurus[,list(parents=strsplit(parents,split="\\|")[[1]]),by=code]
thesaurusGraph <- graph_from_edgelist(as.matrix(relations[,2:1]))

V(thesaurusGraph)$label <- thesaurus[match(V(thesaurusGraph)$name,code),name]

findChildren <- function(term) {
  x <- subcomponent(thesaurusGraph,term,"out") %>% names
  x <- setdiff(x,term)
  return(thesaurus[match(x,code),list(code,name,children)])
}

findParents <- function(term) {
  x <- subcomponent(thesaurusGraph,term,"in") %>% names
  x <- setdiff(x,term)
  return(thesaurus[match(x,code),list(code,name,children)])
}


# Create a data.table with KRAS/BRAF related drugs ----

x <- thesaurus[c(findChildren("C1902")$code,"C62508")][children==0]
x[,Group:="KRAS"]
drugs <- x

# with braf inhibitors
x <- thesaurus[findChildren("C2336")$code]
x <- rbind(x,thesaurus[isDrug][grep("raf kinase|kinases raf",paste(definition,synonyms),ignore.case = T)]) %>% unique
x <- rbind(x,thesaurus[isDrug][grep("pan-raf",paste(name,synonyms),ignore.case = T)]) %>% unique
toadd <- c("C126653","C121958") # radioconjugates
x <- rbind(x,thesaurus[toadd]) %>% unique

x[,Group:="BRAF"]

drugs <- rbind(drugs,x[children==0])

# SOS1 and SHP2

x <- thesaurus[isDrug][grepl("sos[12]",paste(synonyms,definition),ignore.case = T)]
x[,Group:="SOS1"]

drugs <- rbind(drugs,x)

x <- thesaurus[findChildren("C185612")$code]
x[,Group:="SHP2"]

drugs <- rbind(drugs,x)

# RAS TCRs

x <- thesaurus[isDrug][grepl("TCR",synonyms) & grepl("k[-]*ras",synonyms,ignore.case = T)]
x[,Group:="KRAS"]
x[,Type:="TCR"]
drugs <- rbind(drugs,x,fill=T)

# RAS vaccines

x <- thesaurus[grep("K[-]*RAS.+vaccine|KRAS-EphA-2-CAR-DC|ELI-002|GRT-C903|KRAS-EphA-2-CAR-DC",synonyms,ignore.case = T)]
x[,Group:="KRAS"]
x[,Type:="Vaccine"]
drugs <- rbind(drugs,x)

# pan-ras or pan-kras perhaps not included?

x <- thesaurus[isDrug][grep("pan-ras|pan-kras",paste(definition,synonyms),ignore.case = T)]
x$code %in% drugs$code

table(drugs$Group)
drugs[duplicated(code)]

# classification by type and target molecule

drugs$Target <- NULL

#G12C
drugs[match(findChildren("C201853")$code,code),Target:="G12C"]
#BI 3706674 is a pan-KRAS
drugs[grep("BI 3706674",name),Target:="pan-KRAS"]
drugs[grep("G12C",synonyms)]
drugs[grep("fulzerasib|divarasib",synonyms,ignore.case = T),Target:="G12C"]

#G12D
drugs[match(findChildren("C207236")$code,code),Target:="G12D"]
drugs[grep("G12D",synonyms),Target:="G12D"]

#G12V
drugs[grep("G12V",synonyms),Target:="G12V"]

#pan-KRAS (or multi)
drugs[grep("C185876|C188048|C178264|C71146|C206998",code),Target:="pan-KRAS"]
drugs[grep("C209714|C204884|C202009|C204252|C162269|C165638|C2067|C179231|C165514|C200465|C162186",code),Target:="multi-KRAS"]


#SOS1/SHP2
drugs[grep("SOS1|SHP2",Group),Target:="SOS1|SHP2"]

#BRAF
drugs[grep("BRAF",Group),Target:="BRAF"]

#OTHER
drugs[grep("Diazepinomicin",name),`:=`(Type="Inhibitor",Target="pan-KRAS")]

#Classify drug type
drugs[grep("inhibitor",name,ignore.case = T),Type:="Inhibitor"]
drugs[grep("fenib|protafib|rasib|tinib|Paluratide|Rigosertib",name),Type:="Inhibitor"]
drugs[grep("degrader",synonyms,ignore.case = T),Type:="Degrader"]

drugs[grep("Glue",name),Type:="Molecular Glue"]
drugs[grep("Antisense",name),Type:="Antisense"]

drugs[grep("Urea",name),`:=`(synonyms=paste0(synonyms,"|heterocyclic urea"),Type="Inhibitor")]

# pan-ras wrongly identified?

drugs[grep("pan-",paste(name,synonyms,definition),ignore.case=T),list(name,Type,Target)]
drugs[grep("G12C",paste(name,synonyms)),list(name,Type,Target)]
drugs[grep("wild-type",paste(name,synonyms,definition),ignore.case=T),list(name,Type,Target)]
drugs[grep("3706674",name)]

# remove some salt forms of drugs

remove <- "Divarasib Adipate|Rigosertib Sodium|Sorafenib Tosylate|Dabrafenib Mesylate|Lifirafenib Maleate|Regorafenib Anhydrous"

drugs <- drugs[grep(remove,name,invert = T)]


# correct paluratide = LUNA18

drugs[grep("LUNA18",synonyms,ignore.case = T),synonyms:=paste0(synonyms,"|paluratide|PALURATIDE")]
drugs <- drugs[name!="Paluratide"] 

drugs[,From:="NCI Thesaurus"]


# additional manually added drugs

additional.drugs <- fread("data/drug.revision.txt")

for(i in 1:nrow(additional.drugs)) {
  
  hits <- any(grepl(additional.drugs$synonyms[i],drugs$synonyms))
  if(hits) print(additional.drugs[i])
  
}


additional.drugs[,`:=`(From="Manually Added",children=0)]

drugs <- rbind(drugs,additional.drugs,fill=T)

# short names for figures ----

drugs[,short.Name:=name]
drugs[,short.Name := gsub("^.*(Inhibitor|Degrader|Glue|Oligonucleotide|Fluorine|Carbon) ","",short.Name)]
drugs[,short.Name := gsub("mTCR-transduced Autologous Peripheral Blood Lymphocytes","mTCR",short.Name)]
drugs[,short.Name := gsub("^.* Urea","Policyclic Urea",short.Name)]
drugs[,short.Name := gsub("K-RAS","KRAS",short.Name)]

drugs[grep(" NT-112",short.Name),short.Name := "NT-112"]
drugs[grep(" AFNT-211",short.Name),short.Name := "AFNT-211"]
drugs[grep("C204252",code),short.Name := "Anti-KRAS Mutant mTCR"]
drugs[Type=="Vaccine",short.Name := gsub("^.* Vaccine ","",short.Name)]

drugs[name=="Autologous mDC3/8-KRAS Vaccine",short.Name:="Autologous mutant-KRAS Vaccine"]
drugs[name=="Pooled Mutant KRAS-Targeted Long Peptide Vaccine",short.Name:="Pooled mutant-KRAS Vaccine"]
drugs[name=="K-RAS Protein Vaccine",short.Name:="Mutant-KRAS Vaccine"]
drugs[order(nchar(short.Name)),list(name,short.Name)]

# end creating the drugs data.table -------

drugs[,table(From,Group)]

drugs[is.na(code),code:="----"]

x <- drugs[Group!="BRAF",list(code,short.Name,Type,Target)]
x[Target=="G12C",Target:="KRAS G12C"]
x[Target=="G12D",Target:="KRAS G12D"]
x[Target=="G12V",Target:="KRAS G12V"]

x[,Target := factor(Target,c("KRAS G12C","KRAS G12D","KRAS G12V","multi-KRAS","pan-KRAS","SOS1|SHP2","PBR"))]
x <- x[order(Target,Type,short.Name)]

colnames(x) <- c("NCI Thesaurus\nCode","Short Name","Type","Target")

png("output/table1.png",1400,1200)
g1 <- tableGrob(x[1:42],rows=NULL,theme=ttheme_default(base_size=20))
g2 <- tableGrob(x[43:84],rows=NULL,theme=ttheme_default(base_size=20))
grid.arrange(g1,g2,layout_matrix=matrix(1:2,ncol=2))
dev.off()

# Clinical trials for colorectal cancer or solid tumors ----

ctrials <- fread("data/ctg-studies.csv")
setkey(ctrials,"NCT Number")

# correct some misspellings
ctrials[grep("LY4066434.",Interventions),Interventions:=gsub("LY4066434\\.","LY4066434",Interventions)] 
ctrials[grep("5Fluoro",Interventions),Interventions:=gsub("5Fluoro","5-Fluoro",Interventions)] 
ctrials[grep("5FLUORO",Interventions),Interventions:=gsub("5FLUORO","5-FLUORO",Interventions)] 

dim(ctrials)

# search in ctrials

crc <- ctrials[grepl("colon|rectal|gastrointestinal|solid|crc",Conditions,ignore.case = T) & `Study Type`=="INTERVENTIONAL"]
crc <- crc[grep("TREATMENT",`Study Design`)]
crc <- rbind(crc,ctrials[grep("NCT06385925|NCT04330664|NCT03785249",`NCT Number`)])

drug.columns <- c("Drug.Code","Drug.Name","Drug.Type","Drug.Target","Drug.Group")
crc[,(drug.columns) := ""]


for(i in 1:nrow(drugs)) {
  hits <- grep(drugs$synonyms[i],crc$Interventions,ignore.case = T)
  crc[hits,Drug.Code:=paste0(Drug.Code,drugs$code[i],"|")]
  crc[hits,Drug.Name:=paste0(Drug.Name,drugs$short.Name[i],"|")]
  crc[hits,Drug.Type:=paste0(Drug.Type,drugs$Type[i],"|")]
  crc[hits,Drug.Target:=paste0(Drug.Target,drugs$Target[i],"|")]
  crc[hits,Drug.Group:=paste0(Drug.Group,drugs$Group[i],"|")]
  }

# remove the last |

crc[,(drug.columns) := lapply(.SD,function(x) gsub("\\|$","",x)),.SDcols=drug.columns]

dim(crc)

kras.trials <- crc[grep("KRAS|SHP2|SOS1",Drug.Group)]
kras.trials[,(drug.columns) := ""]

drugs2 <- drugs[Group!="BRAF"]

for(i in 1:nrow(drugs2)) {
  hits <- grep(drugs2$synonyms[i],kras.trials$Interventions,ignore.case = T)
  kras.trials[hits,Drug.Code:=paste0(Drug.Code,drugs2$code[i],"|")]
  kras.trials[hits,Drug.Name:=paste0(Drug.Name,drugs2$short.Name[i],"|")]
  kras.trials[hits,Drug.Type:=paste0(Drug.Type,drugs2$Type[i],"|")]
  kras.trials[hits,Drug.Target:=paste0(Drug.Target,drugs2$Target[i],"|")]
  kras.trials[hits,Drug.Group:=paste0(Drug.Group,drugs2$Group[i],"|")]
}

kras.trials[,(drug.columns) := lapply(.SD,function(x) gsub("\\|$","",x)),.SDcols=drug.columns]

kras.trials[Drug.Target!=""]


# trials in review 2024
inReview <- c("NCT05485974","NCT05462717","NCT06244771","NCT04121286","NCT06385925","NCT05737706","NCT06364696","NCT06403735","NCT06040541","NCT06412198","NCT05194995","NCT04330664","NCT05288205","NCT06039384","NCT05198934","NCT04956640","NCT05578092","NCT04699188","NCT04449874","NCT05358249","NCT06026410","NCT06252649","NCT06586515","NCT04793958","NCT03785249","NCT05722327","NCT06599502","NCT06078800","NCT06607185","NCT06447662","NCT06585488","NCT06445062","NCT05379985","NCT05786924","NCT06194877","NCT05200442","NCT06270082","NCT05163028","NCT04117087","NCT04853017","NCT06411691","NCT05726864","NCT06105021","NCT06253520","NCT06487377","NCT06218914")
inReview <- ctrials[inReview]

setdiff(inReview$`NCT Number`,crc[Drug.Type!="",`NCT Number`]) %>% ctrials[.] -> missed
dim(missed)

# check missed trials

ctrials[grepl("k[-]*ras",Interventions,ignore.case = T) & ! `NCT Number` %in% kras.trials$`NCT Number`]










# review trials with more than one identified drug

kras.trials[grep("\\|",Drug.Name),list(`NCT Number`,Interventions,Drug.Name)]



# save for manual revision 

fwrite(kras.trials,"data/kras.studies.csv",sep=",")

# fread the revised version

try(kras.trials <- fread("data/kras.studies.revised.csv"),silent=T)

# add bibligraphic information (takes several min)

# library(easyPubMed)
# 
# sapply(kras.trials$`NCT Number`,function(x) {
#   get_pubmed_ids(x) %>% fetch_all_pubmed_ids() -> ids
#   return(paste(sort(ids),collapse="|"))
# }) -> kras.trials$pmids
# 
# kras.trials[pmids!="",`Study Type`]
# 
# 
# 
# 
# # try to identify acronyms
# 
# acronym <- "[A-Z]{3,}[^ ]*"
# kras.trials[grep(acronym,`Study Title`) ,`Study Title`]
# 
# kras.trials[,str_extract(`Study Title`,acronym) %>% paste(collapse=", "),`NCT Number`][1:100]


# plots of the Clinical Trials -----

# adjust dates for ggplot representation

kras.trials[nchar(`Start Date`) > 8,start:=as.Date(`Start Date`)]
kras.trials[nchar(`Start Date`) == 7,start:=as.Date(paste0(`Start Date`,"-01"))]

kras.trials[nchar(`Completion Date`) > 8,completion:=as.Date(`Completion Date`)]
kras.trials[nchar(`Completion Date`) == 7,completion:=as.Date(paste0(`Completion Date`,"-01"))]
kras.trials[is.na(completion)]

my.palette <- list()
my.palette$kras.type <- c("#3588C7","#1A4EA0","#8FBA6F","#4F8B7A","#405D02","#F8CB1C","#C65684","#EF0060")

kras.trials$Drug.Type %>% split1 %>% lapply(FUN=function(x) x %>% unique %>% sort %>% join1) %>% unlist -> kras.trials$Drug.Type

library(ggh4x)
library(ggnewscale)
library(ggtext)
library(gridExtra)


my.palette$targets <- c("G12C"="#03396c",
                        "G12D"="#03396c",
                        "G12V"="#03396c",
                        "multi-KRAS"="#ff4d00",
                        "pan-KRAS"="#aa2000",
                        "SOS1|SHP2"="#aa2000",
                        "BRAF"="#454545",
                        "RAF"="#454545")

my.palette$targets[drugs$Target] -> drugs$color


plot.trials <- list(aes(x=start,y=`NCT Number`),
                    scale_x_date(date_breaks = "year", date_labels = "%Y",limits = as.Date(c("2000-01-01","2040-12-31"))),
                    geom_vline(xintercept = as.Date("2024-11-20"),color="grey",lwd=3),
                    geom_segment(aes(x=start,xend=completion),lwd=1.5,show.legend = F),
                    scale_color_manual("Drug Type",values=my.palette$kras.type),
                    new_scale_color(),
                    geom_point(size=3,pch=15,aes(color=`Study Status`),show.legend = F),
                    scale_color_manual("Status",
                                       values=c(NOT_YET_RECRUITING="#94C58C",
                                                RECRUITING="#1A8828",
                                                ENROLLING_BY_INVITATION="#0A6921",
                                                ACTIVE_NOT_RECRUITING="#3170DE",
                                                COMPLETED="#092B9C",
                                                WITHDRAWN="#de894a",
                                                SUSPENDED="#cb5323",
                                                TERMINATED="#B22222",
                                                UNKNOWN="#888888")),
                    
                    geom_point(aes(x=as.Date(`Last Update Posted`)),pch=21,lwd=0.2,fill="white"),
                    
                    scale_y_discrete(limits=rev),
                    geom_richtext(aes(label=labels),hjust=1,nudge_x = -365,size=3,fill = NA,label.colour=NA),
                    
                    
                    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)),
                    ylab(""),
                    xlab(""),
                    mytheme,
                    theme(panel.background = element_rect(fill="white")))




kras.trials$Drug.Name %>% strsplit(split="\\|") %>% 
  sapply(FUN=function(x) {
    
    color <- drugs[match(x,short.Name),color] 
    sprintf("<span style='color:%s'>%s</span>",color,x) %>% paste(.,collapse=" | ")
    
  }) -> kras.trials$labels



g1 <- ggplot(kras.trials[`Study Status`=="NOT_YET_RECRUITING"]) + 
  plot.trials +
  ggtitle("Not yet recruiting")

g2 <- ggplot(kras.trials[`Study Status`=="RECRUITING"]) + 
  plot.trials +
  ggtitle("Recruiting")

g3 <- ggplot(kras.trials[`Study Status`=="ENROLLING_BY_INVITATION"]) + 
  plot.trials +
  ggtitle("Enrolling by invitation")

g4 <- ggplot(kras.trials[`Study Status`=="ACTIVE_NOT_RECRUITING"]) + 
  plot.trials +
  ggtitle("Active, not recruiting")

g5 <- ggplot(kras.trials[`Study Status`=="COMPLETED"]) + 
  plot.trials +
  ggtitle("Completed")


library(gridExtra)
grid.arrange(g1,g2,g3,g4,g5,layout_matrix=matrix(c(1,2,2,3,4,5),ncol=2),heights=c(1.5,5,2))




