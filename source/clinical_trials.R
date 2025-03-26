library(data.table)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(gridExtra)
library(igraph)
library(stringr)


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


thesaurus <- fread("data/Thesaurus.FLAT.zip",quote="")
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


# Create a data.table with KRAS/BRAF related drugs ------

# Add Ras inhibitors (code C19902)

drugs <- thesaurus[findChildren("C1902")$code][children==0][,Group:="KRAS"] %>% data.table

# Also include Diazopinomicin

x <- thesaurus["C62508"][,Group:="KRAS"] %>% data.table
drugs <- rbind(drugs,x)

# add braf inhibitors
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

x <- thesaurus[isDrug][grepl("TCR",paste(synonyms,definition)) & grepl("k[-]*ras",paste(synonyms,definition),ignore.case = T)]
x[,Group:="KRAS"]
x[,Type:="TCR"]
drugs <- rbind(drugs,x,fill=T)

# RAS vaccines

x <- thesaurus[grepl("vaccine",paste(synonyms,definition),ignore.case = T) &
                 grepl("K[-]*RAS",paste(synonyms,definition),ignore.case = T)]

# VSV-GP154 is not clearly a KRAS vaccine
x <- x[!grepl("VSV-GP154",synonyms)]

x[,Group:="KRAS"]
x[,Type:="Vaccine"]
drugs <- rbind(drugs,x)

# pan-ras or pan-kras not included?

x <- thesaurus[isDrug][grep("pan-ras|pan-kras|pan-k-ras",paste(definition,synonyms),ignore.case = T)]
x$code %in% drugs$code # all included

table(drugs$Group)
drugs[duplicated(code)]

# classification by type and target molecule

# mutations named in synonyms or definition field

drugs[,named_mutations:=str_extract_all(paste(synonyms,definition) %>% toupper,"(G[0-9]+[A-Z]|V600[A-Z]|BRAF|B-RAF|PAN-[K]*RAS|PAN-[B]*RAF)")[[1]] %>% 
        gsub("B-RAF","BRAF",.) %>% unique %>% paste(.,collapse=":"),by=code]



drugs$Target <- NULL

#G12C

drugs[match(findChildren("C201853")$code,code),Target:="G12C"]
drugs[grep("fulzerasib|divarasib|talorasib",synonyms,ignore.case = T),Target:="G12C"]

drugs[grep("G12C",paste(synonyms,definition),ignore.case = T)]

#G12D
drugs[match(findChildren("C207236")$code,code),Target:="G12D"]
drugs[grep("G12D",synonyms),Target:="G12D"]

drugs[grep("G12D",paste(synonyms,definition),ignore.case = T)]

#SOS1/SHP2
drugs[grep("SOS1|SHP2",Group),Target:="SOS1:SHP2"]

#BRAF
drugs[grep("BRAF",Group),Target:="BRAF"]

#OTHER
drugs[grep("Diazepinomicin",name),`:=`(Type="Inhibitor",Target="pan-RAS")]

#Classify drug type
drugs[grep("inhibitor",name,ignore.case = T),Type:="Inhibitor"]
drugs[grep("(nib|fib|sib|tib|ratide)( |$)",name),Type:="Inhibitor"]

drugs[grep("degrader",synonyms,ignore.case = T),Type:="Degrader"]

drugs[grep("Glue",name,ignore.case = T),Type:="Molecular Glue"]
drugs[grep("Antisense",name,ignore.case = T),Type:="Antisense"]

drugs[grep("Urea",name),`:=`(synonyms=paste0(synonyms,"|heterocyclic urea"),Type="Inhibitor")]

# remove some salt forms of drugs
remove <- "Divarasib Adipate|Rigosertib Sodium|Sorafenib Tosylate|Dabrafenib Mesylate|Lifirafenib Maleate|Regorafenib Anhydrous"
drugs <- drugs[grep(remove,name,invert = T)]


# correct paluratide = LUNA18
# https://pubchem.ncbi.nlm.nih.gov/compound/Paluratide

drugs[grep("LUNA18",synonyms,ignore.case = T),synonyms:=paste0(synonyms,"|paluratide|PALURATIDE")]
drugs <- drugs[name!="Paluratide"] 
drugs[,From:="NCI Thesaurus"]


# Correct discrepancies ----

drugs[grep("G12C",paste(name,synonyms)),list(name,Type,Target,named_mutations)] # ok
drugs[grep("G12D",paste(name,synonyms)),list(code,Type,Target,named_mutations)] # not ok, C212035 also targets G12V
drugs[code=="C212035",Target:="multi-KRAS"]

# pan-KRAS vs multi-KRAS


# Vaccines and TCRs are multi-RAS (they do not target the wild-type)

drugs[Type %in% c("Vaccine","TCR") & is.na(Target)] -> x
drugs[Type %in% c("Vaccine","TCR") & named_mutations=="G12V",Target:="G12V"]
drugs[Type %in% c("Vaccine","TCR") & grepl("-12",name),Target:=c("multi-KRAS","G12V","G12C")]
drugs[code %in% c("C204252","C162186","C162269","C165638","C179231","C200465","C2067","C202009","C204884"),Target:="multi-KRAS"]
drugs[code == "C29136",Target:="G13N"]

# Pan-Ras

drugs[is.na(Target),list(name,definition,code)] -> x

drugs[code %in% c("C178264","C185876","C188048","C71146","C209714","C213206"),Target:="pan-RAS"]





# short names for figures (manually curated) ----

drugs[,short.Name := NULL]
drugs <- merge(drugs,fread("additional_files/short.names.csv")[,list(code,short.Name)],all.x=T)



# additional drugs, manually curated ----

# Note: it is unclear whether LY4066434 is G12D or pan-KRAS

additional.drugs <- fread("additional_files/drug.revision.csv")

intersect(additional.drugs$synonyms %>% strsplit(.,split="\\|") %>% unlist(),
          drugs$synonyms %>% strsplit(.,split="\\|") %>% unlist()) %>% paste(.,collapse="|") -> already.included

additional.drugs[,`:=`(From="Manually Added",children=0)]

drugs <- rbind(drugs,additional.drugs,fill=T)



# end creating the drugs data.table -------

drugs[,table(From,Group)]
drugs[is.na(code),code:="N.A."]

x <- drugs[Target!="BRAF"][order(Target,Type),list("Thesaurus\ncode"=code,"Name"=short.Name,Target,Type)]
dim(x)
png("output/drugs_table.png",1200,1400)
g1 <- tableGrob(x[1:45],rows=NULL,theme=ttheme_default(base_size=20))
g2 <- tableGrob(x[46:90],rows=NULL,theme=ttheme_default(base_size=20))
grid.arrange(g1,g2,layout_matrix=matrix(1:2,ncol=2))
dev.off()

# Clinical trials for colorectal cancer or solid tumors ----
# I haven't found a more elegant way of downloading the clincal trials data
# Downloading clinical trials data requires the input from the user...

browseURL("https://clinicaltrials.gov/search?cond=((colon%20OR%20rectal%20OR%20colorectal)%20AND%20cancer)%20OR%20(solid%20tumors)")
# download ctg-studies 
file.copy(from = "~/Downloads/ctg-studies.csv",to = "additional_files/ctg-studies.csv",overwrite = T)

ctrials <- fread("additional_files/ctg-studies.csv")
setkey(ctrials,"NCT Number")

# correct some misspellings
ctrials[grep("LY4066434.",Interventions),Interventions:=gsub("LY4066434\\.","LY4066434",Interventions)] 
ctrials[grep("5Fluoro",Interventions,ignore.case = T),Interventions:=gsub("5Fluoro","5-Fluoro",Interventions,ignore.case = T)] 

message("Loaded ",dim(ctrials)[1]," clinical trials") 


# search in clinical trials

# Search for clinical trials on CRC (or gastrointestinal) cancers

crc <- ctrials[grepl("colon|rectal|gastrointestinal|crc|solid",Conditions,ignore.case = T) & `Study Type`=="INTERVENTIONAL"]
crc <- crc[grep("TREATMENT",`Study Design`)]
crc <- rbind(crc,ctrials[grep("NCT06385925|NCT04330664|NCT03785249",`NCT Number`)])

drug.columns <- c("Drug.Code","Drug.Name","Drug.Type","Drug.Target","Drug.Group")
crc[,(drug.columns) := ""]

for(i in 1:nrow(drugs)) {
  hits <- grep(drugs$synonyms[i],crc$Interventions,ignore.case = T)
  crc[hits,Drug.Code:=paste(Drug.Code,drugs$code[i],sep="|")]
  crc[hits,Drug.Name:=paste(Drug.Name,drugs$short.Name[i],sep="|")]
  crc[hits,Drug.Type:=paste(Drug.Type,drugs$Type[i],sep="|")]
  crc[hits,Drug.Target:=paste(Drug.Target,drugs$Target[i],sep="|")]
  crc[hits,Drug.Group:=paste(Drug.Group,drugs$Group[i],sep="|")]
  }

# remove the initial |

crc[,(drug.columns) := lapply(.SD,function(x) gsub("^\\|","",x)),.SDcols=drug.columns]

dim(crc)

kras.trials <- crc[grep("KRAS|SHP2|SOS1",Drug.Group)]
message("KRAS Trials n=",nrow(kras.trials))

kras.trials[,(drug.columns) := ""]
drugs2 <- drugs[Target!="BRAF"]

for(i in 1:nrow(drugs2)) {
  hits <- grep(drugs2$synonyms[i],kras.trials$Interventions,ignore.case = T)
  kras.trials[hits,Drug.Code:=paste(Drug.Code,drugs2$code[i],sep="|")]
  kras.trials[hits,Drug.Name:=paste(Drug.Name,drugs2$short.Name[i],sep="|")]
  kras.trials[hits,Drug.Type:=paste(Drug.Type,drugs2$Type[i],sep="|")]
  kras.trials[hits,Drug.Target:=paste(Drug.Target,drugs2$Target[i],sep="|")]
  kras.trials[hits,Drug.Group:=paste(Drug.Group,drugs2$Group[i],sep="|")]
}

kras.trials[,(drug.columns) := lapply(.SD,function(x) gsub("^\\|","",x)),.SDcols=drug.columns]

kras.trials[Drug.Target!=""]

# trials in review 2024
inReview <- c("NCT05485974","NCT05462717","NCT06244771","NCT04121286","NCT06385925","NCT05737706","NCT06364696","NCT06403735","NCT06040541","NCT06412198","NCT05194995","NCT04330664","NCT05288205","NCT06039384","NCT05198934","NCT04956640","NCT05578092","NCT04699188","NCT04449874","NCT05358249","NCT06026410","NCT06252649","NCT06586515","NCT04793958","NCT03785249","NCT05722327","NCT06599502","NCT06078800","NCT06607185","NCT06447662","NCT06585488","NCT06445062","NCT05379985","NCT05786924","NCT06194877","NCT05200442","NCT06270082","NCT05163028","NCT04117087","NCT04853017","NCT06411691","NCT05726864","NCT06105021","NCT06253520","NCT06487377","NCT06218914")
inReview <- ctrials[inReview]

setdiff(inReview$`NCT Number`,crc[Drug.Type!="",`NCT Number`]) %>% ctrials[.] -> missed
dim(missed) # two trials were not retrieved by using the previous approach

# check missed trials

x <- ctrials[grepl("k[-]*ras",Interventions,ignore.case = T) & ! `NCT Number` %in% kras.trials$`NCT Number` &
               grepl("TREATMENT",`Study Design`)]
dim(x) # 7 trials putatively missed

# Targeting KRAS: 1,2,4,5,6,7

x <- rbind(kras.trials,x[c(1,2,4,5,6,7)],fill=T)
x2 <- fread("data/curated.kras.trials.csv") # manually curated

x$`NCT Number` == x2$`NCT Number`

x$Drug.Code <- x2$Drug.Code
x$Drug.Name <- x2$Drug.Name
x$Drug.Type <- x2$Drug.Type
x$Drug.Target <- x2$Drug.Target
x$Drug.Group <- x2$Drug.Group

kras.trials <- x
rm(x,x2)

table(kras.trials$Drug.Type)
table(kras.trials$Drug.Target)
# review trials with more than one identified drug


# save for manual revision 

# fwrite(kras.trials,file="data/curated.kras.trials.csv")

# fread the revised version

try(kras.trials <- fread("data/curated.kras.trials.csv"),silent=F)

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

table(kras.trials$Phases)


x <- kras.trials
x[completion > `Last Update Posted`,last:=completion]
x[completion <= `Last Update Posted`,last:=`Last Update Posted`]
x[,target := Drug.Target %>% split1 %>% unlist %>% sort %>% unique %>% join1,by=`NCT Number`]
x <- x[order(target,start)]
x[,`NCT Number`:=factor(`NCT Number`,levels=as.character(`NCT Number`))]



trial.plot <- list(
  aes(start,`NCT Number`),
  geom_vline(xintercept = as.Date("2025-03-15"),lwd=3,color="grey"),
  geom_segment(aes(yend=`NCT Number`,xend=completion),lwd=2),
  geom_point(aes(color=`Study Status`)),
  xlab(""),
  ylab(""),
  scale_x_date(limits=as.Date(c("2010-01-01","2045-12-31"))),
  scale_y_discrete(limits=rev),
  geom_point(aes(x=as.Date(`Last Update Posted`)),pch=21,lwd=0.2,fill="white"),
  geom_text(aes(label=target,x=start),hjust=1,nudge_x=-365,size=3),
  theme(strip.text.x.top = element_text(size=12,face = "bold"),
        panel.background = element_rect(fill="white"))
  
)

active <- ! x$`Study Status` %in% c("COMPLETED","TERMINATED")

g0 <- ggplot(x[!active]) + trial.plot +
  facet_wrap(~"Completed or terminated") +   
  scale_x_date(limits=as.Date(c("2006-01-01","2030-12-31"))) +
  theme(legend.position = "none")


g1 <- ggplot(x[active & Phases=="PHASE1"]) + trial.plot + facet_wrap(~ "Phase I") + theme(legend.position = "none")


g2 <- ggplot(x[active & grepl("PHASE2",Phases)]) + trial.plot + facet_wrap(~ "Phase I/II or II") + theme(legend.position = "none")

g3 <- ggplot(x[active & Phases=="PHASE3"]) + trial.plot + facet_wrap(~ "Phase III") + 
  geom_text(aes(label=sprintf("%s (%s)",Acronym,Drug.Name),x=last),hjust=0,nudge_x = 365) 

grid.arrange(g0,g1,g2,g3,layout_matrix=matrix(c(2,3,1,2,4,4),byrow=T,ncol=3),heights=c(9,2),widths=c(3,3,2))



