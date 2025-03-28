library(data.table)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(gridExtra)
library(igraph)
library(stringr)
library(httr)
library(editData)
library(googlesheets4)

# useful functions to expand and contract terms

join1 <- function(x) paste(x,collapse="|")
split1 <- function(x) strsplit(x,split="\\|") 


# read the Thesaurus database

if(!file.exists("data/Thesaurus.txt")) {
  
  curl::curl_download("https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus.FLAT.zip",destfile = "data/Thesaurus.zip")
  
}

thesaurus <- fread("data/Thesaurus.zip",quote="")
names(thesaurus) <- c("code","concept IRI","parents","synonyms","definition","display name","concept status","semantic type","concept in subset")
thesaurus$synonyms %>% split1 %>% sapply(FUN=function(x) x[1]) -> thesaurus$name # the preferred name
setkey(thesaurus,"code")
tsearch <- function(x) {
  a <- thesaurus[tolower(x)==tolower(name)]
  if(nrow(a)==0) a <- thesaurus[grep(x,synonyms,ignore.case = T)]
  return(a)
}
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
drugs[grep("SOS1",Group),Target:="SOS1"]
drugs[grep("SHP2",Group),Target:="SHP2"]


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

drugs[grep("ZG2001",synonyms),synonyms:=paste0(synonyms,"|ZG2001|ZG-2001")]


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

additional.drugs <- fread("additional_files/manually_added.csv")
additional.drugs[Code=="",Code:="NA"]

names(drugs)
names(additional.drugs) <- c("code","From","name","short.Name","synonyms","Target","Group","Type")

drugs[! code %in% additional.drugs$code,list(code,Group)]

drugs <- rbind(drugs,additional.drugs[From=="Manually Added"],fill=T)

# end creating the drugs data.table -------

drugs[,table(From,Group)]
x <- drugs[Target!="BRAF"][order(Target,Type),list("NCI\ncode"=code,"Name"=short.Name,Target,Type)]
dim(x)
png("output/drugs_table.png",1200,1600)
g1 <- tableGrob(x[1:52],rows=NULL,theme=ttheme_default(base_size=20))
g2 <- tableGrob(x[53:103],rows=NULL,theme=ttheme_default(base_size=20))
grid.arrange(g1,g2,layout_matrix=matrix(1:2,ncol=2))
dev.off()

# Clinical trials for colorectal cancer or solid tumors ----

if(!file.exists("data/ctg-studies.csv")) {
  
  url_base <- "https://clinicaltrials.gov/api/v2/studies?"
  url_query <- "query.cond=((colon+OR+rectal+OR+colorectal+OR+CRC)+AND+cancer)+OR+(solid+tumors)"
  query_url <- paste0(url_base,url_query,"&format=csv&pageSize=500")
  
  message("Downloading registries 1 to 500")
  n <- 0
  query <- GET(query_url)
  ctrials <- fread(content(query,as="text",encoding = "UTF-8"))
  colnames_trials <- colnames(ctrials)
  nextToken <- query$headers$`x-next-page-token`
  
  while(!is.null(nextToken)) {
    n <- n + 1
    message("Downloading registries ",500*n+1," to ",500*(n+1))
    query_url <- paste0(url_base,url_query,"&format=csv&pageSize=500&pageToken=",nextToken)
    query <- GET(query_url)
    x <- fread(content(query,as="text",encoding = "UTF-8"))
    colnames(x) <- colnames_trials
    ctrials <- rbind(ctrials,x)
    nextToken <- query$headers$`x-next-page-token`
  }
  
  fwrite(ctrials,file="data/ctg-studies.csv")
  
}

ctrials <- fread("data/ctg-studies.csv")
setkey(ctrials,"NCT Number")

# correct some misspellings

ctrials[grep("5Fluoro",Interventions,ignore.case = T),Interventions:=gsub("5Fluoro","5-Fluoro",Interventions,ignore.case = T)] 

message("Loaded ",nrow(ctrials)," clinical trials") 

# Search for clinical trials on CRC (or gastrointestinal) cancers

crc <- ctrials[grepl("colon|rectal|gastrointestinal|crc|solid",Conditions,ignore.case = T) & `Study Type`=="INTERVENTIONAL"]
crc <- crc[grep("TREATMENT",`Study Design`)]

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

# select trials targeting KRAS, SHP2 or SOS1

kras.trials <- crc[grep("KRAS|SHP2|SOS1",Drug.Group)]
message("KRAS Trials n=",nrow(kras.trials))

# trials in review 2024
inReview <- c("NCT05485974","NCT05462717","NCT06244771","NCT04121286","NCT06385925","NCT05737706","NCT06364696","NCT06403735","NCT06040541","NCT06412198","NCT05194995","NCT04330664","NCT05288205","NCT06039384","NCT05198934","NCT04956640","NCT05578092","NCT04699188","NCT04449874","NCT05358249","NCT06026410","NCT06252649","NCT06586515","NCT04793958","NCT03785249","NCT05722327","NCT06599502","NCT06078800","NCT06607185","NCT06447662","NCT06585488","NCT06445062","NCT05379985","NCT05786924","NCT06194877","NCT05200442","NCT06270082","NCT05163028","NCT04117087","NCT04853017","NCT06411691","NCT05726864","NCT06105021","NCT06253520","NCT06487377","NCT06218914")
inReview <- ctrials[inReview]

setdiff(inReview$`NCT Number`,crc[Drug.Type!="",`NCT Number`]) %>% ctrials[.] -> missed

message(nrow(missed)," clinical trials not retrieved by the automatic method, manually added")

kras.trials <- rbind(kras.trials,missed,fill=TRUE)

# reassign the drug information using only drugs targeting KRAS, SOS1 or SHP2
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

kras.trials[Drug.Target=="",list(`NCT Number`,Interventions)]

# check missed trials

x <- ctrials[grepl("k[-]*ras",paste(Interventions),ignore.case = T) & ! `NCT Number` %in% kras.trials$`NCT Number` &
               grepl("TREATMENT",`Study Design`)]


x[,list(`NCT Number`,Interventions)] # NCT04745130 is not targeting KRAS: Sintilimab，regofinib，cetuximab

x <- x[`NCT Number`!="NCT04745130"]

dim(x)

kras.trials <- rbind(kras.trials,x,fill=T)

# all the selected KRAS-related drugs

x <- kras.trials
x[,Synonyms:=drugs2[match(split1(Drug.Code) %>% unlist,code),synonyms %>% join1],by=`NCT Number`]
x$Synonyms


# Checking other possible missed trials

crc[grepl("k[-]*ras",paste(`Brief Summary`),ignore.case = T) & 
      !(`NCT Number` %in% kras.trials$`NCT Number`) &
      !grepl("wild[ -]*type",`Brief Summary`,ignore.case = T)] -> toReview

fwrite(drugs2[,list(Code=code,From,Name=name,"Short Name"=short.Name,synonyms,Target,Group,Type)],file="additional_files/drugs.kras.csv")

lapply(drugs2$synonyms %>% strsplit(split="\\|") %>% unlist %>% toupper %>% unique,function(i) {
  print(i)
  drug <- drugs2[grep(i,synonyms,ignore.case = T),]
  selected <- kras.trials[grep(i,Interventions,ignore.case = T),list(`NCT Number`)]
  data.table("Identified term"=i,drug[,list("Drug Name"=name,Target,Type,Group)],selected)
}) %>% rbindlist %>% na.omit -> foo

merge(kras.trials[,list(`NCT Number`,`Study Title`,`Study URL`,Phases,Interventions)],foo,all.x=T) -> foo
foo[,`Identified term` := sort(`Identified term`)[1],by=list(`NCT Number`,`Drug Name`)]
foo <- unique(foo)

fwrite(foo,file="additional_files/included.csv")

toReview[,strsplit(Interventions,split="(: |\\|)") %>% unlist %>% gsub("(DRUG|PROCEDURE|  )","",.) %>% paste(collapse=", "),by=`NCT Number`]

fwrite(toReview[,list(Include="TBD",`NCT Number`,`Study Title`,Interventions,`Study URL`)],file = "additional_files/toReview.csv")



# manually review the clinical trial table, to verify the information automatically 
# added and the last 6 clinical trials added

fwrite(kras.trials[,list(`NCT Number`,`Study Title`,Acronym,Interventions,Drug.Code,Drug.Name,Drug.Target,Drug.Group)],
       file="data/toReview.kras.trials.csv")

# 


message("Please review ")

# fread the revised version

try(kras.trials <- fread("data/curated.kras.trials.csv"),silent=F)
kras.trials <- merge(ctrials,kras.trials[,list(`NCT Number`,Drug.Code,Drug.Name,Drug.Target,Drug.Group)])

# fix Acronyms for Phase III studies

kras.trials[match(c("NCT04793958","NCT05198934","NCT6252649"),`NCT Number`),Acronym:=c("KRYSTAL-10","CodeBreaK 300","CodeBreaK 301")]

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

kras.trials[nchar(`Primary Completion Date`) > 8,completion:=as.Date(`Primary Completion Date`)]
kras.trials[nchar(`Primary Completion Date`) == 7,completion:=as.Date(paste0(`Primary Completion Date`,"-01"))]
kras.trials[is.na(completion)]

my.palette <- list()
my.palette$study.status <- c(TERMINATED="#9b2226",
                             COMPLETED="#E73D23",
                             UNKNOWN="#555555",
                             ACTIVE_NOT_RECRUITING="#438BC7",
                             RECRUITING="#3AA691",
                             ENROLLING_BY_INVITATION="#59C5AF",
                             NOT_YET_RECRUITING="#95DACC")


x <- kras.trials
x[completion > `Last Update Posted`,last:=completion]
x[completion <= `Last Update Posted`,last:=`Last Update Posted`]
x[,target := Drug.Target %>% split1 %>% unlist %>% sort %>% unique %>% join1,by=`NCT Number`]
x <- x[order(Drug.Target,start)]
x[,`NCT Number`:=factor(`NCT Number`,levels=as.character(`NCT Number`))]

x$`Study Status` <- factor(x$`Study Status`,names(my.palette$study.status))

trial.plot <- list(
  aes(start,`NCT Number`),
  geom_vline(xintercept = as.Date("2025-03-15"),lwd=3,color="lightgrey"),
  geom_segment(aes(yend=`NCT Number`,xend=`Last Update Posted`),lwd=0.5,lty=5),
  geom_segment(aes(yend=`NCT Number`,xend=completion,color=`Study Status`),lwd=2,show.legend=T),
  xlab(""),
  ylab(""),
  scale_x_date(limits=as.Date(c("2012-01-01","2045-12-31"))),
  scale_y_discrete(limits=rev),
  geom_point(aes(x=as.Date(`Last Update Posted`)),pch=21,fill="white"),
  geom_text(aes(label=target,x=start),hjust=1,nudge_x=-60,size=3),
  scale_color_manual("Study Status",values=my.palette$study.status,drop=F),
  theme(strip.text.x.top = element_text(size=16,face = "bold"),
        panel.background = element_rect(fill="white",color="black"),
        panel.grid = element_line(color="#EFEFEF"),
        strip.background = element_rect(color="black"),
        legend.key.height = unit(10, "pt"))
  
)

active <- ! x$`Study Status` %in% c("COMPLETED","TERMINATED")

g0 <- ggplot(x[!active]) + trial.plot +
  facet_wrap(~"Completed or terminated") +   
  scale_x_date(limits=as.Date(c("2000-01-01","2035-12-31"))) +
  theme(legend.position = "none") +
  geom_text(aes(label=sprintf("%s",Drug.Name),x=last),hjust=0,nudge_x = 240,size=3.2) 

g1 <- ggplot(x[active & Phases=="PHASE1"]) + trial.plot + facet_wrap(~ "Phase I") + 
  theme(legend.position = "none") +
  geom_text(aes(label=sprintf("%s",Drug.Name),x=last),hjust=0,nudge_x = 60,size=3.2) 


g2 <- ggplot(x[active & grepl("PHASE2",Phases)]) + trial.plot + facet_wrap(~ "Phase I/II or II") + 
  theme(legend.position = "none") +
  geom_text(aes(label=sprintf("%s",Drug.Name),x=last),hjust=0,nudge_x = 60,size=3.2) 


g3 <- ggplot(x[active & Phases=="PHASE3"]) + trial.plot + facet_wrap(~ "Phase III") + 
  geom_text(aes(label=sprintf("%s (%s)",Drug.Name,Acronym),x=last),hjust=0,nudge_x = 60)  

g4 <- grid.arrange(g0,g1,g2,g3,layout_matrix=matrix(c(2,3,1,2,4,4),byrow=T,ncol=3),heights=c(9,2),widths=c(3,3,2.5))

g3 <- ggplot(x[active & Phases=="PHASE3"]) + trial.plot + facet_wrap(~ "Phase III") + 
  geom_text(aes(label=sprintf("%s (%s)",Drug.Name,Acronym),x=last),hjust=0,nudge_x = 60) +
  scale_color_manual("Study Status",values=my.palette$study.status[3:7],drop=F) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(direction = "vertical"))

g5 <- grid.arrange(g1,g2,g3,layout_matrix=matrix(c(1,2,1,3),byrow=T,ncol=2),heights=c(10,3))

ggsave(filename = "output/ClinicalTrials.tiff",plot = g4,width = 18,height = 10,dpi=200)
ggsave(filename = "output/ClinicalTrials.png",plot = g4,width = 18,height = 10,dpi=200)
ggsave(filename = "output/ClinicalTrials.pdf",plot = g4,width = 18,height = 10)

ggsave(filename = "output/ClinicalTrials2.tiff",plot = g5,width = 12,height = 14,dpi=200)
ggsave(filename = "output/ClinicalTrials2.png",plot = g5,width = 12,height = 14,dpi=200)
ggsave(filename = "output/ClinicalTrials2.pdf",plot = g5,width = 12,height = 14)


