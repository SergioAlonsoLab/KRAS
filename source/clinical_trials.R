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


# Create a data.table with KRAS related drugs ------

# add some missing synonyms for RMC drugs

thesaurus[grep("RMC-9805",synonyms),synonyms:=paste0("Zoldonrasib",synonyms)]
thesaurus[grep("RMC-6236",synonyms),synonyms:=paste0("Daraxonrasib",synonyms)]
thesaurus[grep("RMC-6291",synonyms),synonyms:=paste0("Elironrasib",synonyms)]

# Add Ras inhibitors (code C19902)

drugs <- thesaurus[findChildren("C1902")$code][children==0][,Group:="KRAS"] %>% data.table

# Include Diazopinomicin

x <- thesaurus["C62508"][,Group:="KRAS"] %>% data.table
drugs <- rbind(drugs,x)

# SOS1 and SHP2

x <- thesaurus[isDrug][grepl("sos[-]*[12]",paste(synonyms,definition),ignore.case = T)]
x[,Group:="SOS1"]

drugs <- rbind(drugs,x)

x <- thesaurus[findChildren("C185612")$code]
x[,Group:="SHP2"]

drugs <- rbind(drugs,x)

# RAS TCRs

x <- thesaurus[isDrug][grepl("TCR",paste(synonyms,definition)) & grepl("[ NKH][-]*RAS",paste(synonyms,definition),ignore.case = T)]
x[,Group:="KRAS"]
x[,Type:="TCR"]
drugs <- rbind(drugs,x,fill=T)

# RAS vaccines

x <- thesaurus[grepl("vaccine",paste(synonyms,definition),ignore.case = T) &
                 grepl("[ KNH][-]*RAS",paste(synonyms,definition),ignore.case = T)][name!="Measles" & children==0]

x[,Group:="KRAS"]
x[,Type:="Vaccine"]

# this retrieves also peptides and proteins that "may" be used as vaccines, but are not actual vaccines

x[grepl("Peptide",name) & !grepl("Vaccine",name),Type:="Mutant Peptide"]
x[grepl("Protein",name) & !grepl("Vaccine",name),Type:="Mutant Protein"]

x
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

# Add the target for every compound

drugs[,Target:=""]

#G12C

drugs[match(findChildren("C201853")$code,code),Target:="G12C"]
drugs[grep("fulzerasib|divarasib|talorasib",synonyms,ignore.case = T),Target:="G12C"]

drugs[grep("G12C",paste(synonyms,definition),ignore.case = T)]

#G12D
drugs[match(findChildren("C207236")$code,code),Target:="G12D"]
drugs[named_mutations=="G12D" & Target=="",Target:="G12D"]

#SOS1 or SHP2
drugs[grep("SOS1",Group),Target:="SOS1"]
drugs[grep("SHP2",Group),Target:="SHP2"]

#OTHER
drugs[grep("Diazepinomicin",name),`:=`(Type="Inhibitor",Target="pan-RAS")]

# Check drugs without assigned target
drugs[Target==""]


#Classify drug type
drugs[grep("inhibitor",name,ignore.case = T),Type:="Inhibitor"]
drugs[grep("(nib|fib|sib|tib|ratide)( |$)",name),Type:="Inhibitor"]

drugs[grep("degrader",paste(name,synonyms,sep = "|"),ignore.case = T),Type:="Degrader"]

drugs[grep("Glue",paste(name,synonyms,sep = "|"),ignore.case = T),Type:="Molecular Glue"]
drugs[grep("Antisense",paste(name,synonyms,sep = "|"),ignore.case = T),Type:="Antisense"]

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

drugs[Type=="Vaccine",paste(code,name,definition,synonyms,named_mutations,sep = "\n")] %>% unlist %>% cat(sep="\n--------------\n")


# Vaccines and TCRs are multi-RAS (they do not target the wild-type)

drugs[Type %in% c("Vaccine","TCR") & is.na(Target)] -> x
drugs[Type %in% c("Vaccine","TCR") & named_mutations=="G12V",Target:="G12V"]
drugs[Type %in% c("Vaccine","TCR") & grepl("-12",name),Target:=c("multi-KRAS","G12V","G12C")]
drugs[code %in% c("C204252","C162186","C162269","C165638","C179231","C200465","C2067","C202009","C204884"),Target:="multi-KRAS"]
drugs[code == "C29136",Target:="G13N"]

# Pan-Ras
drugs[Target=="",list(name,definition,code)] -> x
drugs[code %in% c("C178264","C185876","C188048","C71146","C209714","C213206"),Target:="pan-RAS"]
drugs[Target==""]

# short names for figures (manually curated) ----

drugs[,short.Name := NULL]
drugs <- merge(drugs,fread("additional_files/short.names.csv")[,list(code,short.Name)],all.x=T)

# additional drugs, manually curated ----
# Note: it is unclear whether LY4066434 is G12D or pan-KRAS

additional.drugs <- googlesheets4::read_sheet("1NOOz5opFwZg4AAYhkCLXoGnj0EUnmpx1pffT5DhqdAs",1) %>% as.data.table
additional.drugs[Code=="",Code:="NA"]
names(drugs)
names(additional.drugs) <- c("code","From","name","short.Name","synonyms","Target","Group","Type")

additional.drugs[! short.Name %in% drugs$short.Name]
drugs[! code %in% additional.drugs$code,list(code,Group)]

drugs <- rbind(drugs,additional.drugs[From=="Manually Added"],fill=T)

# add data from PubChem -----

lapply(dir("additional_files/",pattern="PubChem_compound*"),
       function(x) {
         d0 <- fread(paste0("additional_files/",x))
         d0[,File:=gsub("PubChem_compound_text_","",x) %>% gsub(".csv","",.)]
         return(d0)
       }) %>% rbindlist() -> PubChem

# Remove duplicated entries in PubChem

PubChem <- PubChem[!duplicated(`Compound CID`)]
PubChem[,Synonyms:=paste(Name,Synonyms,sep="|"),by=`Compound CID`]

# Remove "Carbon"

PubChem <- PubChem[Name!="Carbon"]

# Cross pubchem with drugs table

x <- PubChem[,list(Synonym=strsplit(Synonyms,split="\\|") %>% unlist),by=list(`Compound CID`,Name,`Molecular Weight`)] %>% unique

# Do different compound share a synonym?

x[,Compounds:=paste(`Compound CID`,collapse="|"),by=Synonym]
x[grep("\\|",Compounds),list(`Compound CID`,Synonym,Compounds)][order(Synonym)]

# Remove PubChem compounds already included in the drugs table

x[,Hits:=drugs[grep(tolower(Synonym),tolower(synonyms),fixed = T),short.Name] %>% paste(collapse="|"),by=Synonym]
x[Hits!="",list(`Compound CID`,Name,Hits)] %>% unique -> included
included
PubChem <- PubChem[! `Compound CID` %in% included$`Compound CID`]

# Clinical trials for colorectal cancer or solid tumors --------

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

# correct interventions in some clinical trials

ctrials["NCT06128551",Interventions:="DRUG: RMC-6291|DRUG: RMC-6236"]
ctrials["NCT06619587",Interventions:="DRUG: GDC-7035"]

# Search for clinical trials on CRC (or gastrointestinal) cancers

crc <- ctrials[grepl("colon|rectal|gastrointestinal|crc|solid",Conditions,ignore.case = T) & `Study Type`=="INTERVENTIONAL"]
crc <- crc[grep("TREATMENT",`Study Design`)]

# trials in review 2024
inReview <- c("NCT05485974","NCT05462717","NCT06244771","NCT04121286","NCT06385925","NCT05737706","NCT06364696","NCT06403735","NCT06040541","NCT06412198","NCT05194995","NCT04330664","NCT05288205","NCT06039384","NCT05198934","NCT04956640","NCT05578092","NCT04699188","NCT04449874","NCT05358249","NCT06026410","NCT06252649","NCT06586515","NCT04793958","NCT03785249","NCT05722327","NCT06599502","NCT06078800","NCT06607185","NCT06447662","NCT06585488","NCT06445062","NCT05379985","NCT05786924","NCT06194877","NCT05200442","NCT06270082","NCT05163028","NCT04117087","NCT04853017","NCT06411691","NCT05726864","NCT06105021","NCT06253520","NCT06487377","NCT06218914")
missed <- setdiff(inReview,crc$`NCT Number`)
crc <- rbind(crc,ctrials[missed])

# add CTrials to compounds tables 

x <- drugs[,list(Synonym=strsplit(synonyms,split="\\|") %>% unlist %>% tolower),by=name] %>% unique
x <- x[nchar(Synonym) >= 4]
x[,CTrials:=crc[grep(Synonym,tolower(Interventions),fixed = T),paste(`NCT Number`,collapse="|")],by=Synonym]
x[order(nchar(Synonym))]
x <- x[,list(CTrials=paste(CTrials[CTrials!=""],collapse="|")),by=name]
x


drugs <- merge(drugs,x,all.x=T)
dim(drugs)


x <- PubChem[,list(Synonym=strsplit(Synonyms,split="\\|") %>% unlist %>% tolower),by=`Compound CID`] %>% unique
x <- x[nchar(Synonym) >= 4]
x[,CTrials:=crc[grep(Synonym,tolower(Interventions),fixed = T),paste(`NCT Number`,collapse="|")],by=Synonym]
x[order(nchar(Synonym))]

dim(x[CTrials!=""]) # NONE!!


# ADD COMPOUNDS TO CTrials ------

x <- drugs[CTrials!="",list(short.Name,Type,Group,Target,`NCT Number`=strsplit(CTrials,split="\\|") %>% unlist),by=name] %>% unique
x <- x[,list(Drug=paste(short.Name,collapse="|"),Type=paste(Type,collapse="|"),Group=paste(Group,collapse="|"),Target=paste(Target,collapse="|")),by=`NCT Number`] 
x[grep("\\|",Drug)]
kras.trials <- merge(crc,x)

message("KRAS Trials n=",nrow(kras.trials))

missed <- setdiff(inReview,kras.trials$`NCT Number`)
ctrials[missed,list(`NCT Number`,`Study Title`,Interventions)]

message("The missed clincal trials are testing BRAF / MEK inhibitors: BDTX-4933, BGB-3245, VS-6766, and IK-595")

# check other possibly missed trials (KRAS,SOS1 or SHP2 in interventions)

x <- ctrials[grepl("(k[-]*ras|sos[-]*1|shp[-]*2)",paste(Interventions),ignore.case = T) & ! `NCT Number` %in% kras.trials$`NCT Number` &
               grepl("TREATMENT",`Study Design`)]


x[,list(`NCT Number`,Interventions)] # NCT04745130 is not targeting KRAS: Sintilimab，regofinib，cetuximab

x <- x[`NCT Number`!="NCT04745130"]

message("Missed with K-ras or KRAS in interventions: ",nrow(x))

# check other possible missed trials

crc[grepl("k[-]*ras",paste(`Brief Summary`),ignore.case = T) & 
      !(`NCT Number` %in% kras.trials$`NCT Number`) &
      !grepl("wild[ -]*type",`Brief Summary`,ignore.case = T)] -> toReview

reviewed <- googlesheets4::read_sheet("1NOOz5opFwZg4AAYhkCLXoGnj0EUnmpx1pffT5DhqdAs",4) %>% as.data.table
reviewed <- reviewed[Include=="YES" & !(`NCT Number` %in% kras.trials$`NCT Number`),list(`NCT Number`,Interventions)]

dim(reviewed)

kras.trials <- rbind(kras.trials,crc[match(reviewed$`NCT Number`,`NCT Number`)],fill=T)

# Manually fix
reviewed$`NCT Number`
kras.trials[`NCT Number`=="NCT05202561",`:=`(Type="Vaccine",Group="KRAS",Drug="RNA Vaccine",Target="multi-KRAS")]
kras.trials[`NCT Number`=="NCT06043713",`:=`(Type="TCR",Group="KRAS",Drug="FH-A11KRASG12V-TCR",Target="G12V")]


# fix Acronyms for Phase III studies

kras.trials[match(c("NCT04793958","NCT05198934","NCT6252649"),`NCT Number`),Acronym:=c("KRYSTAL-10","CodeBreaK 300","CodeBreaK 301")]

# add bibligraphic information 
library(rentrez)

entrez_search("pmc","NCT03600883",retmax=500)$ids

sapply(kras.trials$`NCT Number`,function(x) {
  entrez_search("pubmed",x,retmax=500)$ids %>% unlist -> ids
  message(x," hits in PubMed: ",length(ids))
  return(paste(sort(ids),collapse="|"))
}) -> pmids

kras.trials[,pmids:=pmids]

sapply(kras.trials$`NCT Number`,function(x) {
  entrez_search("pmc",sprintf("%s[Abstract] OR %s[Title] OR %s[Body - Key Terms]",x,x,x),retmax=500)$ids %>% unlist -> ids
  if(length(ids) > 0) ids <- paste0("PMC",ids)
  message(x," hits in PMC: ",length(ids))
  return(paste(sort(ids),collapse="|"))
}) -> pmcids

kras.trials[,pmcids:=pmcids]

kras.trials[,CRC_in_Conditions:=grepl("colon|rectal|crc",Conditions,ignore.case = T)]
kras.trials[,CRC_in_Summary:=grepl("colon|rectal|crc",`Brief Summary`,ignore.case = T)]

fwrite(drugs,"data/kras.drugs.csv")
fwrite(kras.trials,"data/kras.clinical.trials.csv")
fwrite(PubChem,"data/pubchem.csv")
