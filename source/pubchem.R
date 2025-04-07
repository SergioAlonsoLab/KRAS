lapply(dir("additional_files/",pattern="PubChem_*"),
       function(x) {
         d0 <- fread(paste0("additional_files/",x))
         d0[,File:=x]
         return(d0)
       }) %>% rbindlist() -> PubChem

# Remove duplicated entries

PubChem <- PubChem[!duplicated(`Compound CID`)]
PubChem[,Synonyms:=paste(Name,Synonyms,sep="|"),by=`Compound CID`]

# cross with crc clinical trials

escapeALL <- function(x) {
  gsub("(\\(|\\)|\\[|\\]|\\{|\\})", ".", x,fixed = F)
}

shortSynonyms <- function(x) {
  x <- strsplit(x,split="\\|") %>% unlist
  x <- x[nchar(x) < 25 & nchar(x) > 4]
  x2 <- gsub("[-]+","",x)
  
  paste(unique(toupper(c(x,x2))),collapse="|")
}

PubChem[,Synonyms2 := Synonyms %>% escapeALL %>% shortSynonyms,by=`Compound CID`]
PubChem[,CTrials := crc[grep(Synonyms2,Interventions,ignore.case = T),paste(`NCT Number`,collapse="|")],by=`Compound CID`]
PubChem[,NTrials := strsplit(CTrials,split="\\|") %>% unlist %>% length,by=`Compound CID`]

PubChem[NTrials>=20] # "carbon
PubChem <- PubChem[NTrials < 20]

fwrite(PubChem,file="additional_files/PubChem.csv")


