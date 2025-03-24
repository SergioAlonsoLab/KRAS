
x <- bowel[order(patientId,SAMPLE_TYPE)][SAMPLE_COUNT > 1,list("Sample ID"=paste(sampleId,collapse="\n"),
                                                               "Sample Type"=paste(SAMPLE_TYPE,collapse="\n"),
                                                               "MSI status"=paste(GROUP,collapse="\n"),
                                                               KRAS=paste(KRAS,collapse="\n"),
                                                               NRAS=paste(NRAS,collapse="\n"),
                                                               HRAS=paste(HRAS,collapse="\n"),
                                                               BRAF=paste(BRAF,collapse="\n"),
                                                               APC=paste(APC,collapse="\n"),
                                                               TP53=paste(TP53,collapse="\n")),by=patientId] 
names(x)[1] <- "Patient ID"

tableS1 <- tableGrob(x,rows=NULL,theme=ttheme_default(base_size=10))
figuresHR(tableS1,11,16,fileBaseName = "TableS1")
