#question1
#E.COLI GENE SEQUENCE
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
#DECOMPRESS
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)
#Blast database
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")
#Sequences present
str("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
seqinr::getLength("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
#question2
#download and reading sequence length festa file
E.COLI <-fasta.seqlengths("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa")
str(E.COLI)
#read the festa file
E.COLI<-read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
E.COLI <- readDNAStringSet("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
#selecting gene of interest
geneofinterest <-extractROWS(E.COLI,10)
str(geneofinterest)
as.list(geneofinterest)
#determining the length
str(E.COLI,10)
DNAStringSet(E.COLI,10)
as.list(geneofinterest)
#proportion of GC
extractseqs(E.COLIseq)
#question3
#performing blast
myblastn_tab
res <- rBLAST::makeblastdb ("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
#sequence matching
db <- read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(db[1:10])
head(names(db))
myseqs <- db[which(names(db) %in% Hits)]
myseqs <- c(myseqs,E.COLI) 

      

