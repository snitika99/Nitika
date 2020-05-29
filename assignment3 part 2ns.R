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
#selecting gene of interest
geneofinterest <-extractROWS(E.COLI,10)
str(geneofinterest)
as.list(geneofinterest)
#determining the length
seqinr:: getLength(E.COLI)
#proportion of GC
extractseqs (geneofinterest)
#question3
#performing blast
myblastn_tab
res <- rBLAST::makeblastdb ("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
#sequence matching
res <- makeblastdb(myseqs=E.COLI, db="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(res)
res
head(res)
#determing first 3 hits
Hits <- as.character(res$seqid[1:3])
Hits
#best matching sequences
db <- read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(db)
str(db[1:6])
head(names(db))
#extraction of top hits
myseqs <- db[which(names(db) %in% Hits)]
myseqs <-c(myseqs,E.COLI)
seqinr:: write.fasta (myseq,names = names(myseqs), file.out="myseqs.fa")
str(myseqs)  
#names of sequence of first top last hit

