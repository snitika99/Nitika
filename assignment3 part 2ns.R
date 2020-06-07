library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")
library("GenomicRanges")
source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")

#question1

#E.COLI GENE SEQUENCE
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

#DECOMPRESS
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)

#Blast database
library("rBLAST")
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype = "nucl","-parse_seqids")

#Sequences present
str("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
seqinr::getLength("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

#question2

#download and reading sequence length festa file
E.COLI <-seqinr::read.fasta("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa")
str(E.COLI)

#read the festa file
E.COLI<-seqinr::read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")

#selecting gene of interest
geneofinterest <- E.COLI[[10]]

#determining the length
seqinr:: getLength(geneofinterest)

#proportion of GC
seqinr::GC(geneofinterest)

#question3
myblastn_tab
res <- myblastn_tab(myseq = E.COLI, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(res)
head(res)


#determing first 3 hits
Hits <- as.character(res$sseqid[1:3])
Hits

#best matching sequences
db <- read.fasta("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
str(db)
str(db[1:6])
head(names(db))

#extraction of top hits
myseqs <- db[which(names(db) %in% Hits)]
myseqs <-c(myseqs,E.COLI)
seqinr:: write.fasta (E.COLI,names = names(myseqs), file.out="myseqs.fa")
str(myseqs)

#extract the names of top hits 
tophit <- db[which(names(db) %in% Hits[1])]
tophit[1:3]

#bit scores and percent
seqinr::write.fasta(tophit,names=names(tophit),file.out = "tophit.fa")
makeblastdb("tophit.fa",dbtype="nucl", "-parse_seqids")
res <- myblastn(myseq = E.COLI, db = "tophit.fa")
cat(res,fill=TRUE)

#question 4 
#read the fasta file
E.coliismysequence <- readDNAStringSet("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
#coordinates of ORFs 
E.coliismysequence_orf <- ORFik::findORFs(E.coliismysequence,startCodon = "ATG", minimumLength = 100)
#ORFs by location
E.coliismysequence_orf <- GenomicRanges::sort(E.coliismysequence_orf)
#for tophits
tophit <- readDNAStringSet("tophit.fa")
tophit_orf <- ORFik::findORFs(tophit,startCodon = "ATG", minimumLength = 100)
tophit_orf <- GenomicRanges::sort(tophit_orf)
tophit_orf
#Extraction of sequences 
mystart = start(E.coliismysequence_orf)[[1]][1]
myend = end(E.coliismysequence_orf)[[1]][1]
#ORFs sequence  
ORF3a <- DNAStringSet(tophit,start = mystart[4], end = myend[4])
ORF3a <- toString(ORF3a)
ORF3a <- s2c(ORF3a)
ORF3a
#creating mutated copy 
ORF3a_mut <- mutator(myseq=ORF3a,100)
ORF3a_mut_ <- DNAString(c2s(ORF3a_mut))
#Conversion string to character
ORF3a_ <- DNAString(c2s(ORF3a))
#Pairwise alignment 
aln <- Biostrings::pairwiseAlignment(ORF3a_,ORF3a_mut_)
pid(aln)
#mismatched sequences
nmismatch(aln)



