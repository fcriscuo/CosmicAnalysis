install.packages("grid", "RColorBrewer","quantmod")
install.packages("seqinr")
source('C:/Rdev/CosmicAnalysis/readdata.R')
source('C:/Rdev/CosmicAnalysis/readdata.R')
dim(micro.df)
names(micro.df)
micro.df(1:2,)
micro.df[1;2,]
micro.df[1:2,]
source('C:/Rdev/CosmicAnalysis/readdata.R')
dim(cosmic.df)
install.packages("stringr")
?grep
?subset
names(cosmic.df)
cosmic.df[100:200,16]
grep("silent",cosmic.d[88,16])
grep("silent",cosmic.df[88,16])
grep("silxxx",cosmic.df[88,16])
x <- grep("silent",cosmic.df[88,16])
x
len(x)
length(x)
class(x)
y <- grep("silxxt",cosmic.df[88,16])
y
length(y)
x
x <- grep("silent",cosmic.df[88,16],fixed=TRUE)
x
class(x)
x == 1
test.df <- subset(cosmic.df, grep("silent",Mutation.Description")>0)
)
test.df <- subset(cosmic.df, grep("silent",Mutation.Description",fixed=TRUE)>0)
?subset
test.df <- subset(cosmic.df, grep("silent",Mutation.Description",fixed=TRUE)>0, select =names(cosmic.df))
names <- names.cosmic.df
names <- names(cosmic.df)
names
class9names
class(names)
cosmic.df[100:200,16]
sm <- "Substitution - coding silent"
test.df <- subset(cosmic.df, Mutation.Description == sm)
dim(test.df)
dim(cosmic.df)
names(test.df)
test.df[1;10,]
test.df[1:10,]
test.df[100:110,]
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
seq     <- getSequence(id = "A2M", type="hgnc_symbol", mart = ensembl, seqType = "transcript_exon_intron")
install.packages("RCurl")
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
seq     <- getSequence(id = "A2M", type="hgnc_symbol", mart = ensembl, seqType = "transcript_exon_intron")
seq
class(seq)
names(seq)
mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="refseq", mart = mart)
seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="peptide", mart = mart)
show(seq)
\
class(seq)
names(seq)
dim(seq)
seq = getSequence(id="BRCA1", type="hgnc_symbol", seqType="coding", mart = mart)
show(seq)
mart
datasets <- listDatasets(mart)
datasets
seq = getSequence(id="BRCA1", type="refseq", seqType="coding", mart = mart)
names(cosmic.df)
cosmic.df[,2]
seq = getSequence(id="ENST00000359318", type="ensembl", seqType="coding", mart = mart)
seq = getSequence(id="ENST00000359318", type="embl", seqType="coding", mart = mart)
show(seq)
x <-useMaart("CosmicMart")
x <-useMart("CosmicMart")
datasets <- listDatasets(x)
datasets
x <-useMart("CosmicMart",dataset="COSMIC66")
listAtrributes(x)
listFilters(x)
seq = getSequence(id="NG_023225.1", type="refseq", seqType="coding", mart = mart)
mart.listFilters()
class(mart)
listFilters(mart)
seq = getSequence(id="NG_023225.1", type="refseq_mrna", seqType="coding", mart = mart)
seq
seq = getSequence(id="NG_023225.1", seqType="coding", mart = mart)
<- org.Hs.egREFSEQ
# Get the entrez gene identifiers that are mapped to any RefSeq ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
# Get the REFSEQ for the first five genes
xx[1:5]
seq = getSequence(id="ENST00000380152", type="ensembl", seqType="coding", mart = mart)
seq = getSequence(id="ENST00000380152", type="ensembl", seqType="coding", mart = mart)
listFilters(mart)
seq = getSequence(id="BRCA1", type="refseq", seqType="coding", mart = mart)
library(seqinr)
?uco
seq <- "atggatttatctgctcttcgcgttgaagaagtacaaaatgtcattaatgctatgcagaaaatcttagagtgtcccatctgtctggagttgatcaaggaacctgtctccacaaagtgtgaccacatattttgcaaattttgcatgctgaaacttctcaaccagaagaaagggccttcacagtgtcctttatgtaagaatgatataaccaaaaggagcctacaagaaagtacgagatttagtcaacttgttgaagagctattgaaaatcatttgtgcttttcagcttgacacaggtttggagtatgcaaacagctataattttgcaaaaaaggaaaataactctcctgaacatctaaaagatgaagtttctatcatccaaagtatgggctacagaaaccgtgccaaaagacttctacagagtgaacccgaaaatccttccttgcaggaaaccagtctcagtgtccaactctctaaccttggaactgtgagaactctgaggacaaagca"
length(seq)
seq
result <- uco(seq)
class(result)
result
?read.fasta
all.genes <- read.fasta(file=".data/all_cosmic_genes.fasta", seqtype="DNA")
pwd()
cwd()
getwd()
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA")
class(all.genes)
all.genes[1:10]
result <- uco(seq,index="rscu")
result
names(all.genes)
all.genes[1]
rm(all.genes)
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA", as.string=TRUE)
names(all.genes)[1:10]
all.genes[1]
all.genes$ENSG00000197490
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA", as.string=TRUE,seqonly=TRUE)
all.genes[1]
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA", as.string=TRUE,aet.attributes=FALSE)
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA", as.string=TRUE,set.attributes=FALSE)
all.genes[1]
all.genes$BRCA1
seq <- all.genes$BRCA1
seq
seq[1]
class(seq)
result <- uco(seq,index=c("rscu"),frame=0)
result
result <- uco(seq,index=c("eff"),frame=0)
result
result <- uco(seq,index=c("freq"),frame=0)
result
seq
?grep
grep("atg",seq)
grep("ttt",seq)
words()
cds <- s2c(paste(words(), collapse = "")))
cds <- s2c(paste(words(), collapse = ""))
cds
uco(cds, index = "rscu")
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA", attributes=FALSE)
all.genes <- read.fasta(file="./data/all_cosmic_genes.fasta", seqtype="DNA", set.attributes=FALSE)
seq <- all.genes$BRCA1
seq
result <- uco(seq,index=c("rscu"),frame=0)
result
result <- uco(seq,index=c("freq"),frame=0)
result
result <- uco(seq,index=c("eff"),frame=0)
result
class(result)
result <- uco(seq,index=c("eff"),frame=0, as.data.frame=TRUE)
result
dim(test.df)
test.df[1]
test.df[1,]
seq
result
seq_g <-seq
seq[1353]
seq_g[1353]<-g
seq_g[1353]<-"g"
result_g <- uco(seq_g,index=c("eff"),frame=0, as.data.frame=TRUE)
result_g
1353%3
1353 %% 3
seq[1351:1353]
seq_g[1351:1353]
result
test.df[2,]
mut <- "c.1353A>G"
mut
str_extract(mut,"[0-9]")
library(stringr)
str_extract(mut,"[0-9]")
str_extract(mut,"[0-9]+")
as.integer(str_extract(mut,"[0-9]+"))
str_replace(mut,"c.","")
mut2 <- str_replace(mut,"c.","")
mut2
mut3 <-str_replace(mut2,">"," ")
mut3
mut3 <-str_replace(mut2,"[0-9]+,>","")
mut3
mut3 <-str_replace(mut2,"[0-9],"")
mut3 <-str_replace(mut2,"[0-9]","")
mut3
mut3 <-str_replace(mut2,"[0-9]+","")
mut3
mut3[1]
mut4 <- str_split(mut3,">")
mut4
mut4[1]
mut4[1,1]
mut4[1][1]
mut4[[1]]
mut4[[1]][1]
mut4[[1]][2]
savehistory("C:/Rdev/CosmicAnalysis/history.R")
