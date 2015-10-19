install.packages("grid", "RColorBrewer","quantmod")
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
source('C:/Rdev/CosmicAnalysis/readdata.R')
source('C:/Rdev/CosmicAnalysis/readdata.R')
source('C:/Rdev/CosmicAnalysis/readdata.R')
source('C:/Rdev/CosmicAnalysis/readdata.R')
nrow <- 2000
skip.rows <- 0
sm <- "Substitution - coding silent"
block.df <- read.delim("./data/cosmicmutantexport.tsv", sep="\t", header=TRUE,
stringsAsFactors=FALSE,nrows=nrow,skip=skip.rows)
silent.df <- subset(block.df, Mutation.Description == sm)
skip.rows <- skip.rows + nrow
names <- names(block.df)
test.row.limit <-20000
while(skip.rows <test.row.limit) {
cat(sprintf("\nStarting wirh row %d", skip.rows))
block.df <- read.delim("./data/cosmicmutantexport.tsv", sep="\t", header=FALSE,                        stringsAsFactors=FALSE,nrows=nrow,skip=skip.rows)
names(block.df)<-names
subset.df <- subset(block.df, Mutation.Description == sm)
silent.df <- rbind(silent.df, subset.df)
skip.rows <-skip.rows + nrow
}
cat(sprintf("\n%d rows have been read there are %d silent mutations",
test.limit.rows, nrow(silent.df)))
#write out silent data frame to excel file
fileName <- "./data/cosmic_silent_mutations.xlsx"
write.xlsx(silent.df,fileName)
cat(sprintf("\n%d rows have been read there are %d silent mutations",
test.row.limit, nrow(silent.df)))
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
names(trna.df)
names(freq.df)
freq.df[1:10]
freq.df[1:10,]
silent.df[100,]
3216 %% 3
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
3216 %% 3
?readFasta
readfasta
?readfasta
library(seqinr)
library(stringr)
library(xlsx)
all.genes <- read.fasta(file = system.file("data/all_cosmic_genes.fasta", package = "seqinr"),
seqtype = c("DNA"), as.string = FALSE
)
all.genes <- read.fasta(file = system.file("data/all_cosmic_genes.fasta", package = "seqinr"),
seqtype = c("DNA"), as.string = FALSE)
all.genes <- read.fasta(file="data/all_cosmic_genes.fasta",
seqtype = c("DNA"), as.string = FALSE)
names(all.genes)
class(all.genes)
length(all.genes)
all.genes["A2ML1"]
all.genes["A2ML1"][3216]
x<-all.genes["A2ML1"]
x[[3216]]
x
x[[1]]
y <<-x[[1]]
class(y)
?SeqFastadna
y$object
y
length(y)
y[3216]
y[3216:3218]
y[3214:3216]
biocLite("Biostrings")
x <- DNAString("ACGT-YN-")
reverseComplement(x)
library(Biostrings)
x <- DNAString("ACGT-YN-")
reverseComplement(x)
x <- DNAString("GCT")
y <- reverseComplement(x)
y
class(y)
y[[1]]
?DNAString
length(d)
length(y)
y
attributres(y).seq
attributres(y)$seq
attributes(y)$seq
y
y$seq
alphabet(y)
y
seq(y)
subseq(y)
y@seq
y
slotNames(y)
as.character(y)
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
class(all.genes)
x<-all.genes["A2ML1"]
x
y <- x[[1]]
y
y[3214:3216]
y
k <- as.character(y)
k
?uco
uco.df <- uco(k, as.data.frame=TRUE)
uco.df
?DNAString
y
k
k[3216]
k[3216] <-"g"
uco_g.df <- uco(k, as.data.frame=TRUE)
uco_df
uco.df
uco.df[uco.df$AA=='Pro',]
uco_g.df[uco.df$AA=='Pro',]
frwq.df[aa=="pro",]
freq.df[aa=="pro",]
names(freq.df)
freq.df[freq.df$aa=="pro",]
freq.df[freq.df$aa=="Pro",]
freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
0.174/ 0.282
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
warnings()
max.freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
warnings()
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
max.fre.df
max.freg.df
max.freq.df
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
x <- "Leu"
freq.df[freq.df$aa ==x,]
freq.df[freq.df$aa ==x,]$codon
l1 <- freq.df[freq.df$aa ==x,]$codon
li[1]
l1[1]
l2 <-  freq.df[freq.df$codon ==l1[1],]
l2
freq.df[freq.df$codon ==l1[1],]$rate <- 1.00
freq.df[freq.df$codon ==l1[1],]
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
freq.df
x <- "asdfgtyuu"
strsplit(x)
strsplit(x,split="")
all.genes <- read.fasta(file="data/all_cosmic_genes.fasta",
seqtype = c("DNA"), as.string = TRUE)
length(all.genes)
all.genes[1]
y <- all.genes[[1]]
y
as.character(y)
y
attributes(y)$name
ng <- length(all.genes)
all.genes.df <- data.frame(name=rep("",ng), seq=rep("",ng)
,stringsAsFactors=FALSE)
for(i in seq(ng)){
g <- all.genes[[i]]
all.genes.df$seq[i] <- as.character(g)
all.genes.df$name[i] <- attributes(g)$name
}
dim(all.genes.df)
all.genes.df[100:105,]
x <- all.genes.df$seq[1000]
x
strsplit(x,strsplit="")
strsplit(x,split="")
?read.xlsx
fileName <- "./data/cosmic_all_genes.xlsx"
write.xlsx(all.genes.df,fileName)
setwd("C:/Rdev/CosmicAnalysis")
write.xlsx(all.genes.df,fileName)
source('C:/Rdev/CosmicAnalysis/trnaFreq.R')
savehistory("C:/Rdev/CosmicAnalysis/history2.R")
