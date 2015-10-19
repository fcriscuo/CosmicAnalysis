#R script to calculate human tRNA frequency based on the number of tRNA genes
# in the tRNA database
setwd("C:/Rdev/CosmicAnalysis")
r <- getOption("repos")             # hard code the US repo for CRAN
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)
rm(r)

library(seqinr)
library(stringr)
library(xlsx)
#biocLite("Biostrings")
library(Biostrings)
# read in the human tRNA genes

trna.df <- read.delim("./data/trna_freq.csv", sep=",", header=TRUE, 
                       stringsAsFactors=FALSE)

aa.list <- unique(trna.df$aa)
nn <- length(unique(trna.df$codon)) 

freq.df <- data.frame(aa=rep("", nn), total=rep(0,nn), codon=rep("", nn),
                      genes=rep(0,nn),   freq=rep(0.0,nn),
                      rate=rep(0.0,nn)
                      ,stringsAsFactors=FALSE )
max.freq.df <- data.frame(aa=rep("",length(aa.list)), 
                          max.freq = rep(0.0,length(aa.list)),stringsAsFactors=FALSE)

#the codons in the database are the anti-codons used for tRNA base
# pairing to mRNA 
# we need to match the genome codon - persist a reverse compliment
revTranscribe <- function(codon){
  c <- DNAString(codon)
  rc <- reverseComplement(c)
  as.character(rc)
}
index <-1
for(i in seq(length(aa.list))){
  AA <- aa.list[i]
  max.freq.df$aa[i] <- aa.list[i]
  aa.df <- subset(trna.df, aa == AA)
  ncodon = sum(aa.df$genes)
  codon.list <- unique(aa.df$codon)
  for(j in seq(length(codon.list))){
    CODON <- codon.list[j]
    codon.df <- subset(aa.df, codon == CODON)
    ngene = codon.df$genes
    freq.df$aa[index] <- AA
    freq.df$total[index] <- freq.df$total[index] + ngene
    freq.df$codon[index] <-revTranscribe(CODON)
    freq.df$genes[index] <- ngene
    freq <- (ngene / ncodon)
    freq.df$freq[index]<-  freq
    if ( freq > max.freq.df$max.freq[i]){
      max.freq.df$max.freq[i] <- freq
    }
      
    cat(sprintf(" \nAA  %s  total %d codon %s genes %d freq %f", AA, ncodon, 
                revTranscribe(CODON), ngene, ngene/ncodon))
    index <- index +1
  }
  
  
}

# calculate the relative rate of each codon within a aa group


for ( i in seq(nrow(max.freq.df))){
  aa <- max.freq.df$aa[i]
  max <- max.freq.df$max.freq[i]
  #look up the codons for this aa
  codons <- freq.df[freq.df$aa ==aa,]$codon
  #set the rate
  for (j in seq(length(codons))) {
    freq.df[freq.df$codon ==codons[j],]$rate <-
      freq.df[freq.df$codon ==codons[j],]$freq /max
  }
  
}

#export to excel file
fileName <- "./data/trna_freq.xlsx"
write.xlsx(freq.df,fileName)

