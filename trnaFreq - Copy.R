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
all.trna <- read.fasta(file="./data/human-trnas.fa", seqtype="DNA", as.string=TRUE)
n <- length(all.trna)
#create a data frame to hold these data

temp.df <- data.frame(aa=rep("", n),codon=rep("", n),
                      desc=rep("",n),
                      stringsAsFactors=FALSE)
for ( i in seq(n)){
  #get the entry's Annot attribute
  annot <- attributes(all.trna[[i]])$Annot
  temp.df$aa[i] <- word(annot,4)
  
  temp.df$codon[i] <- str_match(word(annot,5),"[A-Z]+")
  temp.df$desc[i] <- word(annot,1)
}

 #filter out duplicates
 trna.df <- unique(temp.df)
#filter out undetermined and Sec
  trna.df <- subset(trna.df, aa != "Undet")
  trna.df <- subset(trna.df, aa != "SeC")
  trna.df <- subset(trna.df, aa != "SeC(e)")

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
  ncodon = nrow(aa.df)
  codon.list <- unique(aa.df$codon)
  for(j in seq(length(codon.list))){
    CODON <- codon.list[j]
    codon.df <- subset(aa.df, codon == CODON)
    ngene = nrow(codon.df)
    freq.df$aa[index] <- AA
    freq.df$total[index] <- ncodon
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

