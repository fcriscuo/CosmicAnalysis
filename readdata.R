setwd("C:/Rdev/CosmicAnalysis")
r <- getOption("repos")             # hard code the US repo for CRAN
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)
rm(r)

if ("seqinr" %in% rownames(installed.packages()) == FALSE){
  install.packages("seqinr",dep=TRUE)
}

if ("xlsx" %in% rownames(installed.packages()) == FALSE){
  install.packages("xlsx",dep=TRUE)
}

if ("stringr" %in% rownames(installed.packages()) == FALSE){
  install.packages("stringr",dep=TRUE)
}

if ("sqldf" %in% rownames(installed.packages()) == FALSE){
  install.packages("sqldf",dep=TRUE)
}

source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("biomaRt")
library(seqinr)
library(stringr)
library(xlsx)

#process large mutation file on 2000 row chunks 
# filter out synonmous mutations
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
            test.row.limit, nrow(silent.df)))
 #write out silent data frame to excel file
    fileName <- "./data/cosmic_silent_mutations.xlsx"
    write.xlsx(silent.df,fileName)
    
    