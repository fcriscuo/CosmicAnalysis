
# script to read fata file containing all human genes in cosmic
# database
# returns a data frame with the gene name and reference sequence
setwd("C:/Rdev/CosmicAnalysis")
all.genes <- read.fasta(file="data/all_cosmic_genes.fasta", 
           seqtype = c("DNA"), as.string = TRUE)
ng <- length(all.genes)
all.genes.df <- data.frame(name=rep("",ng), seq=rep("",ng)
                          ,stringsAsFactors=FALSE)
for(i in seq(ng)){
  g <- all.genes[[i]]
  all.genes.df$seq[i] <- as.character(g)
  all.genes.df$name[i] <- attributes(g)$name
}
 
#all.genes.df
fileName <- "./data/cosmic_all_genes.xlsx"
write.xlsx(all.genes.df,fileName)
