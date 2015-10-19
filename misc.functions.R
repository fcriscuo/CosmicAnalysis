#represents a collection of utility functions used throughout the
# application. 
# working directory determined by function callers
#
require(seqinr)
require(Biostrings)
require(stringr)
require(xlsx)

#read in trna frequency data
filename <- "./data/trna_freq.xlsx"
sheet <- 1
trna.df <- read.xlsx(filename,sheet, 
                     stringsAsFactors=FALSE )

#function to return the aa.position as a numeric value from
# a COSMIC aa mutation location (e.g. p.A3663A) or a COSMIC CDS 
# mutation ( e.g. c.1089C>A)
mut.position <- function(cosmic.mutation) {
  s1 <- str_extract(cosmic.mutation,"[0-9]+")
  as.numeric(s1)
}

#function to determine codon from a sequence
codon <- function(seq, aa.mutation){
  # TODO: need to confirm DNA sequence
  # TODO: need to confirm aa mutation
  codon.number <- mut.position(aa.mutation)
  codons <- splitseq(s2c(seq))
  codons[codon.number] 
}

#function to find the mutated aa or nucleotide from COSMIC
# notation

parse.mutation <- function(mutation){
   #x <- str_replace(mutation,"[a-z]+.[0-9]+","") 
   x <- str_replace(mutation,"[a-z]+.","") 
   x <- str_replace(x,"[0-9]+","") 
   x  <- str_replace(x,">","")
   substring(x,2)
  
}



#determine the mutated codon using the wild type sequence
# and the cds mutation
mutated.codon <- function(wt.seq, cds.mutation, aa.mutation) {
  #TODO: validate inputs
  #replace the wild type nucleotide
  # with the cds mutation, the return the codon
  cds.mut.position <- mut.position(cds.mutation)
  cds.mut <- parse.mutation(cds.mutation)
  substr(wt.seq,cds.mut.position,cds.mut.position) <- cds.mut
  #return the mutated codon
  toupper(codon(wt.seq,aa.mutation))
 
}

#function to determine the change in codon frquency as a result
# of the silent mutation
delta.freq <- function(wt.seq, cds.mutation, aa.mutation ){
  #TODO: validate the function arguments
  wt.codon <- toupper(codon(wt.seq, aa.mutation))
  wt.freq <- trna.df$freq[trna.df$codon==wt.codon]
  mut.codon <- mutated.codon(wt.seq, cds.mutation, aa.mutation)
  mut.freq <- trna.df$freq[trna.df$codon==mut.codon]
  mut.freq-wt.freq
  
}