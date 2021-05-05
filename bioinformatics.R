#Install Bioconductor

#function testing if Bioconductor packages are installed
installBiocManager <- function(x){
if(x %in% rownames(installed.packages())==FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(x,ask = FALSE)
}  
  else{
    paste(x,"is already installed!")
    }
}

#install necessary Bioconductor packages
required_packages  <- c("BiocManager", "sangerseqR", "Biostrings", "ORFik")
lapply(required_packages,installBiocManager)

#function testing if R packages are installed
installRpackages <- function(x){
  if(x %in% rownames(installed.packages())==FALSE) {
      install.packages(x)
  }  
  else{
    paste(x,"is already installed!")
  }
}

required_packages  <- c("seqinr", "htmltools", "files")
lapply(required_packages,installRpackages)

#load packages
library(tidyverse)
library(printr)
library(sangerseqR)
library(Biostrings)
library(seqinr)
library(htmltools)
library(files)

#set working directory
setwd("~/Dropbox/R/bioinformatics")

#read sequence file
sequence <- readsangerseq("data/78_4212_91_EQ-CHF3_4239_SQ_F10_035.ab1")
#save the chromatogram to a pdf file
chromatogram(sequence, width = 100, height = 2, 
             filename = "data/sequence_chromatogram.pdf")
#save the sequence as text
seq_text <- primarySeq(sequence, string = TRUE)
#write the sequence to a tex file
write_file(seq_text, "./data/seq_text.txt")
write.fasta(seq_text, names="sequence", file.out ="./data/seq.fasta", open = "w", nbchar = 40, as.string = TRUE)
s <- read.fasta("./data/seq.fasta",seqonly = TRUE)
extractseqs(s)
print(s)
findORFs(seq_txt)

