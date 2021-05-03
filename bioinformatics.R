#Install Bioconductor

#function testing if packages are installed
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
required_packages  <- c("BiocManager", "sangerseqR", "Biostrings")
lapply(required_packages,installBiocManager)

#load packages
library(tidyverse)
library(printr)
library(sangerseqR)
library(Biostrings)

#set working directory
setwd("~/Dropbox/R/bioinformatics")

#read sequence file
sequence <- readsangerseq("data/78_4212_91_EQ-CHF3_4239_SQ_F10_035.ab1")
#save the chromatogram to a pdf file
chromatogram(sequence, width = 100, height = 2, 
             filename = "data/sequence_chromatogram.pdf")
seq_text <- primarySeq(sequence, string = TRUE)
write_file(seq_text, "./data/seq_text.txt")
