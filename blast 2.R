#BLAST search in R

#set working directory
setwd("~/Dropbox/R/bioinformatics")

#function testing if Bioconductor packages are installed
installBiocManager <- function(x)
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
required_packages  <- c("BiocManager", "sangerseqR", "Biostrings", "annotate",
                       "msa", "ggmsa")
lapply(required_packages,installBiocManager)

#load packages
library(Biostrings)
library(annotate)
library(msa)
library(ggmsa)
library(ggplot2)

#read sequence file to BLAST

query <- readDNAStringSet("data/seq.fasta")

blast_output <- blastSequences(query, as="data.frame")
fasta <- DNAStringSet(c(blast_output$Hsp_qseq[1], blast_output$Hsp_hseq[c(1,3,5,7,9)]))
names(fasta) <- c(rep("query",1), blast_output$Hit_id[c(1,3,5,7,9)])

fasta <- readDNAStringSet("Equisetum_blast.fasta")

aln <- msa(fasta)

DNAStringSet(aln)

writeXStringSet(fasta, file="Equisetum_blast.fasta", format="fasta", width=80)

ggmsa(fasta, color = "Chemistry_NT", font = "DroidSansMono", start = 1, end = 400,
      char_width = 0.5, seq_name = TRUE)  + 
      facet_msa(field = 80)

ggsave("Equisetum_alignment.pdf", dpi = "print", width = 8.5, height = 11, units = "in")

