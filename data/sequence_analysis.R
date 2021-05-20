#install latest version of Bioconductor by entering the commands below (the lines starting with "#>"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#install required packages
BiocManager::install(c("Biostrings", "sangerseqR",
                       "annotate", "msa", "ggmsa"))

#load packages
library(sangerseqR)
library(Biostrings)
library(annotate)
library(msa)
library(ggmsa)
library(ggplot2)


#read the ab1 sequence file
sequence <- readsangerseq("./data/78_4212_91_EQ-CHF3_4239_SQ_F10_035.ab1")

#call the bases from the data provided in the file and display the result on the screen
#type "makeBaseCalls into the search field of the "Help" pane on the lower left of your workspace and read about what this funtion does
makeBaseCalls(sequence, ratio = 0.33)

#save the chromatogram to a pdf file
chromatogram(sequence,0,0, width = 100, height = 2, filename = "./data/sequence_chromatogram.pdf")

#save the sequence as text
seq_text <- primarySeq(sequence, string = TRUE)

#write the sequence to file in FASTA format
seq_fasta <- DNAStringSet(seq_text)
writeXStringSet(seq_fasta, file="./data/seq.fasta", format="fasta", width=80)

#load the file that contains the query sequence for the BLAST
query <- readDNAStringSet("data/seq.fasta")

#BLAST the sequence
blast_output <- blastSequences(query, as="data.frame")

#Inspect the result
detail(blast_output)
print(blast_output, show = "complete")

#save the returned matches as a DNA string set class
matching_fasta <- DNAStringSet(c(blast_output$Hsp_qseq[1], blast_output$Hsp_hseq[c(1,3,5,7,9)]))

#add the names of the returned sequences
names(matching_fasta) <- c(rep("query",1), blast_output$Hit_id[c(1,3,5,7,9)])

#write the alignment to a file
writeXStringSet(matching_fasta, file="./data/Equisetum_blast.fasta", format="fasta", width=80)

#plot the alignment
ggmsa(matching_fasta, color = "Chemistry_NT", font = "DroidSansMono", start = 1, end = 600,
      char_width = 0.5, seq_name = TRUE)  +
  facet_msa(field = 80)
ggsave("./data/Equisetum_alignment.png", device = "png", dpi = "screen", width = 8.5, height = 11, units = "in")