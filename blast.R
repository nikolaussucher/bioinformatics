# BLAST search in R
BiocManager::install("annotate")
BiocManager::install("msa")
BiocManager::install("seqinr")
install.packages("ggmsa")
install.packages("bio3d")
library(annotate)
library(Biostrings)
library(msa)
library(seqinr)
library(msa)
library(ggmsa)
library(alignfigR)
library(bio3d)

blast_output <- blastSequences("CTTCTTTTTTACATTTAACAAGATTTTTTNTACATGAATATAAAGAAAGGTCTAAAAATGAAAAATTGATTTATTCTTTTAATAAAAGAAATAAATTTGTCACATTTTCTTGGAATTTTTTTTATATTTTTGAATTGGAATTTTTTTTAACTTCTCTTTTGACTAGGTTTNTAAATAATTTAGTTCTATCCAACTTATTCGATCAAATTAATCTATTAGAAAAAATAAATAATAATAATAAAAGTCAATATTTTTTATCAGAAAAAAGGATCTATAACCAAAATTCCTGTATTCATTATGTACGGTATCAAAATCGTTGTATTATGGCTTCGGAAGGCTTTTATTTTCATGATACGAACTGGATATATTTTATATTGAATATTTGGCAATTTTTTATGCATTTATGGATTCAACCTTTCAGGTTCTCAACAAAGCATTTTCAAAAGCAGAGTTTTTTTTTTCTGGGTTATCAATTTGGCAGAGAATCCAAACTGCTAAAAGTTCGATCTATTTCACTTGATAAATCACCTACAATTTATTCACGTTTAAAAAAAAATCTTTTAAAAACACAAATTGTATATCCAATAGATTTTCTAGCAAAAGAGGGTTTTTGTGATATTTCAGGTTATCCTATTAGTAGATCGACTTGGACTACNTCGACNGATGAAGAANTTTTATTGAATTTTAATAAAATTTGNAAAAGCTTTTATTTTTATTA", as="data.frame")

fasta <- DNAStringSet(blast_output$Hsp_hseq)
names(fasta) = as.character(blast_output$Hit_id)
writeXStringSet (fasta, file="Equisetum.fasta", format="fasta", width=80)
fasta$`gi|1988363571|gb|MT197578.1|`
names(fasta$`gi|1988363571|gb|MT197578.1|`)
head(blast_output)
blast_output$Hit_def[1]
blast_output$Iteration_hits

blast_output$Hsp_qseq

query_seq_fasta <- DNAStringSet(blast_output$Hsp_qseq)
writeXStringSet (query_seq_fasta, file="My_Equisetum.fasta", format="fasta", width=80)


read.alignment("Equisetum.fasta","fasta")
mySequenceFile <- readDNAStringSet("Equisetum.fasta")
myAlignment <- msa(mySequenceFile, type="dna")

print(myAlignment, show = "complete")

aln <- read.fasta("Equisetum.fasta")
plot(aln, labels=basename.pdb(aln$id))
aln2html(aln, file = "align.html")

aln <- Biostrings::readDNAMultipleAlignment(aln)

Biostrings::readDNAMultipleAlignment(mySequenceFile)

ggmsa(myAlignment, color = "Chemistry_NT", start = 50, end = 100, char_width = 0.5) + 
   facet_msa(field = 50)
  
  
  geom_seqlogo(color = "Chemistry_NT") + 
  geom_GC() + 
  theme(legend.position = "none")
  
  class(myAlignment)[[1]] <- "DNAMultipleAlignment""
  factors(myAlignment)
