#!/usr/bin/env Rscript
#
# Taxonomy assignment of amplicon sequences with IDTAXA (DECIPHER), against a training set built
# from the 16S rRNA gene sequences of the GTDB genomes.
#
# Usage: classify.R <sequences.fasta> <training_set.RData> <threshold> <assignments.tsv>
#
# Output: a tab-separated file with the sequence id and its taxonomy (ranks separated by ';').

suppressPackageStartupMessages(library(DECIPHER))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: classify.R <sequences.fasta> <training_set.RData> <threshold> <assignments.tsv>")
}
fasta_file    <- args[1]
training_file <- args[2]
threshold     <- as.numeric(args[3])
out_file      <- args[4]

seqs <- RemoveGaps(readDNAStringSet(fasta_file))

# Loads the `trainingSet` object
load(training_file)

ids <- IdTaxa(seqs,
              trainingSet,
              strand     = "both",     # or "top" if the sequences have the orientation of the training set
              threshold  = threshold,  # 60 (cautious) or 50 (sensible)
              processors = NULL)       # all available processors

assignments <- data.frame(
  seqid    = names(ids),
  taxonomy = vapply(ids, function(x) paste(x$taxon, collapse = ";"), character(1))
)

write.table(assignments, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
