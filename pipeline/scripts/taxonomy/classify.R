#!/usr/bin/env Rscript
#
# Taxonomy assignment of amplicon sequences with IDTAXA (DECIPHER), against a reference (a DECIPHER "training set") built
# from the 16S rRNA gene sequences of the GTDB genomes.
#
# Usage: classify.R <sequences.fasta> <reference.RData> <threshold> <assignments.tsv> [seed]
#
# Output: a tab-separated file with the sequence id and its taxonomy (ranks separated by ';').
#
# Same data, different results? IdTaxa is randomized: for each sequence it repeatedly samples a random subset of
# the sequence's k-mers (bootstrap replicates) and checks which taxon each subset points to; the confidence of a
# rank is the share of replicates that agree on it, and ranks below the threshold are left out (unclassified_...).
# For a sequence whose confidence is close to the threshold, another random sample can put it on the other side,
# so it is assigned one rank deeper or shallower. Without a fixed seed, ~10% of the assignments of the test data
# (99 of 996 ASVs) changed between two runs on the same input. With the seed (default 1), a run gives the same
# assignments as another with the same input, reference, threshold, seed and DECIPHER version (checked with the
# same number of processors). Another seed, or an output made without one (e.g. the former microbetag_prep image),
# differs in the same way: in the depth of the assignments near the threshold, not in the lineages.

suppressPackageStartupMessages(library(DECIPHER))

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(4, 5))) {
  stop("Usage: classify.R <sequences.fasta> <reference.RData> <threshold> <assignments.tsv> [seed]")
}
fasta_file    <- args[1]
training_file <- args[2]
threshold     <- as.numeric(args[3])
out_file      <- args[4]
seed          <- if (length(args) == 5) as.integer(args[5]) else 1L

seqs <- RemoveGaps(readDNAStringSet(fasta_file))

# Loads the `trainingSet` object
load(training_file)

# IdTaxa samples k-mers at random (see above): the seed makes the assignments reproducible
set.seed(seed)

ids <- IdTaxa(seqs,
              trainingSet,
              strand     = "both",     # or "top" if the sequences have the orientation of the reference
              threshold  = threshold,  # 60 (cautious) or 50 (sensible)
              processors = NULL)       # all available processors

assignments <- data.frame(
  seqid    = names(ids),
  taxonomy = vapply(ids, function(x) paste(x$taxon, collapse = ";"), character(1))
)

write.table(assignments, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
