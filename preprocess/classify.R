library(DECIPHER)

# [TODO] fix to correct path

# fasta file with the query sequences
fas <-"seqs_tmp.fa"

# load the sequences from the file
seqs <- readDNAStringSet(fas)

# remove any gaps (if needed)
seqs <- RemoveGaps(seqs)

# load training set object (trainingSet)
load("gtdb_16s.RData")

# classify the sequences
ids <- IdTaxa(seqs,
              trainingSet,
              strand="both", # or "top" if same as trainingSet
              threshold=60, # 60 (cautious) or 50 (sensible)
              processors=NULL) # use all available processors

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  merged_taxon <- paste(x$taxon, collapse = ";")
  merged_taxon
}))


taxid_df <- data.frame(
  seqid = colnames(taxid),
  path = as.vector(taxid)
)

# Specify the file path
file_path <- "assignments.tsv"

# Write the data frame to a tab-separated file
write.table(taxid_df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
