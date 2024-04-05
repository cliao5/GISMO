# Assign column species names to unlabeled matrix to generate final one2[x] matrix (combined_matrix.tsv)
d <- read.table("unlabelled_matrix", header = F, sep = "", row.names = 1)
names <- read.table("names", header = F, sep = "")

colnames(d) <- names$V1
# Output moved to GenerateScore/inputs/unmerged-species_combined_matrix_2023-07-26.tsv
write.table(d, file = "combined_matrix.tsv", sep = "\t", quote = F)
