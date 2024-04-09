library(dplyr)

prefix = "partitioned-heritability/"
# output from GenerateScore
d <- read.csv2(paste0(prefix, "inputs/gismo-mis_consensus_genenames_01-08-2023.txt"), sep = "\t")
## Keep only relevant fields
d <- d[, c("genename", "genes", "dec", "mean.comb")]
d$mean.comb <- as.double(d$mean.comb)
dim(d) # 17439 4

## Drop NA's
d <- d[complete.cases(d), ]
dim(d) # 17438 4

## Grab only first occurrences of gene names, since we are only interested in converting gene names to ENSG
### Relaxed due to 100kb flankings
d <- d[match(unique(d$genename), d$genename), ]
dim(d) # 17433 4

# Read map of gene names to ensembl gene names
t <- read.csv2(paste0(prefix, "inputs/gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3.genes"), sep = "", header = T)

# Merge tables
merged <- merge(d, t, by.x = "genename", by.y = "gene_name", all.x = F)
dim(merged) # 19255 8

# Convert seqnames to int to get rid of non-autosome duplicates
merged$seqnames_fixed <- as.integer(merged$seqnames)
merged <- merged[complete.cases(merged$seqnames_fixed), ]
dim(merged) # 16595 9

# Remove duplicate of genes sharing genename
merged <- merged[match(unique(merged$genename), merged$genename), ]
dim(merged) # 16581 9

stopifnot(length(unique(merged$genename)) == dim(merged)[1])

# Count number of dropped
print(paste(sum(!(unique(d$genename) %in% merged$genename)), "dropped out of original", length(unique(d$genename)), "unique genes"))
#"852 dropped out of original 17433 unique genes"

# Compute deciles
merged$decile <- paste0("d", merged$dec)
print(table(merged$decile, useNA="always"))
stopifnot((is.na(merged$decile) == 0) & (merged$decile != ""))





# Iterate through different deciles
for (c in unique(merged$decile)) {
    # Create directory
    dir.create(file.path(paste0(prefix, "inputs/decile/", c)), recursive = T, showWarnings = F)
    m <- merged[(merged$decile == c), ]

    # Create --gene-set-file input
    ## a gene set file with the names of the genes in your gene set, one line per gene name
    write.table(m$gene_id, file = paste0(prefix, "inputs/decile/", gsub(" ", "", c, fixed = T), "/gene-set-file", sep = ""),
                quote = F, row.names = F, col.names = F)


    # Create --gene-coord-file
    ## a gene coordinate file, with columns GENE, CHR, START, and END, where START and END

    ### Pull needed info from merged
    coord <- m[, c("gene_id", "seqnames_fixed", "start", "end")]
    ### Fix column names
    colnames(coord) <- c("GENE", "CHR", "START", "END")
    ### Fix CHR by prepending "chr"
    coord$CHR <- paste("chr", coord$CHR, sep = "")
    write.table(coord, file = paste(prefix, "inputs/decile/", gsub(" ", "", c, fixed = T), "/gene-coord-file", sep = ""),
                quote = F, row.names = F, col.names = T, sep = "\t")
}

# Iterate through different random deciles
merged$rand_decile <- merged$decile[sample(length(merged$decile))]
print(table(merged$decile, useNA="always"))
print(table(merged$rand_decile, useNA="always"))
for (c in unique(merged$rand_decile)) {
  # Create directory
  dir.create(file.path(paste0(prefix, "inputs/decile_rand/", c)), recursive = T, showWarnings = F)
  m <- merged[(merged$rand_decile == c), ]
  
  # Create --gene-set-file input
  ## a gene set file with the names of the genes in your gene set, one line per gene name
  write.table(m$gene_id, file = paste0(prefix, "inputs/decile_rand/", gsub(" ", "", c, fixed = T), "/gene-set-file", sep = ""),
              quote = F, row.names = F, col.names = F)
  
  
  # Create --gene-coord-file
  ## a gene coordinate file, with columns GENE, CHR, START, and END, where START and END
  
  ### Pull needed info from merged
  coord <- m[, c("gene_id", "seqnames_fixed", "start", "end")]
  ### Fix column names
  colnames(coord) <- c("GENE", "CHR", "START", "END")
  ### Fix CHR by prepending "chr"
  coord$CHR <- paste("chr", coord$CHR, sep = "")
  write.table(coord, file = paste(prefix, "inputs/decile_rand/", gsub(" ", "", c, fixed = T), "/gene-coord-file", sep = ""),
              quote = F, row.names = F, col.names = T, sep = "\t")
}


# Write list of deciles
writeLines(paste0("d", c(1:10)), file(paste0(prefix, "inputs/decile-list")), sep = "\n")


