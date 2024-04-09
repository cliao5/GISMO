library(dplyr)

prefix = "/Users/rye/Projects/GISMO/partitioned-heritability/GISMO/"
# Read original data (from Cal)
d <- read.csv2(paste0(prefix, "inputs/GISMOscore_26-07-2023.txt"), sep = "\t")
## Keep only relevant fields
d <- d[, c("gene", "dec")]
dim(d) # 19464 2

## Drop NA's
d <- d[complete.cases(d), ]
dim(d) # 19464 2

## Grab only first occurrences of gene names
d <- d[match(unique(d$gene), d$gene), ]
dim(d) # 19464 4

# Read map of gene names to ensembl gene names
t <- read.csv2("/Users/rye/Projects/cell-type-ldsc/reference/gencode/v44/gencode.v44.chr_patch_hapl_scaff.basic.annotation.gff3.genes", sep = "", header = T)

# Merge tables
merged <- merge(d, t, by.x = "gene", by.y = "gene_id", all.x = F)
dim(merged) # 29438 6

# Convert seqnames to int to get rid of non-autosome duplicates
merged$seqnames_fixed <- as.integer(merged$seqnames)
merged <- merged[complete.cases(merged$seqnames_fixed), ]
dim(merged) # 18548 7

# Remove duplicate of genes sharing genename
merged <- merged[match(unique(merged$gene_name), merged$gene_name), ]
dim(merged) # 18541 7

stopifnot(length(unique(merged$gene_name)) == dim(merged)[1])

# Count number of dropped
print(paste(sum(!(unique(d$gene) %in% merged$gene)), "dropped out of original", length(unique(d$gene)), "unique genes"))
#"923 dropped out of original 19464 unique genes


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
    write.table(m$gene, file = paste0(prefix, "inputs/decile/", gsub(" ", "", c, fixed = T), "/gene-set-file", sep = ""),
                quote = F, row.names = F, col.names = F)


    # Create --gene-coord-file
    ## a gene coordinate file, with columns GENE, CHR, START, and END, where START and END

    ### Pull needed info from merged
    coord <- m[, c("gene", "seqnames_fixed", "start", "end")]
    ### Fix column names
    colnames(coord) <- c("GENE", "CHR", "START", "END")
    ### Fix CHR by prepending "chr"
    coord$CHR <- paste("chr", coord$CHR, sep = "")
    write.table(coord, file = paste(prefix, "inputs/decile/", gsub(" ", "", c, fixed = T), "/gene-coord-file", sep = ""),
                quote = F, row.names = F, col.names = T, sep = "\t")
}

# Iterate through different random deciles
merged$rand_decile <- merged$decile[sample(length(merged$decile))]
for (c in unique(merged$rand_decile)) {
  # Create directory
  dir.create(file.path(paste0(prefix, "inputs/decile_rand/", c)), recursive = T, showWarnings = F)
  m <- merged[(merged$rand_decile == c), ]
  
  # Create --gene-set-file input
  ## a gene set file with the names of the genes in your gene set, one line per gene name
  write.table(m$gene, file = paste0(prefix, "inputs/decile_rand/", gsub(" ", "", c, fixed = T), "/gene-set-file", sep = ""),
              quote = F, row.names = F, col.names = F)
  
  
  # Create --gene-coord-file
  ## a gene coordinate file, with columns GENE, CHR, START, and END, where START and END
  
  ### Pull needed info from merged
  coord <- m[, c("gene", "seqnames_fixed", "start", "end")]
  ### Fix column names
  colnames(coord) <- c("GENE", "CHR", "START", "END")
  ### Fix CHR by prepending "chr"
  coord$CHR <- paste("chr", coord$CHR, sep = "")
  write.table(coord, file = paste(prefix, "inputs/decile_rand/", gsub(" ", "", c, fixed = T), "/gene-coord-file", sep = ""),
              quote = F, row.names = F, col.names = T, sep = "\t")
}


# Write list of deciles
writeLines(unique(as.character(merged$decile)), file(paste0(prefix, "inputs/decile-list")), sep = "\n")


