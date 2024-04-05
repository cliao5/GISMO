#Finalized GISMO generation metric
library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)
library("viridis")
library(ape)# Load
library(phytools)
library(vctrs)

gene_stat <- data.table::fread('inputs/gen.dist.withassembly.txt', data.table = F)  
gene_stat <- gene_stat[order(gene_stat$species, -abs(gene_stat$contig_N50_bp) ), ] #sort by id and drop lower quality genomes for duplicates
gene_stat <- gene_stat[!duplicated(gene_stat$species),]

#check full list
o2o <- read.table('inputs/unmerged-species_combined_matrix_2023-07-26.tsv', header=T, stringsAsF=F)
colnames(o2o) <- gsub(".*_","",colnames(o2o)) 
o2o[!o2o == "one2zero"] <- 0
o2o[o2o == "one2zero"] <- 1
o2o <- mutate_all(o2o, function(x) as.numeric(as.character(x)))
species_use <- intersect(gene_stat$assembly_name, colnames(o2o))

deciles <- max(gene_stat$lin.num)
gene_stat <- gene_stat[gene_stat$assembly_name %in% species_use,]
gene_stat <- gene_stat[order(gene_stat$lin.num),]
#gene_stat$dec <- dplyr::ntile(gene_stat$lin.num, n=deciles) #don't need deciles unless evo.time 

o2o <- o2o[,gene_stat$assembly_name]
o2o <- t(apply(o2o, 1, as.vector))
dim(o2o) # 19464   462

#jfu 
#deciles <- 10
#gene_stat <- gene_stat[gene_stat$species %in% species_use,]
#gene_stat <- gene_stat[order(gene_stat$div.time),]
#gene_stat$dec <- dplyr::ntile(gene_stat$div.time, n=deciles)

#o2o <- o2o[,gene_stat$species]
#o2o <- t(apply(o2o, 1, as.vector))
#dim(o2o) # 19204 genes, 462 species
#vec_decile <- gene_stat$dec

### Estimate p_hat (each gene by decile)
vec_decile <- gene_stat$lin.num #divided ******8
mat_decile <- matrix(rep(vec_decile, each=dim(o2o)[1]), nrow=dim(o2o)[1], byrow=FALSE)
mat_decile2 <- mat_decile + matrix(rep(((1:dim(o2o)[1])-1)*deciles, each=dim(o2o)[2]), ncol=dim(o2o)[2], byrow=TRUE)

p_hat <- by(as.vector(o2o), as.vector(mat_decile2), mean)
mat_p_hat <- matrix(p_hat[match(mat_decile2, names(p_hat))], ncol=dim(o2o)[2], byrow=FALSE)

### Generate weights matrix
#(1/gene_stat$weights)*(gene_stat$div.time)/(max(gene_stat$div.time))*gene_stat$order.weight
vec_weights <- gene_stat$div.time/max(gene_stat$div.time)
#vec_weights <-  1-((gene_stat$sum.branch.dist)/(max(gene_stat$sum.branch.dist))) #(gene_stat$div.time^2)/(gene_stat$div.time^2) 
mat_weights <- matrix(rep(vec_weights, each=dim(o2o)[1]), ncol=dim(o2o)[2], byrow=FALSE)

### Generate observed stat
stat_obs <- apply(o2o*mat_weights, 1, mean)
stat_obs2 <- as.data.frame(stat_obs)
stat_obs2$gene <- rownames(stat_obs2)

### Generate RV
set.seed(1337)
iter <- 10000

stat_rv_iter <- matrix(NA, nrow=dim(o2o)[1], ncol=iter)
for(i in 1:iter){
  message(i)
  rv <- matrix(rbinom(length(mat_p_hat), 1, as.vector(mat_p_hat)), byrow=FALSE, ncol=dim(o2o)[2])
  # plot(apply(rv, 1, mean)~apply(o2o, 1, mean)) ### Check that the dimensions are correct
  #stat_rv_iter[,i] <- rowSums(rv*mat_weights)/mu2
  stat_rv_iter[,i] <- apply(rv*mat_weights, 1, mean)  
}
stat_ci <- t(apply(stat_rv_iter, 1, function(x) quantile(x, c(0.025,0.5, 0.975))))

df <- as.data.frame(stat_ci)
df$gene <- rownames(o2o)
colnames(df)[1:3] <- c("lower.ci","median","upper.ci")
#write.table(df, "GISMOscore_26-07-2023.txt", quote=F, row.names=F, sep="\t")
#df <- read.table("one2zero.upperci_20quantile.1000perm_divtime2_2023-03-09.txt", header=T)
df$dec <- ntile(df$upper.ci,10)
df$dec_median <- ntile(df$median,10)


## Generate GISMO-mis scores ##
library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)
library("viridis")   

#ortholog missense test
genedist <- data.table::fread('inputs/genestat_26-07-2023.txt', data.table = F)  
x <- fread("inputs/consensus_v2_missense-counts.gz", stringsAsFactors = F, data.table=F)
rownames(x) <- x$species
x <- x[x$species %in% gene_stat$assembly_name,]
x <- x[order(x$species),]
syn <- fread("inputs/consensus_v2_synonymous-counts.gz", stringsAsFactors = F,data.table=F)
syn <- syn[syn$species %in% gene_stat$assembly_name,]
syn <- syn[order(syn$species),]
x <- x[,colnames(x) %in% colnames(syn)]
syn.tmp <- syn[which(rownames(syn) %in% rownames(x)),]

#syn2 <- fread("~/bd/ortholog/all_synonymous-counts.gz")
syn <- as.data.frame(syn)
syn[syn == "None"] <- NA
rownames(syn) <- syn$species
syn <- mutate_all(syn[,-1], function(x) as.numeric(as.character(x)))
syn.means <- as.data.frame(colMeans(syn, na.rm = T))
syn.means$genes <- rownames(syn.means)
colnames(syn.means)[1] <- "mean.syn"

x[x == "None"] <- NA
df2 <- mutate_all(x[,-1], function(x) as.numeric(as.character(x)))
rownames(df2) <- rownames(x)
means <- as.data.frame(colMeans(df2, na.rm = T))
means$genes <- rownames(means)
colnames(means)[1] <- "mean.mis"

#sort df2, syn, and gene_stat
genedist <- genedist[genedist$assembly_name %in% rownames(df2),]
df2 <- df2[order(rownames(df2)),]
syn <- syn[order(rownames(syn)),]
genedist <- genedist[order(genedist$assembly_name),]

#Calculate GISMO-mis 
combined <- (df2)/(syn)
combined[combined == Inf] <- NA

comb.means <- as.data.frame(colMeans(combined, na.rm = T))
comb.means$genes <- rownames(comb.means)
colnames(comb.means)[1] <- "mean.comb"

#mis/syn
combined.means <- merge(means, syn.means, by="genes")
combined.means$mis.cal <- combined.means$mean.mis/combined.means$mean.syn
combined.means <- merge(combined.means, comb.means, by="genes")

#write.table(combined.means, "gismo-mis_consensus_01-08-2023.txt",quote=F, row.names=F, sep="\t")



#write.table(merged, "gismo-mis.sneath_2023-06-25.txt",quote=F, row.names=F, sep="\t")


