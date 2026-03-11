### thermo_capacity_rnaseq_wgcna_only.R
# This script runs the WGCNA analyses presented in Velotta et al. Proc Roy Soc B
# It is written to be run from within an RStudio Project (.Rproj).
#
# folder structure:
# thermogenic_capacity_rnaseq/
#   thermogenic_capacity_rnaseq.Rproj
#   scripts/
#     thermo_capacity_rnaseq_wgcna_only.R
#   data/
#     thermo_hisat2_counts052023.txt
#     thermo_capacity_rnaseq_samples.csv
#     thermo_capacity_metadata_9Feb25.csv
#   results/
#

# Install and load required packages ---------------------------------------
cran_packages <- c(
  "dplyr"
)

bioc_packages <- c(
  "edgeR",
  "WGCNA"
)

installed <- rownames(installed.packages())

for(p in cran_packages){
  if(!(p %in% installed)){
    install.packages(p)
  }
}

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

installed <- rownames(installed.packages())

for(p in bioc_packages){
  if(!(p %in% installed)){
    BiocManager::install(p, ask = FALSE)
  }
}

library(stats)
library(dplyr)
library(edgeR) #bioconductor
library(WGCNA) #bioconductor

options(stringsAsFactors = FALSE)
enableWGCNAThreads()
options(scipen = 999)

# File paths ----------------------------------------------------------------
counts.file <- "data/thermo_hisat2_counts052023.txt"
gastroc.id.file <- "data/thermo_capacity_rnaseq_samples.csv"
lung.metadata.file <- "data/thermo_capacity_metadata_9Feb25.csv"
results.dir <- "results"

if (!dir.exists(results.dir)) {
  dir.create(results.dir, recursive = TRUE)
}

#######################################################################################################################################################
#######################################################################################################################################################
# GASTROC
#######################################################################################################################################################
#######################################################################################################################################################

# read in count and metadata ------------------------------------------------
thermo0 <- read.table(counts.file, header=TRUE, sep="\t")
colnames(thermo0)
thermo <- thermo0[, -c(2:6)] #omit chromosome, start/end, and length data
rownames(thermo) <- thermo$Geneid #replace rownames with GeneIDs
thermo$Geneid <- NULL #rownames MUST be GeneIDs. Columns must ONLY be individuals
colnames(thermo)

gastroc.id <- read.csv(gastroc.id.file)
gastroc.id <- gastroc.id[order(gastroc.id$mouse_id),]
colnames(thermo)
gastroc.id$ID
dim(thermo)
dim(gastroc.id)
colnames(thermo) <- gastroc.id$ID #WARNING. This replaces the individuals IDs, so make sure they are in order
head(thermo)
dim(thermo)

# gastroc analysis ----------------------------------------------------------
gastroc.id <- subset(gastroc.id, tissue=="gas")
gastroc.id <- gastroc.id[order(gastroc.id$mouse_id),]
gas.id <- gastroc.id$ID
gastroc <- thermo[,colnames(thermo) %in% gas.id]
colnames(gastroc)
dim(gastroc)
match(gastroc.id$ID, colnames(gastroc))
colnames(gastroc) <- gastroc.id$mouse_id
match(gastroc.id$mouse_id, colnames(gastroc))
identical(gastroc.id$mouse_id, colnames(gastroc))

# organize metadata ---------------------------------------------------------
population <- gastroc.id$population
acclimation <- gastroc.id$acclimation
gas.table <- data.frame(population, acclimation)
group <- factor(paste(gas.table$population, gas.table$acclimation, sep="_"))
cbind(gas.table, group=group)

gas.table$population = as.factor(gas.table$population)
gas.table$acclimation <- as.factor(gas.table$acclimation)
gas.table$acclimation <- relevel(gas.table$acclimation, ref="N")
gas.table$population <- relevel(gas.table$population, ref="LN")
design <- model.matrix(~population*acclimation, data=gas.table)

# filter and normalize ------------------------------------------------------
y0 <- DGEList(counts=gastroc, group=group)
keep <- filterByExpr(y0, design)
keep_gastroc <- y0[keep,,keep.lib.sizes=FALSE]
dim(keep_gastroc)

y0 <- DGEList(counts=keep_gastroc, group=group)
y <- calcNormFactors(y0)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
gastroc.norm <- cpm(y, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)
sqrt(y$common.dispersion)

# WGCNA ---------------------------------------------------------------------
head(gastroc.norm)
dim(gastroc.norm)
gastroc.Expr0 = as.data.frame(t(gastroc.norm))
rownames(gastroc.Expr0) = colnames(gastroc.norm)
gastroc.gsg = goodSamplesGenes(gastroc.Expr0, verbose = 3)
gastroc.gsg$allOK
gastroc.Expr = gastroc.Expr0

gastroc.cor <- WGCNA::cor
gastroc.Net <- blockwiseModules(
  gastroc.Expr,
  power = 7,
  maxBlockSize = dim(gastroc.Expr)[2],
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(results.dir, "ExprTOM_gastroc"),
  verbose = 3
)

save(gastroc.Net, file = file.path(results.dir, "thermo_capacity_gastroc_Network.RData"))

load(file = file.path(results.dir, "thermo_capacity_gastroc_Network.RData"))
table(gastroc.Net$colors)
gastroc.moduleLabels = gastroc.Net$colors
gastroc.moduleColors = labels2colors(gastroc.Net$colors)
length(gastroc.moduleColors)
gastroc.MEs = gastroc.Net$MEs
gastroc.geneTree = gastroc.Net$dendrograms[[1]]
table(gastroc.moduleColors)
dim(table(gastroc.moduleColors))

#######################################################################################################################################################
#######################################################################################################################################################
# LUNG
#######################################################################################################################################################
#######################################################################################################################################################

# read in lung metadata -----------------------------------------------------
lung.sample_info <- read.csv(lung.metadata.file, sep = ",", strip.white = TRUE)

lung.sample_info$group <- paste(lung.sample_info$population, lung.sample_info$acclimation, sep="_")
lung.sample_info$group <- factor(lung.sample_info$group)
lung.sample_info$family <- factor(lung.sample_info$family)
lung.sample_info$lung_sequenced <- gsub("_S[0-9]+.*", "", lung.sample_info$lung_sequenced)
lung.sample_info$lung_sequenced <- gsub("-", ".", lung.sample_info$lung_sequenced)
lung.sample_info$lung_sequenced <- gsub("_S[0-9]+.*", "", lung.sample_info$lung_sequenced)
lung.sample_info$lung_sequenced <- gsub("-", ".", lung.sample_info$lung_sequenced)

# rebuild count table for lung exactly as in lung script --------------------
thermo_counts <- read.table(counts.file, header=TRUE, sep="\t")
colnames(thermo_counts)
rownames(thermo_counts) <- thermo_counts$Geneid
thermo_counts <- thermo_counts[grep("bam$", colnames(thermo_counts))]
colnames(thermo_counts)
colnames(thermo_counts) <- gsub("^bam.", "", colnames(thermo_counts))
colnames(thermo_counts) <- gsub("_S[0-9]+.bam$", "", colnames(thermo_counts))
head(thermo_counts)
dim(thermo_counts)

# restrict count and metadata files for lung analysis -----------------------
thermo_counts_lung <- thermo_counts[, colnames(thermo_counts) %in% lung.sample_info$lung_sequenced]
lung.id <- lung.sample_info[lung.sample_info$lung_sequenced %in% colnames(thermo_counts), ]
lung.id$Sample_ID <- lung.id$lung_sequenced
rownames(lung.id) <- lung.id$Sample_ID

thermo_counts_lung <- thermo_counts_lung[, order(colnames(thermo_counts_lung))]
lung.id <- lung.id[order(rownames(lung.id)), ]

# remove known lung outliers ------------------------------------------------
lung_outliers <- c("LN.F1.113.am.lung",
                   "LN.F1.132.af.lung",
                   "LN.F1.113.b.lung",
                   "ME.F1.136.j.lung")

thermo_counts_lung <- thermo_counts_lung[, !colnames(thermo_counts_lung) %in% lung_outliers]
lung.id <- lung.id[!lung.id$Sample_ID %in% lung_outliers, ]

# formatting checks ---------------------------------------------------------
colnames(thermo_counts_lung)
head(lung.id$Sample_ID)
dim(thermo_counts_lung)
dim(lung.id)
identical(colnames(thermo_counts_lung), rownames(lung.id))

# organize metadata ---------------------------------------------------------
population <- lung.id$population
acclimation <- lung.id$acclimation
lung.table <- data.frame(population, acclimation)
group <- factor(paste(lung.table$population, lung.table$acclimation, sep="_"))
cbind(lung.table, group=group)

lung.table$population <- as.factor(lung.table$population)
lung.table$acclimation <- as.factor(lung.table$acclimation)
lung.table$acclimation <- relevel(lung.table$acclimation, ref="N")
lung.table$population <- relevel(lung.table$population, ref="LN")
design <- model.matrix(~population*acclimation, data=lung.table)
rownames(design) <- rownames(lung.id)

# filter and normalize ------------------------------------------------------
y0 <- DGEList(counts=thermo_counts_lung, group=group)
keep <- filterByExpr(y0, design)
keep_lung <- y0[keep,,keep.lib.sizes=FALSE]
dim(keep_lung)

y0 <- DGEList(counts=keep_lung, group=group)
y <- calcNormFactors(y0)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
lung.norm <- cpm(y, log=TRUE, prior.count=1, normalized.lib.sizes=TRUE)
sqrt(y$common.dispersion)

# WGCNA ---------------------------------------------------------------------
head(lung.norm)
dim(lung.norm)
lung.Expr0 = as.data.frame(t(lung.norm))
rownames(lung.Expr0) = colnames(lung.norm)
lung.gsg = goodSamplesGenes(lung.Expr0, verbose = 3)
lung.gsg$allOK
lung.Expr = lung.Expr0

lung.cor <- WGCNA::cor
lung.Net <- blockwiseModules(
  lung.Expr,
  power = 12,
  maxBlockSize = dim(lung.Expr)[2],
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(results.dir, "ExprTOM_lung"),
  verbose = 3
)

save(lung.Net, file = file.path(results.dir, "thermo_capacity_lung_Network.RData"))

load(file = file.path(results.dir, "thermo_capacity_lung_Network.RData"))
table(lung.Net$colors)
lung.moduleLabels = lung.Net$colors
lung.moduleColors = labels2colors(lung.Net$colors)
length(lung.moduleColors)
lung.MEs = lung.Net$MEs
lung.geneTree = lung.Net$dendrograms[[1]]
table(lung.moduleColors)
dim(table(lung.moduleColors))

# session info --------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(results.dir, "sessionInfo.txt"))