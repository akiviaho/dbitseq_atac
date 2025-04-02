library(ArchR)
library(parallel)
set.seed(1)


# PARAMETERS
n_threads <- 16
data_dir <- "/lustre/scratch/kiviaho/dbitseq_atac"
inputFiles <- c("/lustre/scratch/kiviaho/dbitseq_atac/sample_A_Sorted.bed.gz")
names(inputFiles) <- c("sample_A")
proj_name <- "test_project"


# Configure
addArchRThreads(threads = 15) 
addArchRGenome('hg38')


# ################################################################################################
# Arrow files and project 

setwd(sprintf("%s/archr_outs", data_dir))

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS =0, # No QC
  minFrags =0, # No QC
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
#  force=TRUE
)


# Create project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = proj_name,
  copyArrows = TRUE
)


# ################################################################################################
# Preprocessing

# SVD, Clustering, UMAP
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                       name = "IterativeLSI", force=TRUE)


# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features <- proj@reducedDims[["IterativeLSI"]]$LSIFeatures
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)


# Peaks using NFR fragments
proj <- addClusters(input = proj, reducedDims = "IterativeLSI") # Call clusters

# Call peaks samplewise instead of clusterwise
proj <- addGroupCoverages(proj,groupBy='Sample')
proj <- addReproduciblePeakSet(proj,groupBy='Sample')
# Counts
proj <- addPeakMatrix(proj)

# Save 
proj <- saveArchRProject(ArchRProj = proj)


# ################################################################################################

# Export
dir.create(sprintf("%s/export", proj_name))
write.csv(getReducedDims(proj), sprintf('%s/export/svd.csv', proj_name), quote=FALSE)
write.csv(getCellColData(proj), sprintf('%s/export/cell_metadata.csv', proj_name), quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
scores <- as.matrix(scores)
rownames(scores) <- rowData(gene.scores)$name
write.csv(scores, sprintf('%s/export/gene_scores.csv', proj_name), quote=FALSE)

# Peak counts
peaks <- getPeakSet(proj)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')

# Reorder peaks 
# Chromosome order
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order)
    reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
reordered_features <- Reduce("c", reordered_features)

# Export counts
dir.create(sprintf("%s/export/peak_counts", proj_name))
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, sprintf('%s/export/peak_counts/counts.mtx', proj_name))
write.csv(colnames(peak.counts), sprintf('%s/export/peak_counts/cells.csv', proj_name), quote=FALSE)
names(reordered_features) <- sprintf("Peak%d", 1:length(reordered_features))
write.csv(as.data.frame(reordered_features), sprintf('%s/export/peak_counts/peaks.csv', proj_name), quote=FALSE)


# Plot fragment size and TSS enrichment plots
p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
