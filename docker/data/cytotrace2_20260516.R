library(Seurat)
library(CytoTRACE2)
library(dplyr)
library(ggplot2)
library(getopt)

command = matrix(c("rds","r",1,"character",
                   "outpath","o",1,"character",
                   "species","s",2,"character",
                   "ncore","n",2,"numeric",
                   "prefix","h",2,"character",
                   "slot_type","l",2,"character")
                 ,byrow = TRUE, ncol = 4)

args = getopt(command)
if (is.null(args$rds)) stop("missing -r/--rds")
if (is.null(args$outpath)) stop("missing -o/--outpath")

output <- args$outpath
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
prefix <- if (!is.null(args$prefix)) args$prefix else "CytoTRACE2"
species <- if (!is.null(args$species)) tolower(args$species) else "human"
core <- if (!is.null(args$ncore)) as.numeric(args$ncore) else 8
slot_type <- if (!is.null(args$slot_type)) args$slot_type else "counts"

if (!species %in% c("human", "mouse")) {
  stop("species must be 'human' or 'mouse'")
}

message("[Start] CytoTRACE2")
message("input=", args$rds)
message("output=", output)
message("species=", species)
message("ncore=", core)

seuset <- readRDS(args$rds)
DefaultAssay(seuset) <- "RNA"
if ("RNA" %in% names(seuset@assays) && inherits(seuset[["RNA"]], "Assay5") && exists("JoinLayers")) {
  seuset <- tryCatch(
    JoinLayers(seuset, assay = "RNA"),
    error = function(e) {
      warning("JoinLayers skipped: ", conditionMessage(e))
      seuset
    }
  )
}

expr <- tryCatch(
  GetAssayData(seuset, assay = "RNA", layer = slot_type),
  error = function(e) {
    warning("GetAssayData layer='", slot_type, "' failed, fallback to slot='", slot_type, "': ", conditionMessage(e))
    GetAssayData(seuset, assay = "RNA", slot = slot_type)
  }
)
expr <- as.matrix(expr)
message("matrix genes=", nrow(expr), " cells=", ncol(expr))

results <- CytoTRACE2::cytotrace2(
  expr,
  species = species,
  is_seurat = FALSE,
  slot_type = slot_type,
  ncores = core
)
results <- as.data.frame(results)
if (is.null(rownames(results)) || any(rownames(results) == "")) {
  stop("CytoTRACE2 result has no cell rownames")
}

write.table(
  cbind(cell = rownames(results), results),
  file = file.path(output, paste0(prefix, ".table.txt")),
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
saveRDS(results, file.path(output, paste0(prefix, ".result.rds")))

common_cells <- intersect(colnames(seuset), rownames(results))
if (length(common_cells) == 0) {
  stop("No CytoTRACE2 result cells matched the Seurat object")
}
seuset <- AddMetaData(seuset, metadata = results[common_cells, , drop = FALSE])
saveRDS(seuset, file.path(output, paste0(prefix, ".seurat.rds")))

plots <- CytoTRACE2::plotData(seuset, is_seurat = TRUE, slot_type = slot_type)
for (nm in names(plots)) {
  ggsave(file.path(output, paste0(prefix, ".", nm, ".pdf")), plots[[nm]], width = 8.5, height = 8)
  ggsave(file.path(output, paste0(prefix, ".", nm, ".png")), plots[[nm]], width = 8.5, height = 8)
}

if ("CytoTRACE2_Score" %in% colnames(seuset@meta.data)) {
  p <- FeaturePlot(seuset, "CytoTRACE2_Score", raster = FALSE, reduction = "umap")
  ggsave(file.path(output, paste0(prefix, ".score.FeaturePlot.pdf")), p, width = 8.5, height = 8)
  ggsave(file.path(output, paste0(prefix, ".score.FeaturePlot.png")), p, width = 8.5, height = 8)
}
if ("CytoTRACE2_Relative" %in% colnames(seuset@meta.data)) {
  p <- FeaturePlot(seuset, "CytoTRACE2_Relative", raster = FALSE, reduction = "umap")
  ggsave(file.path(output, paste0(prefix, ".relative.FeaturePlot.pdf")), p, width = 8.5, height = 8)
  ggsave(file.path(output, paste0(prefix, ".relative.FeaturePlot.png")), p, width = 8.5, height = 8)
}
if ("CytoTRACE2_Potency" %in% colnames(seuset@meta.data)) {
  p <- DimPlot(seuset, reduction = "umap", group.by = "CytoTRACE2_Potency")
  ggsave(file.path(output, paste0(prefix, ".potency.DimPlot.pdf")), p, width = 8.5, height = 8)
  ggsave(file.path(output, paste0(prefix, ".potency.DimPlot.png")), p, width = 8.5, height = 8)
}

cluster_col <- if ("seurat_clusters" %in% colnames(seuset@meta.data)) "seurat_clusters" else NULL
if (!is.null(cluster_col) && "CytoTRACE2_Score" %in% colnames(seuset@meta.data)) {
  cyto_ident <- data.frame(
    Cluster = as.character(seuset@meta.data[[cluster_col]]),
    CytoTRACE2_Score = seuset@meta.data$CytoTRACE2_Score,
    CytoTRACE2_Relative = seuset@meta.data$CytoTRACE2_Relative,
    CytoTRACE2_Potency = seuset@meta.data$CytoTRACE2_Potency,
    row.names = rownames(seuset@meta.data)
  )
  cyto_order <- cyto_ident %>%
    group_by(Cluster) %>%
    summarise(score_median = median(CytoTRACE2_Score, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score_median))
  cyto_ident$Cluster <- factor(cyto_ident$Cluster, levels = cyto_order$Cluster)
  write.table(
    cbind(cell = rownames(cyto_ident), cyto_ident),
    file = file.path(output, paste0(prefix, ".cluster_table.txt")),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
  p <- ggplot(cyto_ident, aes(x = Cluster, y = CytoTRACE2_Score, fill = Cluster, color = Cluster)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.35, size = 0.4) +
    labs(x = "Cluster", y = "CytoTRACE2 score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")
  ggsave(file.path(output, paste0(prefix, ".score.boxplot_cluster.pdf")), p, width = max(8, length(unique(cyto_ident$Cluster)) * 0.45), height = 8)
  ggsave(file.path(output, paste0(prefix, ".score.boxplot_cluster.png")), p, width = max(8, length(unique(cyto_ident$Cluster)) * 0.45), height = 8)
}

message("[End] CytoTRACE2")
