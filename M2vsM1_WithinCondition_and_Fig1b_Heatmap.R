################################################################################
## M2-vs-M1 contrasts computed within each stimulation condition, and the
## Fig. 1b log2FC heatmap rebuilt from them.
##
## Supersedes two steps of DGE_and_Pathway_Analysis.R:
##
##   1. The "Differential Expression M2 Vs M1" loop subset the DESeqDataSet as
##      dds[metadata$condition == condition, ].  The first index of a
##      DESeqDataSet selects GENES, not samples, so a 24-element logical was
##      recycled across ~18,000 genes.  The four DGE/M2vsM1/*.csv files are
##      consequently four disjoint quarters of the gene list, and each carries
##      the cell_type main effect of the interaction model -- the M2-vs-M1
##      difference at the REFERENCE level of condition (ADT), not at the
##      condition named in the file.
##      Here each condition gets its own DESeqDataSet built from that
##      condition's six samples with design ~ cell_type.
##
##   2. The heatmap block ran combined_data[is.na(combined_data)] <- 0, which
##      turned a missing padj into 0 and therefore into "***" on a cell whose
##      log2FC had also been set to 0.  Here genes absent from a contrast stay
##      NA: drawn grey, never starred.
##
## In:  Raw_Data/raw_combined_featureCounts.txt
## Out: DGE/M2vsM1_WithinCondition/M2vsM1.{U,ADT,NLRP3,Combo}.DGE.csv
##      Figures/panels/Fig1b_DGE_heatmap.png
##      Figures/source_data/Fig1b_log2FC_padj.csv
################################################################################

library(DESeq2)
library(pheatmap)

repo <- Sys.getenv("REPO_ROOT", unset = ".")

dir.create(file.path(repo, "DGE/M2vsM1_WithinCondition"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo, "Figures/panels"),             recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo, "Figures/source_data"),        recursive = TRUE, showWarnings = FALSE)

## ---- counts and sample table ------------------------------------------------

feature_counts <- read.table(file.path(repo, "Raw_Data/raw_combined_featureCounts.txt"),
                             header = TRUE, sep = "\t", comment.char = "#")
count_data <- feature_counts[, 7:ncol(feature_counts)]
rownames(count_data) <- feature_counts$Geneid

metadata <- data.frame(
  sample = c("M1_U_1", "M1_U_2", "M1_ADT_3", "M1_ADT_1", "M1_ADT_2", "M1_U_3",
             "M1_NLRP3_1", "M1_NLRP3_2", "M1_Combo_3", "M1_Combo_1", "M1_Combo_2", "M1_NLRP3_3",
             "M2_U_1", "M2_U_2", "M2_U_3", "M2_ADT_1", "M2_ADT_2", "M2_ADT_3",
             "M2_NLRP3_1", "M2_NLRP3_2", "M2_NLRP3_3", "M2_Combo_1", "M2_Combo_2", "M2_Combo_3"),
  condition = c("M1_U", "M1_U", "M1_ADT", "M1_ADT", "M1_ADT", "M1_U",
                "M1_NLRP3", "M1_NLRP3", "M1_Combo", "M1_Combo", "M1_Combo", "M1_NLRP3",
                "M2_U", "M2_U", "M2_U", "M2_ADT", "M2_ADT", "M2_ADT",
                "M2_NLRP3", "M2_NLRP3", "M2_NLRP3", "M2_Combo", "M2_Combo", "M2_Combo"),
  file.name = colnames(count_data)
)
metadata$cell_type <- ifelse(grepl("M1", metadata$condition), "M1", "M2")
metadata$condition <- gsub("M1_|M2_", "", metadata$condition)
rownames(metadata) <- metadata$file.name
metadata <- metadata[colnames(count_data), ]

stopifnot(identical(rownames(metadata), colnames(count_data)))

## ---- M2 vs M1 within each condition ----------------------------------------

m2vsm1 <- list()
for (cond in c("U", "ADT", "NLRP3", "Combo")) {
  keep <- metadata$condition == cond
  message(sprintf("M2 vs M1 in %s: %d samples (%s)", cond, sum(keep),
                  paste(metadata$sample[keep], collapse = ", ")))
  dds_c <- DESeqDataSetFromMatrix(countData = count_data[, keep],
                                  colData   = metadata[keep, ],
                                  design    = ~ cell_type)
  dds_c <- DESeq(dds_c)
  res <- as.data.frame(results(dds_c, contrast = c("cell_type", "M2", "M1")))
  res <- na.omit(res)
  m2vsm1[[cond]] <- res
  write.csv(res, file.path(repo, "DGE/M2vsM1_WithinCondition",
                           sprintf("M2vsM1.%s.DGE.csv", cond)))
}

## ---- within-polarisation contrasts (unchanged in substance) -----------------

dds_M1 <- DESeq(DESeqDataSetFromMatrix(count_data[, metadata$cell_type == "M1"],
                                       metadata[metadata$cell_type == "M1", ], ~ condition))
dds_M2 <- DESeq(DESeqDataSetFromMatrix(count_data[, metadata$cell_type == "M2"],
                                       metadata[metadata$cell_type == "M2", ], ~ condition))

res_M1_ADT_vs_U <- na.omit(as.data.frame(results(dds_M1, contrast = c("condition", "ADT", "U"))))
res_M2_ADT_vs_U <- na.omit(as.data.frame(results(dds_M2, contrast = c("condition", "ADT", "U"))))

## agreement with the committed tables, which were subset correctly
for (chk in list(list("DGE/M1/M1.ADT.vs.U.DGE.csv", res_M1_ADT_vs_U),
                 list("DGE/M2/M2.ADT.vs.U.DGE.csv", res_M2_ADT_vs_U))) {
  old <- read.csv(file.path(repo, chk[[1]]), row.names = 1)
  g <- intersect(rownames(old), rownames(chk[[2]]))
  message(sprintf("%s: %d genes shared with recomputation, max |dlog2FC| = %.2e",
                  chk[[1]], length(g),
                  max(abs(old[g, "log2FoldChange"] - chk[[2]][g, "log2FoldChange"]))))
}

## ---- Fig. 1b panel ----------------------------------------------------------

panel_genes <- c("Spp1", "Ikbkb", "Traf6", "Ifnl2", "Mavs", "Eftud2", "Ddx58",
                 "Rela", "Irf3", "Il1a", "Ifih1", "Prkra", "Gsdmd", "Id2",
                 "Tlr2", "Ifnb1", "Ap3b1", "Nlrp3", "Nod2", "Il10rb",
                 "Tlr9", "Il17ra", "Sod1", "Hrh4", "Ifnar2", "Sdf2l1",
                 "Vpreb1", "Pnp", "Igll1", "Rnase7", "Siglec15", "Psmb5",
                 "Cebpe", "Il25", "Cmtm5", "Psme1", "Psme2", "Nfatc4",
                 "Mx1", "Lif", "Ctsg", "Gzmb", "Tcn2", "Tnfrsf11a",
                 "Nfkbia", "Tmem173", "Stat1", "Il22ra1", "Bcl2l1",
                 "Jak1", "Jak2", "Jak3", "Tyk2", "Cd68", "Trem2",
                 "C1qa", "C1qb", "C1qc", "Apoe", "Cd163")

contrasts_used <- list("M2: ARB vs U" = res_M2_ADT_vs_U,
                       "M1: ARB vs U" = res_M1_ADT_vs_U,
                       "U: M2 vs M1"  = m2vsm1[["U"]])

pull <- function(res, col) {
  v <- res[panel_genes, col]
  names(v) <- panel_genes
  v
}
logFC <- sapply(contrasts_used, pull, col = "log2FoldChange")
padj  <- sapply(contrasts_used, pull, col = "padj")
rownames(logFC) <- rownames(padj) <- panel_genes

## genes with no testable value in any contrast carry no information
keep_gene <- rowSums(!is.na(logFC)) > 0
logFC <- logFC[keep_gene, , drop = FALSE]
padj  <- padj[keep_gene, , drop = FALSE]

message(sprintf("panel genes: %d requested, %d testable in >=1 contrast", length(panel_genes), nrow(logFC)))
for (j in colnames(logFC))
  message(sprintf("  %-14s tested %2d/%2d, FDR<0.05 in %2d",
                  j, sum(!is.na(logFC[, j])), nrow(logFC), sum(padj[, j] < 0.05, na.rm = TRUE)))

stars <- ifelse(is.na(padj), "",
         ifelse(padj <= 0.001, "***",
         ifelse(padj <= 0.01,  "**",
         ifelse(padj <= 0.05,  "*", ""))))

## row order: cluster on the tested values, NA treated as no change for distance only
row_order <- hclust(dist(ifelse(is.na(logFC), 0, logFC)))

## Colour: white is anchored at log2FC = 0 by splitting the ramp at zero.  The
## original code built a single ramp over the data range, which placed white at
## about -0.8 and rendered every near-zero cell as saturated green.  The scale
## is clipped at +/- CLIP because the range runs to -14.2 while 91% of cells sit
## within +/- 4; exact values are in Figures/source_data/Fig1b_log2FC_padj.csv.
CLIP <- 4
pal <- c(colorRampPalette(c("darkblue", "purple", "magenta", "lightpink", "white"))(50),
         colorRampPalette(c("white", "green"))(50))
brk <- seq(-CLIP, CLIP, length.out = 101)
logFC_disp <- pmin(pmax(logFC, -CLIP), CLIP)

png(file.path(repo, "Figures/panels/Fig1b_DGE_heatmap.png"),
    width = 780, height = 2300, res = 200)
pheatmap(logFC_disp,
         cluster_rows    = row_order,
         cluster_cols    = FALSE,
         display_numbers = stars,
         fontsize_number = 12,
         number_color    = "red",
         na_col          = "grey85",
         main            = sprintf("log2FC, clipped +/-%d\ngrey = not tested", CLIP),
         color           = pal,
         breaks          = brk,
         scale           = "none")
dev.off()

src <- data.frame(Gene = rownames(logFC),
                  log2FC_M2_ARB_vs_U = logFC[, 1], padj_M2_ARB_vs_U = padj[, 1],
                  log2FC_M1_ARB_vs_U = logFC[, 2], padj_M1_ARB_vs_U = padj[, 2],
                  log2FC_U_M2_vs_M1  = logFC[, 3], padj_U_M2_vs_M1  = padj[, 3],
                  row.names = NULL)
write.csv(src, file.path(repo, "Figures/source_data/Fig1b_log2FC_padj.csv"), row.names = FALSE)

sessionInfo()
