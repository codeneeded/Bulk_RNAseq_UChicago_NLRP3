################################################################################
## Audit of the curated-module GSEA behind Extended Data Fig. 1b.
##
## Checks three things the original ENRICHMENT ANALYSIS block in
## DGE_and_Pathway_Analysis.R does not report:
##
##   1. How many genes of each human module survive conversion to mouse
##      orthologs (mapIds ALIAS -> ENTREZID -> Orthology.eg.db).  Unmapped
##      genes disappear silently.
##   2. Whether a module was dropped by the clusterProfiler::GSEA size filter.
##      GSEA() defaults to minGSSize = 10; Non_Phagocytosis starts with 8 human
##      genes, so it never reaches the plot.
##   3. Whether the enrichment calls depend on the ranking metric.  The original
##      ranks genes by raw log2FoldChange, which at n = 3 puts low-count genes
##      with extreme, poorly estimated fold changes at the head and tail of the
##      list.  Here the same modules are re-tested ranking by the Wald statistic.
##
## Also runs the modules against the sample-matched M2-vs-M1 contrasts from
## M2vsM1_WithinCondition_and_Fig1b_Heatmap.R.
##
## In:  Raw_Data/raw_combined_featureCounts.txt
## Out: Reviewer_Comment/Module_Analysis/Ortholog_Mapping/module_ortholog_map.csv
##      Reviewer_Comment/Module_Analysis/Ortholog_Mapping/module_sizes.csv
##      Reviewer_Comment/Module_Analysis/Ranking_Robustness/module_gsea_by_ranking.csv
################################################################################

suppressPackageStartupMessages({
  library(DESeq2); library(clusterProfiler)
  library(org.Hs.eg.db); library(org.Mm.eg.db); library(Orthology.eg.db)
})
set.seed(42)

repo <- Sys.getenv("REPO_ROOT", unset = ".")
od1 <- file.path(repo, "Reviewer_Comment/Module_Analysis/Ortholog_Mapping")
od2 <- file.path(repo, "Reviewer_Comment/Module_Analysis/Ranking_Robustness")
dir.create(od1, recursive = TRUE, showWarnings = FALSE)
dir.create(od2, recursive = TRUE, showWarnings = FALSE)

## ---- modules exactly as defined in DGE_and_Pathway_Analysis.R ---------------

genes_M1_pol <- c("AK3","APOL1","APOL2","APOL3","APOL6","ATF3","BCL2A1","BIRC3","CCL5","CCL15",
  "CCL19","CCL20","CCR7","CHI3L2","CSPG2","CXCL9","CXCL10","CXCL11","ECGF1","EDN1","FAS","GADD45G",
  "HESX1","HSD11B1","HSXIAPAF1","IGFBP4","IL6","IL12B","IL15","IL2RA","IL7R","IL15RA","INDO","INHBA",
  "IRF1","IRF7","OAS2","OASL","PBEF1","PDGFA","PFKFB3","PFKP","PLA1A","PSMA2","PSMB9","PSME2","PTX3",
  "SLC21A15","SLC2A6","SLC31A2","SLC7A5","SPHK1","TNF","TRAIL")
genes_M2_pol <- c("ADK","ALOX15","CA2","CCL13","CCL18","CCL23","CD36","CERK","CHN2","CLECSF13","CTSC",
  "CXCR4","DCL-1","DCSIGN","DECTIN1","EGR2","FGL2","FN1","GAS7","GPR86","HEXB","HNMT","HRH1","HS3ST1",
  "HS3ST2","IGF1","LIPA","LTA4H","MAF","MRC1","MS4A4A","MS4A6A","MSR1","P2RY14","P2RY5","SEPP1",
  "SLC21A9","SLC38A6","SLC4A7","TGFB1","TGFBR2","TLR5","TPST2")
genes_inflammasome <- c("AIM2","TRAF6","TXNIP","XIAP","TNFAIP3","BCL2","CASP1","CASP12","CASP4","CARD6",
  "CTSB","CTSV","IL18","IL1A","IL1B","IL18BP","IL33","MYD88","IRF1","IRF3","NFkB1","NFKBID","NFKBIB",
  "NLRP3","PSTPIP1","PTGS2","PYCARD","TIRAP","TNFSF14","TNFSF4")
genes_pyroptosis <- c("AIM2","CASP1","CASP3","CASP4","CASP5","CASP6","CASP8","CASP9","ELANE","GPX4",
  "GSDMA","GSDMB","GSDMC","GSDMD","GSDME","IL1A","IL1B","IL6","IL18","NOD1","NOD2","NLRC4","NLRP1",
  "NLRP2","NLRP3","NLRP6","NLRP7","PJVK","PLCG1","PRKACA","PYCARD","SCAF11","TNF")
genes_non_phagocytosis <- c("ATAD1","B2M","FCGR1","IGHM","IL4","MFF","STAP1","SYNE1")
genes_aerobic_glycolysis <- c("ALDOA","BSG","ENO1","ENO2","ENO3","GPI","HK1","LDHA","LDHB","PFKL",
  "PFKM","PFKP","PGAM1","PGAM2","PGK1","PKLR","PKM","SLC16A1","SLC16A3","SLC2A1","TPI1")

gene_sets_human <- list(M1_polarization = genes_M1_pol, M2_polarization = genes_M2_pol,
                        Inflammasome = genes_inflammasome, Pyroptosis = genes_pyroptosis,
                        Non_Phagocytosis = genes_non_phagocytosis,
                        Aerobic_Glycolysis = genes_aerobic_glycolysis)

## ---- ortholog mapping, recording every failure ------------------------------

map_one <- function(symbols) {
  human <- suppressMessages(mapIds(org.Hs.eg.db, symbols, "ENTREZID", "ALIAS"))
  out <- data.frame(Human.Symbol = symbols, Human.Entrez = as.character(human),
                    stringsAsFactors = FALSE)
  ok <- !is.na(out$Human.Entrez)
  out$Mouse.Entrez <- NA_character_; out$Mouse.Symbol <- NA_character_
  if (any(ok)) {
    mapped <- suppressMessages(select(Orthology.eg.db, out$Human.Entrez[ok],
                                      "Mus.musculus", "Homo.sapiens"))
    mapped <- mapped[match(out$Human.Entrez[ok], as.character(mapped[, 1])), ]
    out$Mouse.Entrez[ok] <- as.character(mapped[, 2])
    has <- ok & !is.na(out$Mouse.Entrez)
    if (any(has))
      out$Mouse.Symbol[has] <- suppressMessages(
        mapIds(org.Mm.eg.db, out$Mouse.Entrez[has], "SYMBOL", "ENTREZID"))
  }
  out$failed_at <- ifelse(is.na(out$Human.Entrez), "human symbol not recognised",
                   ifelse(is.na(out$Mouse.Entrez), "no mouse ortholog",
                   ifelse(is.na(out$Mouse.Symbol), "no mouse symbol", "")))
  out
}

map_all <- do.call(rbind, lapply(names(gene_sets_human), function(m) {
  d <- map_one(gene_sets_human[[m]]); d$Module <- m; d
}))
write.csv(map_all, file.path(od1, "module_ortholog_map.csv"), row.names = FALSE)

gene_sets_mouse <- lapply(gene_sets_human, function(g) {
  m <- map_one(g); unique(na.omit(m$Mouse.Symbol))
})
term2gene <- do.call(rbind, lapply(names(gene_sets_mouse), function(s)
  data.frame(term = s, gene = gene_sets_mouse[[s]], stringsAsFactors = FALSE)))

## ---- counts and DESeq2 ------------------------------------------------------

fc <- read.table(file.path(repo, "Raw_Data/raw_combined_featureCounts.txt"),
                 header = TRUE, sep = "\t", comment.char = "#")
count_data <- fc[, 7:ncol(fc)]; rownames(count_data) <- fc$Geneid
md <- data.frame(
  sample = c("M1_U_1","M1_U_2","M1_ADT_3","M1_ADT_1","M1_ADT_2","M1_U_3","M1_NLRP3_1","M1_NLRP3_2",
             "M1_Combo_3","M1_Combo_1","M1_Combo_2","M1_NLRP3_3","M2_U_1","M2_U_2","M2_U_3",
             "M2_ADT_1","M2_ADT_2","M2_ADT_3","M2_NLRP3_1","M2_NLRP3_2","M2_NLRP3_3",
             "M2_Combo_1","M2_Combo_2","M2_Combo_3"),
  file.name = colnames(count_data))
md$cell_type <- substr(md$sample, 1, 2)
md$condition <- sub("^M[12]_", "", sub("_[0-9]$", "", md$sample))
rownames(md) <- md$file.name

dds_M1 <- DESeq(DESeqDataSetFromMatrix(count_data[, md$cell_type == "M1"],
                                       md[md$cell_type == "M1", ], ~ condition))
dds_M2 <- DESeq(DESeqDataSetFromMatrix(count_data[, md$cell_type == "M2"],
                                       md[md$cell_type == "M2", ], ~ condition))

res_list <- list()
for (cmp in list(c("ADT","U"), c("NLRP3","U"), c("Combo","U"), c("Combo","ADT"),
                 c("Combo","NLRP3"), c("NLRP3","ADT"))) {
  res_list[[paste0("M1_", cmp[1], "_vs_", cmp[2])]] <-
    as.data.frame(results(dds_M1, contrast = c("condition", cmp[1], cmp[2])))
  res_list[[paste0("M2_", cmp[1], "_vs_", cmp[2])]] <-
    as.data.frame(results(dds_M2, contrast = c("condition", cmp[1], cmp[2])))
}
for (cond in c("U","ADT","NLRP3","Combo")) {
  k <- md$condition == cond
  d <- DESeq(DESeqDataSetFromMatrix(count_data[, k], md[k, ], ~ cell_type))
  res_list[[paste0("M2vsM1_", cond)]] <-
    as.data.frame(results(d, contrast = c("cell_type", "M2", "M1")))
}

## ---- GSEA under two ranking metrics, minGSSize lowered to 5 -----------------

run_gsea <- function(rank_vec, label, contrast) {
  rank_vec <- sort(rank_vec[is.finite(rank_vec)], decreasing = TRUE)
  g <- GSEA(geneList = rank_vec, TERM2GENE = term2gene, pvalueCutoff = 1,
            minGSSize = 5, maxGSSize = 500, seed = TRUE, verbose = FALSE)
  if (nrow(g@result) == 0) return(NULL)
  data.frame(Contrast = contrast, Ranking = label, Module = g@result$ID,
             setSize = g@result$setSize, NES = g@result$NES,
             pvalue = g@result$pvalue, p.adjust = g@result$p.adjust,
             row.names = NULL)
}

out <- list()
for (nm in names(res_list)) {
  r <- res_list[[nm]]
  out[[paste0(nm, "_lfc")]]  <- run_gsea(setNames(r$log2FoldChange, rownames(r)), "log2FoldChange", nm)
  out[[paste0(nm, "_stat")]] <- run_gsea(setNames(r$stat, rownames(r)), "Wald statistic", nm)
}
res_tab <- do.call(rbind, out)
write.csv(res_tab, file.path(od2, "module_gsea_by_ranking.csv"), row.names = FALSE)

sizes <- data.frame(
  Module = names(gene_sets_human),
  human_genes = sapply(gene_sets_human, length),
  mouse_orthologs = sapply(gene_sets_mouse, length),
  detected_M1_ADT_vs_U = sapply(names(gene_sets_mouse), function(m)
    length(intersect(gene_sets_mouse[[m]], rownames(res_list[["M1_ADT_vs_U"]])))),
  dropped_at_minGSSize_10 = sapply(gene_sets_mouse, function(g) length(g) < 10),
  row.names = NULL)
write.csv(sizes, file.path(od1, "module_sizes.csv"), row.names = FALSE)

cat("\n== module sizes ==\n"); print(sizes)
cat("\n== Ext. Data Fig. 1b contrasts, both rankings ==\n")
print(res_tab[res_tab$Contrast %in% c("M1_ADT_vs_U","M2_ADT_vs_U"),
              c("Contrast","Ranking","Module","setSize","NES","p.adjust")], row.names = FALSE)
cat("\n== sample-matched M2 vs M1, Wald ranking ==\n")
print(res_tab[grepl("^M2vsM1", res_tab$Contrast) & res_tab$Ranking == "Wald statistic",
              c("Contrast","Module","setSize","NES","p.adjust")], row.names = FALSE)

## ---- sensitivity: restore orthologs the automatic mapper missed -------------
## Orthology.eg.db returns no mouse ortholog for several human genes that do
## have an unambiguous same-symbol mouse counterpart, including CASP1 and
## PYCARD -- the two most central inflammasome genes.  NFkB1 additionally fails
## because mapIds is case-sensitive.  Only unambiguous 1:1 same-symbol rescues
## are added here; genuinely absent (APOL1-3, CASP5, NLRP7) or one-to-many
## (GSDMA/C, BCL2A1, MS4A6A, CD209) cases are left out and listed in
## module_ortholog_map.csv for manual curation.

rescue <- list(CASP1 = "Casp1", PYCARD = "Pycard", NFkB1 = "Nfkb1", CXCL10 = "Cxcl10",
               ALOX15 = "Alox15", BIRC3 = "Birc3", OAS2 = "Oas2", IGHM = "Ighm")
detected <- rownames(res_list[["M1_ADT_vs_U"]])
rescue <- rescue[unlist(rescue) %in% detected]
cat("\nrescued orthologs present in the data:", paste(unlist(rescue), collapse = ", "), "\n")

gene_sets_fixed <- lapply(names(gene_sets_human), function(m) {
  add <- unlist(rescue[intersect(names(rescue), gene_sets_human[[m]])])
  unique(c(gene_sets_mouse[[m]], add))
})
names(gene_sets_fixed) <- names(gene_sets_human)
term2gene_fixed <- do.call(rbind, lapply(names(gene_sets_fixed), function(s)
  data.frame(term = s, gene = gene_sets_fixed[[s]], stringsAsFactors = FALSE)))

t2g_backup <- term2gene; term2gene <- term2gene_fixed
sens <- do.call(rbind, lapply(c("M1_ADT_vs_U", "M2_ADT_vs_U", "M2vsM1_U"), function(nm) {
  r <- res_list[[nm]]
  run_gsea(setNames(r$stat, rownames(r)), "Wald statistic, repaired modules", nm)
}))
term2gene <- t2g_backup

write.csv(rbind(res_tab, sens), file.path(od2, "module_gsea_by_ranking.csv"), row.names = FALSE)
cat("\n== repaired modules, Wald ranking ==\n")
print(sens[, c("Contrast","Module","setSize","NES","p.adjust")], row.names = FALSE)
