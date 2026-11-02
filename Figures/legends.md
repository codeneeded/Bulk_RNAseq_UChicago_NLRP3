# Figure legends — bulk RNA-seq panels

Panels contributed by this repository: Fig. 1b and Extended Data Fig. 1b.
Full text, with the methods section, is in Methods_and_Figure_Legends.docx (not tracked here).

## Fig. 1b
Panel: `Figures/panels/Fig1b_DGE_heatmap.png` — source data `Figures/source_data/Fig1b_log2FC_padj.csv`
Built by `M2vsM1_WithinCondition_and_Fig1b_Heatmap.R`.

DESeq2 log2 fold changes for a 50-gene panel across three comparisons: M2 ARB vs U, M1 ARB vs U,
and M2 vs M1 in untreated cells. Grey = no estimable fold change, white = no change,
scale clipped at ±4. Asterisks: Benjamini–Hochberg adjusted p < 0.05 (*), < 0.01 (**), < 0.001 (***).
n = 3 biological replicates per group.

Supersedes `Reviewer_Comment/Raw/DGE_heatmap.png`, whose third column was the M2-vs-M1 difference
under enzalutamide rather than in untreated cells, on a quarter of the genes, and which starred
genes that were never tested. Both files are kept; the panel above is the one that ships.

## Extended Data Fig. 1b
Panels: `Reviewer_Comment/Module_Analysis/GSEA_Plots/GSEA_M1_ADT_vs_U.png` and `GSEA_M2_ADT_vs_U.png`.

Curated-module GSEA, genes ranked by log2 fold change, five modules (non-phagocytosis fell below the
minimum set size). The dotplot does not encode direction: in M1 all three significant modules are
negatively enriched — M1 polarization NES −2.17, inflammasome −1.77, pyroptosis −1.75. Nothing reaches
adjusted p < 0.05 in M2. Robustness under an alternative ranking and a repaired ortholog mapping:
`Reviewer_Comment/Module_Analysis/Ranking_Robustness/module_gsea_by_ranking.csv`.

## QC
`Figures/panels/QC_polarization_and_treatment_response.png` — source data
`Figures/source_data/QC_PCA_and_markers.csv`. Not a manuscript panel; documents that the
enzalutamide-treated M1 arms do not carry an M1 activation programme.
