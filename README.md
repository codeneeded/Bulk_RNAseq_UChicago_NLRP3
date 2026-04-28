# Bulk RNA-seq of Polarized Macrophages under NLRP3 Inflammasome Stimulation

Bulk RNA-seq analysis of **mouse bone marrow–derived macrophages** *(confirm cell source)* polarized to **M1** or **M2** states and challenged with **NLRP3 inflammasome–activating stimuli**. The design tests how prior polarization (M1 vs M2) shapes the transcriptional response to NLRP3 activation, and whether a combined stimulus (Combo) drives a response distinct from either single agonist.

This is a collaboration with the **University of Chicago**.

> *Note on naming.* The condition `ADT` here refers to the lab-internal abbreviation for one of the inflammasome stimuli used in this study and is **unrelated** to "ADT" in the CITE-seq sense (Antibody-Derived Tags). *(Spell out the actual reagent: e.g. ATP, alum, etc.)*

---

## Experimental design

| Factor | Levels | Notes |
|---|---|---|
| **Polarization** | M1, M2 | *(e.g. LPS+IFNγ for M1, IL-4 for M2 — fill in)* |
| **Stimulus** | U (unstimulated), ADT, NLRP3, Combo | U = vehicle; Combo = ADT + NLRP3 *(confirm)* |
| **Replicates** | n = 3 biological replicates per condition × polarization | 24 libraries total |
| **Species** | *Mus musculus* | mm10 / GRCm38 *(confirm assembly)* |
| **Library prep** | *(e.g. polyA-selected stranded mRNA — fill in)* | |
| **Sequencing** | *(platform, read length, single/paired — fill in)* | |

Sample naming follows `{Polarization}_{Condition}_{Replicate}` — e.g. `M1_NLRP3_2`, `M2_Combo_3`.

---

## Repository structure

```
Bulk_RNAseq_UChicago_NLRP3/
├── DGE_and_Pathway_Analysis.R           # End-to-end analysis script (DESeq2 → clusterProfiler → pathview)
│
├── Raw_Data/
│   └── raw_combined_featureCounts.txt   # featureCounts matrix (subread), all 24 samples
│
├── DGE/                                 # DESeq2 differential expression results
│   ├── M1/                              # Within-M1 contrasts
│   │   ├── M1.ADT.vs.U.DGE.csv
│   │   ├── M1.NLRP3.vs.U.DGE.csv
│   │   ├── M1.Combo.vs.U.DGE.csv
│   │   ├── M1.NLRP3.vs.ADT.DGE.csv
│   │   ├── M1.Combo.vs.ADT.DGE.csv
│   │   └── M1.Combo.vs.NLRP3.DGE.csv
│   ├── M2/                              # Within-M2 contrasts (same six)
│   │   └── M2.{...}.DGE.csv
│   └── M2vsM1/                          # M2-vs-M1 within each stimulation condition
│       ├── M2vsM1.U.DGE.csv
│       ├── M2vsM1.ADT.DGE.csv
│       ├── M2vsM1.NLRP3.DGE.csv
│       └── M2vsM1.Combo.DGE.csv
│
├── Pathway_Analysis/                    # clusterProfiler GO + KEGG GSEA per contrast
│   ├── M1/<contrast>/GO/                # Cnetplot, Dotplot, Emmap, Heatplot, GSEA_Enrichment.csv
│   ├── M1/<contrast>/KEGG/              # Same plot suite + Enrichment.csv
│   └── M1/<contrast>/KEGG_Pathview/     # Per-pathway Pathview diagrams (mmu*.pathview.png)
│   └── M2/<contrast>/...                # Mirror structure for M2 contrasts
│
├── Reviewer_Comment/                    # Analyses added during peer review
│   ├── Module_Analysis/
│   │   ├── GSEA_Plots/                  # GSEA enrichment plots per contrast
│   │   ├── GSEA_Results/                # Saved GSEA result objects (.rds)
│   │   └── GSEA_Summary/                # Per-contrast and combined CSV summaries
│   └── Raw/                             # Additional heatmaps (raw + row-scaled)
│
├── Heatmap_U.png, Heatmap_M1.png, Heatmap_M2.png   # Sample-level heatmaps
├── Bulk_RNAseq_UChicago_NLRP3.Rproj
└── .gitignore
```

---

## Pipeline overview

The full analysis lives in [`DGE_and_Pathway_Analysis.R`](DGE_and_Pathway_Analysis.R) and proceeds in five stages.

### 1. Read counting (upstream of this repo)
FASTQs aligned with **Subread** to the mouse genome; gene-level counts produced with **featureCounts** and merged into `Raw_Data/raw_combined_featureCounts.txt` (24 columns, one per BAM).

### 2. Differential gene expression — DESeq2
A `DESeqDataSet` is built from the featureCounts matrix using `~ condition` (or polarization-stratified) designs. Six within-polarization contrasts are computed for each of M1 and M2:

```
{ADT, NLRP3, Combo} vs U
NLRP3 vs ADT
Combo vs ADT
Combo vs NLRP3
```

plus four **M2-vs-M1** contrasts (one at each stimulation level: U, ADT, NLRP3, Combo). Each contrast is exported to `DGE/.../<contrast>.DGE.csv`.

### 3. Pathway / functional enrichment — clusterProfiler
For each contrast, ranked gene lists drive **GSEA** against:
- **GO** (Biological Process) via `gseGO` with `org.Mm.eg.db`
- **KEGG** via `gseKEGG` (organism `mmu`)

Plots (Dotplot, Cnetplot, Emmap/enrichment-map, Heatplot) and enrichment CSVs are written to `Pathway_Analysis/<polarization>/<contrast>/{GO,KEGG}/`.

### 4. KEGG pathway visualization — pathview
For top KEGG hits per contrast, **`pathview`** renders gene-level log2FC onto the canonical KEGG pathway diagram (e.g. `mmu04621` NOD-like receptor signaling, `mmu04668` TNF signaling, `mmu05171` Coronavirus disease, etc.) in `KEGG_Pathview/`.

### 5. Reviewer-driven re-analysis — module GSEA
Added in revision: a complementary GSEA pass against curated gene-set modules (`Reviewer_Comment/Module_Analysis/`), with per-contrast plots, saved `GSEA_*.rds` objects, and a combined summary table `GSEA_Summary_All_Comparisons.csv` plus mean log2FC roll-ups.

---

## Dependencies

R (≥ 4.2) with the following packages:

**Differential expression & data wrangling**
- [`DESeq2`](https://bioconductor.org/packages/DESeq2/) — count modeling and DE testing
- `dplyr`, `data.table`, `magrittr` — tabular wrangling

**Annotation**
- [`org.Mm.eg.db`](https://bioconductor.org/packages/org.Mm.eg.db/) — mouse gene annotation
- [`org.Hs.eg.db`](https://bioconductor.org/packages/org.Hs.eg.db/), [`Orthology.eg.db`](https://bioconductor.org/packages/Orthology.eg.db/) — for any human-mouse ortholog mapping
- [`biomaRt`](https://bioconductor.org/packages/biomaRt/) — gene ID conversion

**Functional enrichment**
- [`clusterProfiler`](https://yulab-smu.top/biomedical-knowledge-mining-book/) — GO / KEGG GSEA
- [`enrichplot`](https://bioconductor.org/packages/enrichplot/), [`DOSE`](https://bioconductor.org/packages/DOSE/) — enrichment visualization
- [`pathview`](https://bioconductor.org/packages/pathview/) — KEGG pathway gene-level overlays

**Visualization**
- [`ComplexHeatmap`](https://bioconductor.org/packages/ComplexHeatmap/), `circlize` — sample / DEG heatmaps
- `ggplot2` — general plotting

---

## Reproducing the analysis

1. Clone the repo and open `Bulk_RNAseq_UChicago_NLRP3.Rproj` in RStudio.
2. Update the path at the top of `DGE_and_Pathway_Analysis.R`:
   ```r
   setwd("…/Bulk_RNAseq_UChicago_NLRP3")
   load.path <- "…/Bulk_RNAseq_UChicago_NLRP3/Raw_Data/"
   ```
3. Ensure `Raw_Data/raw_combined_featureCounts.txt` is in place.
4. Source the script — it will (re)generate `DGE/`, `Pathway_Analysis/`, and the heatmaps.

> Raw FASTQ files are not committed. *(Add GEO/SRA accession or data availability statement here.)*

---

## Key contrasts at a glance

| Question | Contrast(s) |
|---|---|
| Effect of each stimulus on M1 macrophages | `M1.{ADT,NLRP3,Combo}.vs.U` |
| Effect of each stimulus on M2 macrophages | `M2.{ADT,NLRP3,Combo}.vs.U` |
| Does Combo do more than NLRP3 alone? | `{M1,M2}.Combo.vs.NLRP3` |
| Does Combo do more than ADT alone? | `{M1,M2}.Combo.vs.ADT` |
| M1 vs M2 baseline transcriptome | `M2vsM1.U` |
| Polarization-dependent stimulus response | `M2vsM1.{ADT,NLRP3,Combo}` |

---

## Citation

Manuscript in Preperation

## License

*(Add a license — e.g. MIT for code, CC-BY for figures — or state "All rights reserved" if you'd rather decide later.)*
