# ep_cellchat_analysis

Re-analysis of the paediatric ependymoma scRNA-seq cohort **GSE125969** (Gillen et al. 2020, *Cell Reports*) using Seurat v5, CellChat v2, and Hallmark module scoring, framed as an independent characterisation of the comparator TME used in Andrade et al. 2024 (*Nat Commun*).

## Scientific framing

Andrade et al. 2024 established that the K27M H3-mutant paediatric HGG tumour microenvironment is immunosuppressive, lymphoid-depleted, and myeloid-polarised for proliferation / metabolism. Their specificity claim rests on a brief comparison to the ependymoma cohort of Gillen et al. 2020. This project re-characterises the ependymoma side of that comparison at higher resolution: cross-patient integration, ligand–receptor inference focused on tumour-cell communication, and cohort-wide Hallmark programme scoring. The output is a testable biological narrative about what ependymoma's TME actually looks like (structurally stromal, with directionally inverted MIF signalling centred on the hypoxic MEC state) and which therapeutic claims survive the contrast.

The headline biological findings are written up in [`reports/04-cellchat-biological-interpretation.md`](reports/04-cellchat-biological-interpretation.md) (not tracked, see below).

## Pipeline

Each analysis step is a separate Quarto notebook (`.qmd`) under `scripts/`, passing state via `.rds` checkpoints.

| Step | Notebook | Role |
|---|---|---|
| 01 | [`scripts/01-data-import-qc.qmd`](scripts/01-data-import-qc.qmd) | Load the GSE125969 count matrix and metadata. Derive `patient_id` from cell barcodes. Attach `tumor_subtype` and anatomical `location` as biological covariates. Hard-floor QC filtration on `nFeature_RNA`, `nCount_RNA`, `percent.mt`. |
| 02 | [`scripts/02-normalize-integrate.qmd`](scripts/02-normalize-integrate.qmd) | Seurat v5 layered workflow: split `RNA` by `patient_id`, `LogNormalize`, top-2000 HVFs, `ScaleData`, `RunPCA(npcs=50)`. Unintegrated baseline UMAP. Harmony integration on `patient_id`. Louvain clustering, Harmony UMAP. `JoinLayers()` for downstream DE. |
| 03 | [`scripts/03-annotate-shinycell.qmd`](scripts/03-annotate-shinycell.qmd) | Refactor the 19 shipped Gillen `cell_type` labels into a 12-class `cell_class` column (3 non-tumour, 5 PFA neoplastic states, 4 subtype-specific tumour bins). Drop the shipped `Doublets` cluster. Build a ShinyCell2 browser app. |
| 04 | [`scripts/04-cellchat-tumor.qmd`](scripts/04-cellchat-tumor.qmd) | CellChat v2 inference (full `CellChatDB.human`, triMean, `nboot = 100`, `min.cells = 10`). Aggregate network, information-flow ranking, tumour-incoming + tumour-outgoing bubbles, per-pathway chord plots for the six Andrade-contrast axes, signalling-role analyses. 83 pathways and 3300 significant L-R interactions retained on this cohort. |
| 05 | [`scripts/05-hallmark-module-scores.qmd`](scripts/05-hallmark-module-scores.qmd) | MSigDB Hallmark pathways via `msigdbr`, refactored into a named list and scored per cell with `Seurat::AddModuleScore`. Cohort-wide z-scored heatmap, inflammation / metabolism / proliferation violin panels, tumour-state hypoxia ridge, hypoxia UMAP. Rebuilt ShinyCell app with scores exposed as metadata. |

Notebook 04's expensive `computeCommunProb` step is also mirrored into a set of small `.R` driver scripts in `scripts/runners/` (git-ignored; local driver code for shell execution with per-stage logging and dual PDF + PNG outputs). Those are convenience scaffolding, not part of the tracked analysis.

## Directory layout

```
ep_cellchat_analysis/
├── ep-cellchat-analysis.Rproj      # RStudio project file (at root)
├── .gitignore                      # Whitelist approach, only scripts/ and context/ tracked
├── README.md                       # This file
├── scripts/                        # All analysis code, tracked
│   ├── 01-data-import-qc.qmd
│   ├── 02-normalize-integrate.qmd
│   ├── 03-annotate-shinycell.qmd
│   ├── 04-cellchat-tumor.qmd
│   ├── 05-hallmark-module-scores.qmd
│   └── runners/                    # Local .R runners for 04, git-ignored
├── context/                        # Brief + data context, tracked
│   ├── data_context.md
│   └── gse125969_analysis_brief.md
├── read/                           # Raw input data (GSE125969 TSVs), NOT tracked
├── checkpoints/                    # Intermediate .rds files, NOT tracked
├── write/                          # Generated outputs, NOT tracked
│   ├── figures/                    #   PDFs + PNGs for every plot
│   ├── tables/                     #   CSVs exported from each step
│   └── apps/                       #   ShinyCell2 app directories
└── reports/                        # Interpretation drafts, NOT tracked
    └── 04-cellchat-biological-interpretation.md
```

## What is tracked

- `scripts/**` (notebooks), excluding `scripts/runners/`
- `context/**` (analysis brief, data context notes)
- `README.md`, `.gitignore`, `ep-cellchat-analysis.Rproj`
- `.gitkeep` placeholders inside `read/`, `checkpoints/`, `write/figures/`, `write/tables/` (so the expected directory layout is visible in a fresh clone without committing the data itself)

## What is NOT tracked

- Raw data under `read/` (GSE125969 count matrix is ~875 MB)
- Intermediate `.rds` checkpoints under `checkpoints/`
- All generated figures, tables, and ShinyCell apps under `write/`
- Interpretation drafts under `reports/`
- Local runner scripts under `scripts/runners/`

## Replication

```bash
# 1. Clone and enter the project
git clone <remote-url> ep_cellchat_analysis
cd ep_cellchat_analysis

# 2. Pull the GSE125969 files into read/
#    Both files are linked from GSE125969 on GEO.
mkdir -p read
curl -o read/gse125969_count_matrix.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125969/suppl/GSE125969_count_matrix.tsv.gz"
curl -o read/gse125969_cell_metadata.tsv.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125969/suppl/GSE125969_cell_metadata.tsv.gz"
gunzip read/*.gz

# 3. Open the project in RStudio via ep-cellchat-analysis.Rproj
#    Then run the notebooks in order (01 -> 05).
```

Notebook 04's `computeCommunProb` takes ~20–30 minutes on this cohort; everything else completes in minutes on a laptop.

## Dependencies

**Core single-cell stack**

- `Seurat` (v5+), `harmony`, `tidyverse`, `readr`, `Matrix`, `patchwork`

**Communication inference**

- `CellChat` v2 (GitHub: `jinworks/CellChat`), `ComplexHeatmap`, `presto` (GitHub: `immunogenomics/presto`), `future`

**Hallmark scoring**

- `msigdbr` (CRAN)

**Visualisation**

- `BadranSeq` (GitHub: `wolf5996/BadranSeq`) for `EnhancedElbowPlot`, `do_DimPlot`, `do_FeaturePlot`
- `circlize` for heatmap colour ramps

**Interactive browser**

- `ShinyCell2` (GitHub: `the-ouyang-lab/ShinyCell2`)

Installation notes:

- `CellChat` depends on `ComplexHeatmap` from Bioconductor, which `install_github` does not resolve automatically. Install via `BiocManager::install("ComplexHeatmap")` first.
- `presto` is a required (not optional) runtime dependency for CellChat's Wilcoxon step in v2. Install from GitHub.
- `rlang >= 1.2.0` is required by the current stack; if the interactive session has an older version loaded, install in a fresh Rscript process and restart R.

## Code style

- Every QMD chunk is independently runnable: declares its own `library()` calls, loads its input from a checkpoint file, writes its output to a checkpoint or figure / table path.
- All paths are rooted at the project (e.g. `read/...`, `checkpoints/...`, `write/figures/...`); no `../` hops.
- Metadata operations on Seurat objects are done by pulling `seu@meta.data` as a tibble, transforming with `dplyr`, then writing back to the Seurat object. `tidyseurat` is deliberately avoided.
- Plots use `BadranSeq::do_*` and `EnhancedElbowPlot` first, Seurat / SCpubr only as a fallback where BadranSeq does not cover the plot type.
- PDFs are written with `device = cairo_pdf` to handle Unicode cleanly.
- Comments in code are used sparingly and only to explain the *why* of a non-obvious step; `# Section ----------` headers organise each chunk into Libraries / Inputs / Processing / Outputs.

## Author

Badran Elshenawy, University of Oxford (2026).

## Citation

The underlying data is Gillen et al. 2020, *Cell Reports* (PMID: 32783945). The comparator framing is Andrade et al. 2024, *Nature Communications* (DOI: 10.1038/s41467-024-52096-w).
