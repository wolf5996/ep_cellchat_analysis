# GSE125969 Analysis Brief — Pathania Lab Interview Demo

## Purpose

This analysis is a targeted technical demonstration for a meeting tomorrow with Dr Manav Pathania (PI, Pathania Lab, Ludwig Institute for Cancer Research, Oxford). Context: interview-stage evaluation for a Senior Computational Postdoctoral Research Scientist position (Ref 185648). Manav is a pure experimentalist who is hiring his first in-house senior computational postdoc and has historically relied on the Kleinman lab at McGill for computational work. He is evaluating whether the candidate can extract biology from data and generate bench-testable hypotheses, not whether the candidate can write methodologically novel code.

The deliverables are biology-forward. Methods are scaffolding. Every figure and conclusion must be reportable in a ten-minute conversation where the audience cares about what the data says, not how you processed it.

## Scientific framing

Andrade et al. 2024 (*Nature Communications*, DOI 10.1038/s41467-024-52096-w) characterise the immune tumour microenvironment of H3-mutant paediatric gliomas using scRNA-seq and IMC, concluding that these tumours are overwhelmingly myeloid-infiltrated, lymphoid-depleted, and transcriptionally immunosuppressive. Ependymomas (posterior fossa group A and supratentorial) serve as the comparator cohort in that paper, showing mixed myeloid/lymphoid infiltration — a pattern the authors use to argue that the H3-mutant phenotype is specific to HGG, not a generic feature of paediatric brain tumours.

GSE125969 is that ependymoma cohort. Analysing it directly lets us re-establish the comparator phenotype independently, using modern tools (CellChat v2 rather than the CellPhoneDB v2 used in Andrade, plus systematic Hallmark enrichment). The analytical contribution is showing that the ependymoma TME has a fundamentally different cell-cell communication structure and myeloid polarisation profile from the K27M TME described in Andrade figures 2 and 4. This validates the Andrade paper's specificity claim using the same raw data, re-analysed with a more systematic approach.

This framing matters. If Manav asks what the analysis is for, the answer is: *re-demonstrating the comparator phenotype that anchors the paper's main specificity claim, using tools that produce more systematic and bench-testable outputs than the original analysis.*

## The three analytical pillars

Each pillar is one of the three skills to be visibly demonstrated. Each must produce at least one presentable figure.

### Pillar 1 — Data integration

GSE125969 contains multiple ependymoma samples (verify sample count and metadata at download). Demonstrate cross-sample integration using the standard Seurat v5 + Harmony workflow:

- Per-sample QC (mitochondrial content, nUMI, nFeature filtering)
- Normalisation (LogNormalize) and variable feature selection
- PCA on top 2000 variable genes
- Harmony integration with Sample as the batch variable
- UMAP embedding on Harmony-corrected PCs
- Shared nearest neighbour clustering (Louvain, k=20)

The explicit demonstration here is that heterogeneous samples can be harmonised into a single analytical space without batch effects dominating the clustering. Include a pre-integration vs post-integration UMAP side-by-side as one of the figures — this is the standard way to visually demonstrate integration success.

If GSE125969 contains both supratentorial and posterior fossa ependymomas (it should — Andrade mentions 20 PFA and 2 hemispheric), treat location as a biological covariate rather than a batch effect. Preserve biological signal, remove technical signal.

### Pillar 2 — CellChat v2 for cell-cell communication

Use CellChat v2 (Jin et al. 2024, *Nature Protocols*) not v1. The v2 package handles multi-subunit ligand-receptor complexes and has expanded CellChatDB coverage. Run on the fully integrated and annotated dataset.

Specific requirements:

- Use the human CellChatDB with all three interaction types (Secreted Signalling, ECM-Receptor, Cell-Cell Contact)
- Default parameters: 1000 permutations, trimean for expression averaging
- Compute communication probability at both the ligand-receptor pair level and pathway level
- Identify the top signalling pathways ranked by information flow

Biological targets to extract from the CellChat output:

- Which ligand-receptor pairs mediate myeloid-to-tumour communication in ependymoma?
- Are the immunosuppressive axes identified by Andrade in K27M (OSM-LIFR, LGALS9-CD47, MIF-CD74) present or absent in ependymoma?
- What communication axes involving lymphoid cells appear in ependymoma but should be absent in K27M (e.g., chemokine signalling from myeloid to T cells)?

The contrast with Andrade's K27M findings is the biological story. Ependymoma should show active myeloid-to-lymphoid chemokine signalling (CCL3, CCL4, CXCL families) and functional antigen presentation, consistent with its lymphoid-permissive TME. K27M does not. Frame the outputs to make this contrast visible.

Produce these figures:

- Aggregate communication network showing overall signalling strength between cell types (circle plot or chord diagram)
- Signalling pathway rankings with top 15-20 pathways
- Specific ligand-receptor pair heatmap for myeloid-to-lymphoid and myeloid-to-tumour communication

### Pillar 3 — Hallmark signature scoring on myeloid cells

Use single-sample GSEA (ssGSEA) via the `escape` R package on the MSigDB Hallmark gene set collection (v7.5.1 or current). Score at the individual cell level, then aggregate by cell type cluster.

Focus the interpretation on myeloid subsets specifically, mirroring the approach used in Andrade figure 2. Key Hallmark pathways to highlight:

- Inflammatory response
- TNF-α signalling via NF-κB
- Interferon gamma response
- IL6 JAK STAT3 signalling
- Complement
- Hypoxia
- Glycolysis
- Oxidative phosphorylation

The Andrade paper showed that K27M microglia are enriched for proliferation (G2M checkpoint, mitotic spindle) and metabolic (glycolysis) pathways while lacking inflammatory and antigen-presentation pathways. The comparator prediction for ependymoma myeloid cells is the opposite — enrichment for inflammatory and cytokine signalling pathways relative to proliferation/metabolism. Produce the figure in a format that makes this contrast readable at a glance (e.g., a heatmap of NES values with cell types as columns and pathways as rows, split by myeloid subtype).

## Technical specifications

**Environment.** R ≥ 4.3 with Seurat v5, harmony, CellChat v2 (github: jinworks/CellChat), escape, msigdbr, SeuratDisk for data conversion as needed. Python fallback acceptable for data download and format conversion only — the analysis itself should stay in R for CellChat compatibility.

**Data acquisition.** GSE125969 is on GEO. Use `GEOquery` to pull the series metadata, then download the count matrices. If the data is in 10x format (matrix.mtx, barcodes.tsv, features.tsv), use `Read10X`. If it's in h5 format, use `Read10X_h5`. Verify sample count before proceeding — the Andrade paper references up to 22 ependymoma samples, though the scRNA-seq subset used was smaller.

**Cell type annotation.** Two options in order of preference:

1. Reference-based annotation using a published brain immune cell reference (e.g., the Pombo Antunes et al. 2021 *Nature Neuroscience* reference used by Andrade, accession GSE163120) via SingleR or Seurat label transfer.
2. Marker-based manual annotation using canonical markers: PTPRC (CD45) for immune, CSF1R/CX3CR1/P2RY12/TMEM119 for microglia, CD68/CD163 for macrophages, CD3D/CD3E/CD8A/CD4 for T cells, CD19/MS4A1 for B cells, GFAP/OLIG2 for glial/tumour.

If time is tight, option 2 is faster and more defensible in conversation because you can explain every marker. Reference-based annotation is more rigorous but requires the reference to download cleanly.

**Computational resources.** If running on a single machine, CellChat on ~50k cells with 1000 permutations takes ~30 minutes. Harmony integration is fast. Plan for 4-6 hours of total wall-clock time end-to-end.

## Deliverables

Produce these outputs in `/outputs` directory:

1. `figure_1_integration.pdf` — Pre- vs post-Harmony UMAPs, coloured by sample and cell type. Two panels side-by-side.
2. `figure_2_celltype_composition.pdf` — Stacked bar chart of immune cell composition per sample, plus a summary bar for the cohort average. Compare against Andrade figure 1D's ependymoma proportions where possible.
3. `figure_3_cellchat_network.pdf` — Aggregate cell-cell communication network (circle plot) showing interaction strength between all cell types.
4. `figure_4_cellchat_pathways.pdf` — Top 15-20 signalling pathways ranked by information flow, with the axes relevant to the Andrade contrast highlighted.
5. `figure_5_hallmark_myeloid.pdf` — Heatmap of Hallmark NES values across myeloid cell subtypes, with pathways of interest (inflammatory, glycolysis, IFN) visually prominent.
6. `analysis_summary.md` — A 500-1000 word writeup structured as:
   - One-paragraph framing (what was analysed and why it matters for the Andrade paper's claims)
   - Key findings with figure references (not methods descriptions)
   - Biological interpretation (what does this tell us about ependymoma vs K27M TME)
   - Testable hypotheses (what bench experiments would this analysis suggest)
7. `analysis_script.R` — The reproducible pipeline. Clean, commented, runnable end-to-end.

The summary writeup is the single most important deliverable. If anything has to be cut for time, cut the figures. The writeup is what gets shown to Manav.

## Scope and failure modes

**Time budget.** Tonight only. Realistically 4-6 hours of compute and iteration. Do not attempt to build novel methodology. Do not chase edge cases in the data.

**Minimum viable deliverable.** If only one pillar completes, that pillar must be CellChat (Pillar 2) because it's the most analytically distinctive and the most directly aligned with the Andrade paper's central framework. Integration (Pillar 1) can be simplified to per-sample analysis if integration fails. Hallmark scoring (Pillar 3) is the cleanest to produce and should work even if the other two have issues.

**Failure modes to avoid.**

- Do not run CellChat at the broad-category level (all myeloid as one cluster) — it loses the biological nuance the analysis is for. Stratify myeloid into at least microglia, macrophages, and monocytes.
- Do not over-interpret rare populations. If a cell type has <30 cells, note it but don't build conclusions on it.
- Do not produce methodologically virtuosic figures that take time away from interpretation. A simple heatmap with good annotations beats a fancy alluvial plot.
- Do not reproduce the Andrade paper's analysis in full — that's not the point. The point is to extend it with CellChat v2 and Hallmark scoring on the ependymoma data that Andrade used as a comparator but didn't deeply characterise.

## Framing for the conversation

If Manav asks about this analysis, the opening line is: *"I re-analysed the ependymoma scRNA-seq cohort from Andrade 2024 with CellChat v2 and Hallmark scoring to independently characterise the comparator TME. The analysis confirms that ependymomas have active myeloid-to-lymphoid chemokine signalling and inflammatory programme engagement that the K27M tumours lack — which supports the paper's specificity claim but also identifies specific signalling axes that could be experimentally tested for reversal of the K27M phenotype."*

That framing does three things: credits the paper, extends it meaningfully, and hands Manav a bench-testable hypothesis. Those are the three things he's hiring for.

Do not oversell. If results are partial, say they're partial. If the analysis raises questions rather than answers them, present the questions. Manav has been through peer review on his own papers and will read hedged, calibrated claims as more credible than unhedged ones.
