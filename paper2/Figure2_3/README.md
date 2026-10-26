# Transcription Factor Analysis — Tabula Sapiens 2.0

This folder holds the Jupyter notebooks for the transcription factor (TF) analyses in:

> The Tabula Sapiens Consortium. **Tabula Sapiens 2.0: A comprehensive transcriptomic atlas of human cell types.**
> *Cell* 189 (2026). https://doi.org/10.1016/j.cell.2026.08.010

The notebooks measure how cell-type specific each of the **1,635 expressed human TFs** is. They use the τ (tau) statistic on
mean log-normalized expression across **175 cell types from 28 tissues** (10x droplet data only). The outputs are:

- **Figure 2:** TF × cell type expression heatmap and broad cell-type insets.
- **Figure 3:** GO enrichment of the 745 non-specific TFs.
- **Supplementary Figures 8–11:** donor reproducibility of τ, Human Protein Atlas (HPA) protein validation, and three special TF groups.

For the rest of the Tabula Sapiens 2.0 code and data, see the [top-level repository README](../README.md).

---

## Contents

| Notebook | Produces (paper figure / table) |
|---|---|
| `TFonTSMean_09-02-26.ipynb` | **Step 0 – preprocessing.** Reads the full TS2.0 h5ad. Writes TF × cell type, TF × broad cell class and TF × tissue mean log-normalized expression tables, the percent of non-zero cells, and the cleaned metadata (all in `TFonly/`). |
| `TF_Fig2_09-02-26.ipynb` | Computes τ for every TF (`TauValuesCelltype.csv`, `LowestTauTFs100.csv`). Makes **Fig. 2A** (heatmap + `Fig2A_ColorKey.pdf`), **Fig. 2B–F** insets (T cells, fibroblasts, stem cells, endothelial cells, male germ cells) and **Fig. S9A** (τ histogram). Also writes `CellTypeColorPalettes.csv` and `BroadCellTypeColorPalette.csv`. |
| `DonorCorrelegrams_10-19-26.ipynb` | **Fig. S8**: Pearson correlation of per-donor τ for donors TSP2, TSP14, TSP21, TSP25 and TSP27. |
| `V9GoTermAnalysisOnUbiq.ipynb` | GSEApy/Enrichr GO enrichment (`GO_Biological_Process_2021`) of the non-specific TFs (τ < 0.85), using all TFs as background. Output goes to **Table S5** (files in `GoAnalysisV9/`). |
| `Fig3_BubbleChart_02-13-26.ipynb` | **Fig. 3**: GO-term bubble chart (69 terms, adj. p < 0.02, grouped into 4 broad categories). Also writes cell-type / compartment count tables to `Counts/`. |
| `ProteinTF-Tiss_10-19-26.ipynb` | Matches TS2.0 TFs to HPA immunohistochemistry by tissue (`TFproDetTiss.csv`). Makes **Fig. S9B** and **Fig. S9D**. |
| `ProteinTF-CellType_10-19-26.ipynb` | Matches TS2.0 TFs to HPA detection by cell type (`TFproDetCT.csv`). Makes **Fig. S9C** and **Fig. S9E**. |
| `ProteinTF_SubCell_10-19-26.ipynb` | **Fig. S9F**: nuclear localization of the 100 lowest-τ TFs across tissues. |
| `SuppFig10_ThreePartChart.ipynb` | **Fig. S10**: heatmap of 587 TFs in three groups (core essential, no known motif, low-τ zinc finger). Also writes `Fig2A_Tau_Distribution.pdf`. |
| `SuppFig11_Tau_Ubiq_Distributions.ipynb` | **Fig. S11**: τ histograms for the same three TF groups. |

The two PDFs in this folder are the manuscript (`TS2.0_online_090226.pdf`) and the supplementary information
(`TabulaSapiens2_Supplemental_Information.pdf`), included for reference.

## Run order

```
TFonTSMean  ──►  TF_Fig2  (writes TauValuesCelltype.csv)
                   │
                   ├──► DonorCorrelegrams                      (Fig S8)
                   ├──► SuppFig10_ThreePartChart               (Fig S10)
                   ├──► SuppFig11_Tau_Ubiq_Distributions       (Fig S11)
                   ├──► ProteinTF-Tiss ──► ProteinTF_SubCell   (Fig S9B, S9D ─► S9F)
                   ├──► ProteinTF-CellType                     (Fig S9C, S9E)
                   └──► V9GoTermAnalysisOnUbiq ──► [manual curation] ──► Fig3_BubbleChart   (Table S5 ─► Fig 3)
```

Every downstream notebook needs only the small CSVs in `data/` (below), so you can run any of them without first re-running
`TFonTSMean`, which needs the full ~1.1M-cell h5ad.

## Methods summary

**Expression tables.** Only 10x droplet cells are used, so that all organs are treated the same way. For each TF, the notebooks compute the mean
log-normalized expression per cell type (175 of the 182 annotated cell types), per broad cell class (38) and per tissue (28).
The TFs come from the University of Toronto *Human Transcription Factors* list (https://humantfs.ccbr.utoronto.ca/allTFs.php;
Lambert et al., *Cell* 2018), mapped by Ensembl ID. Of the 1,639 TFs, 1,635 are expressed; SHOX and ZBED1 have zero expression.

**τ specificity.** For each TF, the expression vector *x* over *N* cell types is divided by its maximum. Then:

$$\tau = \frac{\sum_{i=1}^{N} (1 - x_i)}{N - 1}$$

The paper cites `tspex` v0.6.3 (`tspex.TissueSpecificity(TF_exp, 'tau', log=False)`). The notebooks compute the same formula inline
(see `TF_Fig2_09-02-26.ipynb`), so `tspex` is **not** a dependency. TFs with τ > 0.85 count as cell-type specific (890 TFs);
TFs with τ ≤ 0.85 count as non-specific (745 TFs).

**Heatmaps.** These use seaborn. In the full-size heatmaps, `vmax` is set to the 99.5th percentile; the Fig. 2B–F insets are not clipped.

**Protein validation.** RNA-detected TFs are compared with HPA immunohistochemistry detection by tissue, cell type and
subcellular location. The analysis uses 699 TFs common to both datasets. TF–tissue / TF–cell-type pairs with "uncertain" antibody staining are
excluded, and so is thymus.

**GO enrichment.** This uses `gseapy.enrichr` with `GO_Biological_Process_2021`, the 745 non-specific TFs as the query and all
1,635 TFs as background. The 69 terms with adjusted p < 0.02 were grouped **by hand** into broad functional categories for Fig. 3.

## Data

Every notebook expects its inputs in one data directory, laid out like this. Files marked **(repo)** are included in `data/`.
Files marked **(external)** must be downloaded.

```
data/
├── HumanTranscriptionFactorsEnsembl.csv          (repo)  TF list with Ensembl IDs
├── TauValuesCelltype.csv                         (repo)  τ per TF  [output of TF_Fig2]
├── LowestTauTFs100.csv                           (repo)  100 lowest-τ TFs  [output of TF_Fig2]
├── CellTypeColorPalettes.csv                     (repo)  cell-type colors  [output of TF_Fig2]
├── BroadCellTypeColorPalette.csv                 (repo)  broad cell-type colors  [output of TF_Fig2]
├── FullListForKonstantin_KK.csv                  (repo)  HPA tissue / cell type / subcellular table (see Acknowledgments)
├── TFproDetTiss.csv, TFproDetCT.csv              (repo)  standardized HPA detection tables  [outputs of ProteinTF-*]
├── TFonly/
│   ├── TSP1_30_by_CellType_Full_Version2d_Nov2024_LogNorm_10XOnly_Mean_TFonly.csv        (repo)
│   ├── TSP1_30_by_BroadCellClass_Full_Version2d_Nov2024_LogNorm_10XOnly_Mean_TFonly.csv  (repo)
│   ├── TSP1_30_by_Tissue_Full_Version2d_Nov2024_LogNorm_10XOnly_Mean_TFonly.csv          (repo)
│   ├── TSP1_30_by_PercentNZ_Full_Version2d_Nov2024_LogNorm_10XOnly_Mean_TFonly.csv       (repo)
│   ├── TSP1_30_by_CellType_TSP{2,14,21,25,27}_Version2d_Nov2024_LogNorm_10XOnly_Mean_TFonly.csv  (repo)  per-donor means for Fig S8
│   └── TSP1_30_version2d_10X_SS_2024_11_12_metadata_CLEAN_TFonly.csv                     (external, ~460 MB; output of TFonTSMean)
├── ClaudeFiles/                                  (repo)  TF group lists for Fig S10 / S11
│   ├── 49_core_essential_TS2_TF.csv              (n = 49)
│   ├── TFs_unknown_motif_codebook.csv            (n = 214)
│   ├── all_zinc_finger_TFs.csv                   (n = 763)
│   └── low_tau_zinc_finger_TFs.csv               (n = 418)
├── GoAnalysisV9/                                 (repo)  GO results (Table S5) and curated term groupings for Fig 3
│   ├── RawTable5V9_UbiqGenesGOBiologicalProcesses2021.csv
│   ├── SuppTable5_V9_UbiqGenesGOBiologicalProcesses2021.xlsx
│   ├── SuppTable5V9_GO_TopPage.csv
│   ├── ShortGoTermsUnedited.csv                  written by Fig3 notebook for hand editing
│   └── ShortGoTermsEdited.csv                    hand-curated version read back by Fig3 notebook
└── Counts/                                       (repo)  cell-type / compartment counts  [outputs of Fig3 notebook]
```

**External data.** These are needed only for `TFonTSMean` and the count tables in `Fig3_BubbleChart`:

- `TSP1_30_min200_2500_decontx_scvi_donorassay_version2d_20241126.h5ad` and
  `TSP1_30_metadata_min200_2500_decontx_scvi_donorassay_version2d_20241112.csv`. These are the Tabula Sapiens 2.0 full object and its
  metadata, available from [figshare](https://figshare.com/articles/dataset/Tabula_Sapiens_v2/27921984),
  [CELLxGENE](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) and GEO
  ([GSE306755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306755)).
- The ~460 MB `..._metadata_CLEAN_TFonly.csv` is too large for GitHub. Recreate it by running `TFonTSMean`.

## Running the notebooks

1. Clone the repository and set up a Python environment with the packages the notebooks import.
2. **Edit the path variables in the first code cells of each notebook.** The notebooks hard-code the authors' cluster
   paths. Set `DATA_IN_DIR`, `DATA_OUT_DIR`, `DATA_IN_MEANS_DIR` (= `data/TFonly/`), `SHARED_DATA_DIR` (location of the
   h5ad/metadata) and, in `SuppFig11`, `DATA_DIR`/`OUT_DIR` to point at your local `data/` folder.
   Paths are built by plain string concatenation, so keep the trailing `/`.
3. Run the notebooks in the order shown in [Run order](#run-order). Figures are saved as PDFs into the output directory. Several notebooks
   (the `ProteinTF-*` ones) write their figures to `DATA_IN_DIR`.

## Notes and known quirks

- **`TFonTSMean`:** the `read_h5ad` line is commented out because the object had already been loaded in the session. Uncomment it
  before running. If you set `ALL_GENES = True`, the notebook writes to an `AllGenes/` subfolder, which you must create first.
  The per-donor tables (`..._TSP2_...` etc.) used by `DonorCorrelegrams` come from the same code run on one donor at a time.
- **`V9GoTermAnalysisOnUbiq`:** this was originally run on a laptop. The `to_csv` call that saves the raw results is commented out to
  avoid overwriting the published table.
- **`Fig3_BubbleChart`:** the GO term short names and category groupings come from the hand-edited
  `GoAnalysisV9/ShortGoTermsEdited.csv`, so Fig. 3 is not fully regenerated from code.
- **`DonorCorrelegrams`:** this notebook was first written against `tspex` and then switched to the inline τ code from `TF_Fig2`.


