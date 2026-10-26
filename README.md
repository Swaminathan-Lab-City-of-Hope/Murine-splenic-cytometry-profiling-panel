# Murine splenic spectral cytometry

Contains the R scripts used for FlowSOM clustering, UMAP visualization, and marker-expression analysis of murine splenic leukocytes, as presented in Figure 3 and Supplementary Figure 7 of *High-dimensional spectral cytometry for phenotyping murine splenic leukocytes in health and disease models*.

## Scripts

| Script | Analysis |
| --- | --- |
| `00_samples.R` | Sample metadata and acquisition batches |
| `01_preprocess.R` | FCS import, arcsinh transformation, and batch correction with cyCombine |
| `02_flowsom.R` | FlowSOM clustering and 25 metaclusters |
| `03_umap.R` | Cell subsampling, PCA, UMAP, cluster labels and colors (Figure 3A), and strain-level and individual-mouse UMAPs (Supplementary Figure 7) |
| `04_heatmap.R` | Marker-expression heatmap (Figure 3B) |

## Requirements

R 4.3.2 and the packages used in the analysis: `flowCore` 2.14.2, `Biobase` 2.62.0, `cyCombine` 0.3.0, `FlowSOM` 2.10.0, `uwot` 0.2.4, and `ggplot2` 4.0.3. Several scripts check package versions explicitly.

## Usage

Place the six live-cell-gated, spectrally unmixed FCS files in `data/fcs/`, using the filenames specified in `scripts/00_samples.R`. From the repository root, run:

```bash
Rscript run_pipeline.R
```

The scripts run in sequence and save intermediate results and figures to `results/`.

## Figures

- **Figure 3A:** `results/03B_K25_numbered_UMAP_spatial_palette.png`
- **Figure 3B:** `results/04_K25_marker_heatmap.png`
- **Supplementary Figure 7:** `results/Supplementary_Figure_7_UMAP.png`

## Data

The FCS files are not included in this repository. Manual gating and spectral unmixing, performed before the R analysis in FlowJo and SpectroFlo, are not covered by these scripts.
