# manteltext

Analysis code for the framing text (*Manteltext*) of a cumulative dissertation. It re-analyses single-cell data from two papers on severe COVID-19:

- **Wendisch D**, Dietrich O, Mari T, *et al.* SARS-CoV-2 infection triggers profibrotic macrophage responses and lung fibrosis. *Cell* 184(26):6243–6261.e27 (2021). [doi:10.1016/j.cell.2021.11.033](https://doi.org/10.1016/j.cell.2021.11.033)
- Schulte-Schrepping J, Reusch N, Paclik D, *et al.* (incl. Wendisch D). Severe COVID-19 is marked by a dysregulated myeloid cell compartment. *Cell* 182(6):1419–1440.e23 (2020). [doi:10.1016/j.cell.2020.08.001](https://doi.org/10.1016/j.cell.2020.08.001)

> **Status: archived.** This is the analysis as it was for the thesis defense, kept for transparency. The input data are not included, so the notebooks will not run end-to-end without them.

## Collaboration

The analysis was developed together with [Oliver Dietrich](https://github.com/OliverDietrich), second author of the 2021 paper, and builds on the work of the [Saliba lab](https://github.com/saliba-lab) (Single Cell Analysis team, [Helmholtz Institute for RNA-based Infection Research](https://www.helmholtz-hiri.de/en/research/organisation/teams/team/single-cell-analysis/), Würzburg; senior author A.-E. Saliba). The original analysis code of the 2021 paper is published there: [saliba-lab/covid19-profibrotic-macrophages](https://github.com/saliba-lab/covid19-profibrotic-macrophages) (MIT license). A related collection of BAL samples from patients with severe COVID-19 is in [saliba-lab/covid19-bal-atlas-integration](https://github.com/saliba-lab/covid19-bal-atlas-integration).

## What is in here

| Folder | Contents |
|---|---|
| `scripts/` | Working notebooks and helper functions (`functions_bal_macrophages.R`, `functions_mono.R`, `initialize_bal_data.R`). Bronchoalveolar lavage (BAL) macrophage gene-set scoring and comparison with other lung fibrosis datasets, monocyte analysis, compositional analysis, and multi-epitope ligand cartography (MELC) imaging. |
| `ordered_scripts/` | The scripts arranged for the defense: figure notebooks, differential-expression models for the Schulte-Schrepping cohort (DESeq2 with severity, sex and donor terms), clustering checks, and a few of the resulting plots and tables. |
| `python/` | scCODA cell-composition analyses (Python code saved as `.txt`) and their result tables. |

`scripts/MELC_old.Rmd` is an earlier, separate MELC analysis. Some notebooks in `scripts/` are earlier versions of ones in `ordered_scripts/`.

## Data

Raw and processed input data are not part of this repository. They come from the two publications above (see each paper's data availability statement). For the 2021 paper, the Saliba lab repository above describes where the count matrices and Seurat objects are stored. Scripts expect them in a sibling folder two levels above the project:

```
<parent>/
├── raw_data/
│   ├── schulte_schrepping_et_al/
│   └── wendisch_et_al_final/
└── <projects>/manteltext/      # this repository
```

Intermediate results and figures are written to `intermediate_data/` and `output/` (git-ignored).

## Requirements

R with Seurat, tidyverse, tidyseurat, DESeq2, here, janitor, ggrepel, ggpubr, effsize, patchwork and readxl. Python with [scCODA](https://github.com/theislab/scCODA), pandas and matplotlib. No lockfile was kept for this project.

Open `manteltext_git.Rproj` so that `here::here()` resolves to the project root.
