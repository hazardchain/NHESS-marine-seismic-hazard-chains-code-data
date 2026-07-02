# Code and data — Marine seismic hazard chains (NHESS)

Reproducibility materials for the paper *Quantitative Modelling of Marine Seismic
Hazard Chains: Topology, Core Amplification Nodes, and Implications* (Shi, Zhang,
Hu; submitted to *Natural Hazards and Earth System Sciences*).

Manuscript (LaTeX) repository: <https://github.com/Chunyu-Hugh/NHESS-marine-seismic-hazard-chains>

**Interactive network viewer (Figs 1 & 2):** <http://chunyuhu.com/NHESS-marine-seismic-hazard-chains-code-data/>
— click any node to trace its couplings; Chinese / English toggle. Source in [`docs/`](docs/).

## Contents

| Folder | What it holds |
|--------|---------------|
| `data/` | The documented cascade tables that define the hazard networks — marine (`data/marine/`) and mainland China (`data/mainland_china/`). See `data/README.md`. |
| `code/` | R scripts that build the networks, compute all reported metrics/figures, and run the robustness analysis. See `code/README.md`. |
| `code/results/` | Example outputs of the robustness script (summary + Fig. 3). |
| `docs/` | Interactive web viewer of the networks (served via GitHub Pages). |

## Reproduce

Requires R (≥4.x) with: `docxtractr, igraph, dplyr, expm, Matrix, ggplot2, reshape2, gridExtra, svglite`.

```r
cd code
# 1) Marine 39-node network + Fig. 1 + A/M heatmaps
#    (set path1 <- "../data/marine/tsunami_Dohmen_2025_EQ.docx" in the script)
Rscript multi_hazard_China_sy.R
# 2) Robustness / sampling-sufficiency (Sect. 3.4 + Fig. 3)
Rscript sensitivity_crossvalidation.R
```

Verified (R 4.5.3 + igraph) to reproduce the paper's **39 nodes / 152 edges** and
every Table 2 and Table 3 value; the four `mainland_china/` files reproduce the
15-class continental network (Table 4).

## License and citation

*TODO(authors): add a license (e.g. CC-BY-4.0 for the data, MIT for the code) and,
after archiving this repository to Zenodo, insert the dataset DOI here and in the
manuscript's Code-and-data-availability statement.*
