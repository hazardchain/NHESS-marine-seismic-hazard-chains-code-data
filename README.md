# Code and data — Marine seismic hazard chains (NHESS)

Reproducibility materials for *Quantitative Modelling of Marine Seismic Hazard
Chains: Topology, Core Amplification Nodes, and Implications* (Shi, Hu, Zhang;
submitted to *Natural Hazards and Earth System Sciences*).

Manuscript (LaTeX) repository: <https://github.com/Chunyu-Hugh/NHESS-marine-seismic-hazard-chains>

**Interactive network viewer (Figs 1 & 2):** <http://chunyuhu.com/NHESS-marine-seismic-hazard-chains-code-data/>
— click any node to trace its couplings; Chinese / English toggle. Source in [`docs/`](docs/).

## One document reproduces everything

**[`reproduce.Rmd`](reproduce.Rmd)** regenerates **every table and figure** in the
paper directly from the raw cascade tables in [`data/`](data/). Knit it
(RStudio → *Knit*, or `rmarkdown::render("reproduce.Rmd")`) to produce
`reproduce.html`, which contains, each next to the exact code that makes it:

- Marine network — **39 nodes / 152 edges**
- **Table 2** node degrees · **Table 3** four-layer net flow (+42/+1/−16/−27) · **Table 4** 15-class marine (15/57) vs mainland (15/47)
- **Figure 1** marine in/out-degree network · **Figure 2** marine vs mainland (15 classes)
- **Sect. 3.4 / Figure 3** robustness — event-level cluster bootstrap + sampling sufficiency

A pre-knitted **`reproduce.html`** is included, so the full reproduction can be
read without running R.

### Requirements & how to run

R (≥ 4.x) with `docxtractr`, `igraph`, `dplyr` (plus `rmarkdown`, `knitr` to knit):

```r
install.packages(c("docxtractr", "igraph", "dplyr", "rmarkdown", "knitr"))
rmarkdown::render("reproduce.Rmd")        # from the repository root
```

Paths are **relative** and use `file.path()`, which is correct on Windows, Linux
and macOS alike — no per-OS configuration is needed (R also accepts `/` on
Windows). The OS the run used is printed in the document's Session information.

## Contents

| Path | What it holds |
|------|---------------|
| `reproduce.Rmd` / `reproduce.html` | the all-in-one reproduction (tables + figures + robustness) |
| `data/marine/tsunami_Dohmen_2025_EQ.docx` | marine cascade table (39-node network) |
| `data/mainland_china/earthquake_table_China_4_{E,N,NW,SW}.docx` | mainland China cascade tables (4 regions) |
| `data/Aggregation of the 39 marine hazard nodes into the 15 common.docx` | 39 → 15 class mapping (paper appendix) |
| `docs/` | interactive web viewer of the networks (GitHub Pages) |

*(Earlier standalone `.R` scripts were consolidated into `reproduce.Rmd`; they
remain recoverable in the git history.)*

## License and citation

*TODO(authors): add a license (e.g. CC-BY-4.0 for data, MIT for code) and, after
archiving this repository to Zenodo, insert the dataset DOI here and in the
manuscript's Code-and-data-availability statement.*
