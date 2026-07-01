# Code — marine & continental seismic hazard-chain networks

R code that builds the directed hazard-chain networks from the documented cascade
tables in `../data/`, computes the graph metrics reported in the paper, and draws
the network figures and the A/M interaction-matrix heatmaps.

## Files
| File | Role |
|------|------|
| `multi_hazard_China_sy.R` | Main script: reads a cascade table, builds the directed graph, computes degree/closeness/betweenness, draws the in/out-degree network figures, and the per-step A and M heatmaps. |
| `different_functions.R` | Helpers: `read_doc()` (extract table 1 from a `.docx`) and `mk_colors()` (bin transition values for the heatmaps). |

## Dependencies (R)
```r
install.packages(c("docxtractr","igraph","dplyr","ggplot2","reshape2",
                   "expm","gridExtra","Matrix","svglite"))
```

## How to run
1. Set the working directory to this `code/` folder.
2. In `multi_hazard_China_sy.R`, point `path1` at the table you want:
   - Marine network: `path1 <- "../data/marine/tsunami_Dohmen_2025_EQ.docx"`
     (this is the canonical file; it reproduces the paper's 39 nodes / 152 edges).
   - Continental network: read the four `../data/mainland_china/earthquake_table_China_4_*.docx`
     files (uncomment the `path2..path4` / `rbind(...)` lines).
3. `source("multi_hazard_China_sy.R")`.

### Outputs
- `fig_graph_indegree.svg`, `fig_graph_outdegree.svg` — the network figures
  (node colour & size = in/out-degree). These are the basis of the manuscript figures.
- `fig_graph_all_disasters.pdf` — both panels in one PDF.
- `EQ_AND_OTHER_CHINA_tau<τ>.pdf` — heatmaps of the one-step coupling matrix `A`
  and the cumulative interaction matrix `M = Σ_{i=1}^{τ} A^i`, for τ = 1…11.

## How the network is built (edge-construction criterion)
Each row of the input table carries a `Proposed Encoding` field listing the
directed links of one documented cascade, e.g. `SEQ,TS;TS,CFL;CFL,PP;`
(a `+` branch in the prose becomes two separate links). The script splits these
into `from,to` pairs over **all** documented chains, then **de-duplicates**: a
directed edge `i→j` is kept **iff the pair (i,j) is documented at least once**.
The graph is therefore the **binarised one-step coupling matrix A** (edge present
/ absent); node in/out-degree, density, centralization, path length and net flow
are all computed on this binary graph with `igraph`.

## Note on `τ` (relevant to the manuscript captions)
The `for (tau in 1:11)` loop only feeds the **A/M heatmaps**: it accumulates
`M = Σ_{i=1}^{τ} A^i` and shows it at each truncation depth (so "τ=10" = the
matrix summed to 10 interaction steps). **The degree-based metrics and the
network figures used in the manuscript do not depend on τ** — they come from the
binary graph above. See the manuscript `%TODO(authors)` on τ=10.

## ✅ Reproducibility — verified
The script reads `tsunami_Dohmen_2025_EQ.docx`. Run against
`../data/marine/tsunami_Dohmen_2025_EQ.docx` (R 4.5.3 + igraph) it reproduces the
paper **exactly: 39 nodes / 152 edges**, every Table 2 degree, and every Table 3
layer net flow. A mathematical cross-check confirms the metric formulas match
igraph and that the catastrophe-dynamics series converges (ρ(A)=0.27<1), with the
`τ=10` truncation matching the closed form `(I−A)⁻¹−I` to ~4×10⁻⁶. One caveat:
the uniform edge weight 0.1 makes the top node (SEQ) row-sum 2.9>1, so **A** is
not substochastic — reword Eq. (8)–(9) or normalise rows. The augmented
`..._EQ_CLS_SLS.docx` (41/162, extra `TC` node) under
`../data/marine/_extended_not_used_in_paper/` is *not* the paper's dataset.
