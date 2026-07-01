# Data — documented seismic hazard-chain records

Input tables for the network analysis. Each is a Word table with the same three
columns; only the third is read by the code.

| Column | Meaning |
|--------|---------|
| `ID` | Event / source identifier for the cascade (e.g. `LTT_Honshu_2024`). |
| `Observed Cascading Effects` | Human-readable chain (e.g. `SEQ > Tsunami > Coastal Flooding > Death`). |
| `Proposed Encoding` | Machine-readable directed links of that chain, `trigger,target;trigger,target;…` (a `+` branch becomes two links). This is the only column the code parses. |

## Inventory

### `marine/`
- `tsunami_Dohmen_2025_EQ.docx` — **the canonical file behind the paper**
  (366 cascade rows). Source basis: Dohmen (2025) tsunami record plus the authors'
  encoding. **Verified (R 4.5.3 + igraph) to reproduce the paper exactly: 39 nodes
  / 152 edges**, every Table 2 degree (CF 14/13, TS 9/12, NF 12/5, …), and every
  Table 3 layer net flow (+42 / +1 / −16 / −27). This is the file the analysis
  script reads.
- `_extended_not_used_in_paper/tsunami_Dohmen_2025_EQ_CLS_SLS.docx` — a later,
  augmented version (396 rows; adds CLS/SLS links and a `TC` node). It yields
  41 nodes / 162 edges and does **not** match the paper; kept here only for
  reference. Do **not** deposit this one as the paper's dataset.

### `mainland_china/`
Continental earthquake hazard-chain records, split by region — the provenance of
the **continental comparison network** (paper Sect. 3.3):
- `earthquake_table_China_4_E.docx` — East China
- `earthquake_table_China_4_N.docx` — North China
- `earthquake_table_China_4_NW.docx` — Northwest China
- `earthquake_table_China_4_SW.docx` — Southwest China

(The code `rbind`s these four into a single continental network.)

## Edge-construction criterion (reusable wording for the paper Methods)
> A directed edge `i→j` is placed in the network if and only if the ordered peril
> pair (i, j) is documented in the `Proposed Encoding` of at least one recorded
> cascade chain; the union over all chains is then de-duplicated (self-loops
> retained). The analysed graph is thus the binarised one-step coupling matrix.

## ✅ Reproducibility — verified
Running the actual `multi_hazard_China_sy.R` graph build (R 4.5.3 + igraph):

- **Marine** (`marine/tsunami_Dohmen_2025_EQ.docx`): **39 nodes / 152 edges**, and
  every Table 2 degree and Table 3 layer net flow reproduce the paper exactly.
- **Mainland** (`mainland_china/` four files merged): **15 nodes / 52 edges**,
  matching the paper's 15-class network of 47 edges once self-loops and duplicate
  edges are removed (as the paper states).

**Mathematical cross-check** on the marine network (paper data): the metric
formulas match igraph; the catastrophe-dynamics series converges, ρ(**A**)=0.27<1,
and the truncation `M(τ=10)` equals the closed form `(I−A)⁻¹−I` to ~4×10⁻⁶, so the
`τ=10` used in the captions is well justified. One caveat: with the code's uniform
edge weight 0.1, the highest-out-degree node (SEQ, 29 out-edges) has row-sum 2.9>1,
so **A** is *not* substochastic — i.e. the 0.1 is a fixed coupling weight, not a
probability. The manuscript text around Eq. (8)–(9) should be reworded accordingly
(or the rows normalised).

> Note: an augmented version (`..._EQ_CLS_SLS.docx`, 41 nodes / 162 edges, with an
> extra `TC` node) exists but is **not** the paper's dataset and is intentionally
> excluded from this repository.
