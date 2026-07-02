## =====================================================================
## reproduce_paper_data.R
## Reproduces ONLY the numbers reported in the NHESS manuscript.
## All figure/plotting code has been removed on purpose; this script
## keeps only the data-producing logic.
##
##   Table 2  node in/out degrees      -> fine 39-node marine network
##   Table 3  four-layer net flow      -> fine 39-node marine network
##   Table 4  15-class comparison      -> marine (coarse) vs mainland
##
## The robustness numbers of Sect. 3.4 are produced separately by
##   sensitivity_crossvalidation.R
##
## ---------------------------------------------------------------------
## INPUT FILES (relative to code/, matching this repo's data/ layout):
##
##   data/marine/tsunami_Dohmen_2025_EQ.docx
##       (366 chains; SEQ/TS/CFL/... fine codes; reproduces Table 2/3/4-marine)
##   data/mainland_china/earthquake_table_China_4_E.docx
##   data/mainland_china/earthquake_table_China_4_N.docx
##   data/mainland_china/earthquake_table_China_4_SW.docx
##   data/mainland_china/earthquake_table_China_4_NW.docx
##   (all four China files must be DISTINCT; verify with md5sum -- a duplicate
##    slot silently drops the mainland edge count from 47 to 46)
##
## read_doc() must be provided by different_functions.R (a docx-table
## reader returning a data.frame with a $Proposed.Encoding column). The
## file currently named different_functions.R holds the robustness summary
## text, NOT this function; a docxtractr-based fallback is defined below.
## =====================================================================

suppressMessages({library(dplyr); library(igraph)})

## --- docx reader (fallback if read_doc() is unavailable) -------------
if (!exists("read_doc")) {
  if (!requireNamespace("docxtractr", quietly = TRUE))
    stop("Need read_doc() from different_functions.R, or install.packages('docxtractr').")
  read_doc <- function(path) {
    tb <- docxtractr::docx_extract_tbl(docxtractr::read_docx(path), tbl_number = 1)
    data.frame(tb, check.names = FALSE)
  }
}

DATA     <- "../data"                          # repo layout: data/marine, data/mainland_china
marine   <- function(f) file.path(DATA, "marine", f)
mainland <- function(f) file.path(DATA, "mainland_china", f)

## --- build a directed hazard-chain graph from one or more docx -------
## Parses $Proposed.Encoding: chains split on ';', pairs split on ','.
## trimws() avoids " TS" and "TS" being treated as different nodes.
build_graph <- function(paths) {
  x   <- do.call(rbind, lapply(paths, read_doc))
  enc <- x$Proposed.Encoding
  from <- character(0); to <- character(0)
  for (s in enc) for (seg in strsplit(s, ";")[[1]]) {
    pr <- trimws(strsplit(seg, ",")[[1]])
    if (length(pr) >= 2 && all(nzchar(pr[1:2]))) { from <- c(from, pr[1]); to <- c(to, pr[2]) }
  }
  final <- distinct(data.frame(from = from, to = to))   # unique directed edges (self-loops kept)
  graph_from_data_frame(final, directed = TRUE)
}

## =====================================================================
## FINE 39-node marine network  ->  Table 2 and Table 3
## =====================================================================
g_marine <- build_graph(marine("tsunami_Dohmen_2025_EQ.docx"))   # fine 39-code encoding, 366 chains

cat(sprintf("\n[marine fine]  nodes=%d  edges=%d\n", vcount(g_marine), ecount(g_marine)))

## Table 2 -- per-node in/out degree (self-loops excluded, as in the paper)
g2 <- simplify(g_marine, remove.multiple = TRUE, remove.loops = TRUE)
tab2 <- data.frame(node = V(g2)$name,
                   out  = degree(g2, mode = "out"),
                   `in` = degree(g2, mode = "in"),
                   net  = degree(g2, mode = "out") - degree(g2, mode = "in"),
                   check.names = FALSE)
cat("\n== Table 2: node degrees ==\n")
print(tab2[order(-tab2$out), ], row.names = FALSE)

## Table 3 -- four-layer net flow (Trigger / Natural / Engineering / Consequence)
layer_map <- c(
  SEQ="Trigger", EQ="Trigger", AS="Trigger", VE="Trigger",
  TS="Natural", CFL="Natural", GU="Natural", SLS="Natural", CLS="Natural",
  LI="Natural", MB="Natural", SUB="Natural", WP="Natural", LS="Natural",
  RO="Natural", WD="Natural", DI="Natural", LB="Natural", FA="Natural",
  FF="Natural", MED="Natural", SA="Natural",
  PHF="Engineering", SDF="Engineering", CF="Engineering", NF="Engineering",
  FI="Engineering", NFA="Engineering", EX="Engineering", SCF="Engineering", SPF="Engineering",
  BI="Consequence", IR="Consequence", PP="Consequence", MTD="Consequence",
  SU="Consequence", HD="Consequence", EC="Consequence", CP="Consequence")
L4 <- c("Trigger", "Natural", "Engineering", "Consequence")

layer_flow <- function(g, lm, levels) {
  lay <- factor(lm[V(g)$name], levels = levels)
  o <- tapply(degree(g, mode = "out"), lay, sum); o[is.na(o)] <- 0
  i <- tapply(degree(g, mode = "in"),  lay, sum); i[is.na(i)] <- 0
  data.frame(layer = levels, nodes = as.integer(table(lay)[levels]),
             out = as.integer(o[levels]), `in` = as.integer(i[levels]),
             net = as.integer(o[levels] - i[levels]),
             net_pct = round((o[levels] - i[levels]) / ecount(g) * 100, 1),
             check.names = FALSE)
}
cat("\n== Table 3: four-layer net flow (fine network) ==\n")
print(layer_flow(g_marine, layer_map, L4), row.names = FALSE)   # expected net: +42 +1 -16 -27

## =====================================================================
## 15-class comparison  ->  Table 4  (marine coarse vs mainland)
## =====================================================================
## 39 -> 15 aggregation, aligned to the mainland vocabulary.
coarse_map <- c(
  SEQ="EQ or AS", AS="EQ or AS", EQ="EQ or AS", VE="EQ or AS",
  GU="Mass Slide", SLS="Mass Slide", CLS="Mass Slide", LI="Mass Slide",
  MB="Mass Slide", SUB="Mass Slide", LS="Mass Slide", RO="Mass Slide",
  LB="Mass Slide", SA="Mass Slide",
  PHF="Critical Infra. Fail.", SDF="Critical Infra. Fail.", CF="Critical Infra. Fail.",
  NFA="Critical Infra. Fail.", SPF="Critical Infra. Fail.",
  BI="Business Interr.", MTD="Business Interr.",
  WP="Water Dil.", WD="Water Dil.",
  NF="Critical NET Fail.", SCF="Critical NET Fail.",
  FA="Flora an Fauna", MED="Flora an Fauna",
  FI="Fire", EX="Fire",
  IR="Healthcare Deg.", HD="Healthcare Deg.",
  TS="Flood", CFL="Flood", FF="Flood",
  EC="Econ. Crisis", DI="Disease", PP="Phy.-Psych. Trauma",
  CP="Culture Perils", SU="Social Unrest")

layer15 <- c(
  "EQ or AS"="Trigger",
  "Mass Slide"="Natural", "Water Dil."="Natural", "Flora an Fauna"="Natural",
  "Flood"="Natural", "Disease"="Natural",
  "Critical Infra. Fail."="Engineering", "Critical NET Fail."="Engineering", "Fire"="Engineering",
  "Business Interr."="Consequence", "Healthcare Deg."="Consequence",
  "Econ. Crisis"="Consequence", "Phy.-Psych. Trauma"="Consequence",
  "Culture Perils"="Consequence", "Social Unrest"="Consequence")

## contract the fine marine graph to 15 classes; mainland is already 15-class
miss <- setdiff(V(g_marine)$name, names(coarse_map))
if (length(miss)) stop("unmapped marine nodes: ", paste(miss, collapse = ", "))
grp        <- factor(coarse_map[V(g_marine)$name])
g_marine_c <- contract(g_marine, as.integer(grp), vertex.attr.comb = "first")
V(g_marine_c)$name <- levels(grp)
g_marine_c   <- simplify(g_marine_c, remove.multiple = TRUE, remove.loops = TRUE)

g_mainland   <- build_graph(c(mainland("earthquake_table_China_4_E.docx"),
                              mainland("earthquake_table_China_4_N.docx"),
                              mainland("earthquake_table_China_4_SW.docx"),
                              mainland("earthquake_table_China_4_NW.docx")))
g_mainland_c <- simplify(g_mainland, remove.multiple = TRUE, remove.loops = TRUE)

## Table 4 -- comparable topology metrics
topo <- function(g, name) data.frame(
  network = name, nodes = vcount(g), edges = ecount(g),
  density   = round(edge_density(g, loops = FALSE), 3),
  max_out   = max(degree(g, mode = "out")),
  Cout      = round(centr_degree(g, mode = "out", loops = FALSE)$centralization, 3),
  Cin       = round(centr_degree(g, mode = "in",  loops = FALSE)$centralization, 3),
  sources   = sum(degree(g, mode = "in")  == 0),
  sinks     = sum(degree(g, mode = "out") == 0),
  mean_path = round(mean_distance(g, directed = TRUE), 2))

cat("\n== Table 4: 15-class comparison ==\n")
print(rbind(topo(g_marine_c, "Marine(coarse)"), topo(g_mainland_c, "Mainland")), row.names = FALSE)
cat("\n-- marine layer net flow --\n"); print(layer_flow(g_marine_c,   layer15, L4), row.names = FALSE)
cat("\n-- mainland layer net flow --\n"); print(layer_flow(g_mainland_c, layer15, L4), row.names = FALSE)
