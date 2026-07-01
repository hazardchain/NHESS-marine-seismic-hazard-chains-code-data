# =====================================================================
# Robustness of the marine seismic hazard network
# Sampling unit = the source EVENT (tsunami). The 366 documented chains
# cluster under 58 events, so all resampling is done at the EVENT level
# (cluster / block bootstrap; Huang 2018; Efron & Tibshirani 1993).
#
#   (1) edge-frequency threshold sensitivity (edges seen in >=k events)
#   (2) event-level cluster bootstrap: 95% CIs + qualitative-conclusion stability
#       (feature-stability across resamples)
#   (3) sampling-sufficiency (edge-accumulation) curve over events
#       (cf. Casas et al. 2018 for network metrics)
#   (4) k-fold completeness: out-of-sample edge coverage (events split into folds)
#
# Outputs: results/robustness_summary.txt, results/fig_robustness.pdf/.png
# Run from code/:  Rscript sensitivity_crossvalidation.R
# =====================================================================
suppressMessages({library(docxtractr); library(igraph)})
set.seed(1)
DATA <- "../data/marine/tsunami_Dohmen_2025_EQ.docx"
OUT  <- "results"; dir.create(OUT, showWarnings = FALSE)

# ---- taxonomy layers (Table 1) --------------------------------------
TRIG <- c("SEQ","EQ","AS","VE")
NAT  <- c("SLS","CLS","LS","RO","LB","SA","LI","GU","SUB","MB","TS","CFL","FF","MED","WP","WD","DI","FA","SS")
ENG  <- c("PHF","SDF","SPF","SCF","NFA","FI","EX","CF","NF")
CON  <- c("MTD","BI","EC","PP","HD","IR","SU","CP")
layer_of <- function(v) ifelse(v%in%TRIG,"T",ifelse(v%in%NAT,"N",ifelse(v%in%ENG,"E",ifelse(v%in%CON,"C","?"))))

# ---- read table; forward-fill event ID; edges per chain -------------
x  <- data.frame(docx_extract_tbl(read_docx(DATA), tbl_number = 1))
id <- x$ID; cur <- NA; ev <- character(nrow(x))
for (i in seq_along(id)) { if (!is.na(id[i]) && trimws(id[i]) != "") cur <- id[i]; ev[i] <- cur }
chain_edges <- lapply(seq_len(nrow(x)), function(i){
  e <- c(); for (p in strsplit(x$Proposed.Encoding[i], ";")[[1]]) {
    pr <- strsplit(p, ",")[[1]]
    if (length(pr) >= 2 && !is.na(pr[1]) && !is.na(pr[2])) e <- c(e, paste(pr[1], pr[2], sep="|")) }
  unique(e)
})
events    <- unique(ev)
ev_edges  <- lapply(events, function(g) unique(unlist(chain_edges[ev == g])))  # edge set of each event
n_ev      <- length(events)
# edge -> number of distinct events it appears in
efreq <- table(unlist(ev_edges))
full  <- length(efreq)

build <- function(edges){ ee <- do.call(rbind, strsplit(edges, "\\|"))
  graph_from_data_frame(data.frame(from = ee[,1], to = ee[,2]), directed = TRUE) }
metrics <- function(g){
  din <- degree(g, mode="in"); dout <- degree(g, mode="out"); nm <- V(g)$name
  gv <- function(v,d) if (v %in% nm) d[v] else 0
  L  <- layer_of(nm); nf <- sapply(c("T","N","E","C"), function(k) sum(dout[L==k]) - sum(din[L==k]))
  list(N=vcount(g), E=ecount(g), SEQout=gv("SEQ",dout), EQout=gv("EQ",dout), CFin=gv("CF",din),
       TSin=gv("TS",din), NFin=gv("NF",din), NFout=gv("NF",dout),
       nfT=nf["T"], nfN=nf["N"], nfE=nf["E"], nfC=nf["C"],
       C1=setequal(nm[din==0], c("SEQ","EQ")),                            # only SEQ,EQ pure sources
       C2=(names(which.max(dout))=="SEQ" && names(which.max(din))=="CF"), # SEQ top-out, CF top-in
       C3=(gv("NF",din) > gv("NF",dout)),                                 # NF convergence node
       C4=(nf["T"]>nf["N"] && nf["N"]>nf["E"] && nf["E"]>nf["C"] && nf["T"]>0 && nf["C"]<0)) # monotone cascade
}

con <- file(file.path(OUT,"robustness_summary.txt"), open="wt")
w <- function(...) { cat(..., "\n"); cat(..., "\n", file=con) }
w(sprintf("Chains=%d  Events(tsunamis)=%d  Unique edges=%d\n", length(chain_edges), n_ev, full))

# ---- (1) edge-frequency (by event) threshold ------------------------
w("== (1) Edge-frequency threshold (edge seen in >=k distinct events) ==")
for (t in 1:3){ m <- metrics(build(names(efreq)[efreq>=t]))
  w(sprintf("k>=%d: N=%d E=%d | SEQout=%d CFin=%d NFin=%d NFout=%d | netT/N/E/C=%+d/%+d/%+d/%+d | C1=%s C2=%s C3=%s C4=%s",
            t,m$N,m$E,m$SEQout,m$CFin,m$NFin,m$NFout,m$nfT,m$nfN,m$nfE,m$nfC,m$C1,m$C2,m$C3,m$C4)) }

# ---- (2) event-level cluster bootstrap ------------------------------
w("\n== (2) Event-level cluster bootstrap (resample 58 events, B=1000) ==")
B <- 1000
res <- lapply(1:B, function(b){ idx <- sample(n_ev, n_ev, replace=TRUE)
  metrics(build(unique(unlist(ev_edges[idx])))) })
grab <- function(k) sapply(res, function(r) r[[k]])
ci   <- function(k){ v <- grab(k); sprintf("%.1f [%.0f, %.0f]", mean(v), quantile(v,.025), quantile(v,.975)) }
for (k in c("SEQout","EQout","CFin","TSin","NFin","NFout","nfT","nfN","nfE","nfC")) w(sprintf("  %-7s %s", k, ci(k)))
stab <- sapply(c("C1","C2","C3","C4"), function(k) 100*mean(grab(k)))
for (k in names(stab)) w(sprintf("  %s holds in %.1f%% of resamples", k, stab[k]))

# ---- (3) sampling-sufficiency (edge accumulation over events) -------
w("\n== (3) Sampling-sufficiency curve (add events) ==")
fr <- seq(0.05, 1, by=0.05)
acc <- sapply(fr, function(f){ k <- max(1, round(f*n_ev))
  mean(sapply(1:200, function(s) length(unique(unlist(ev_edges[sample(n_ev,k)]))))) })
for (i in c(5,10,15,20)) w(sprintf("  %3.0f%% events -> %.0f edges (%.0f%%)", 100*fr[i], acc[i], 100*acc[i]/full))

# ---- (4) k-fold completeness (out-of-sample edge coverage) ----------
w("\n== (4) 5-fold completeness: out-of-sample edge coverage (events in folds) ==")
fold <- sample(rep(1:5, length.out=n_ev))
recov <- function(sel) sapply(1:5, function(f){
  train <- unique(unlist(ev_edges[fold!=f])); test <- intersect(sel, unique(unlist(ev_edges[fold==f])))
  if (length(test)==0) NA else 100*sum(test %in% train)/length(test) })
covAll  <- mean(recov(names(efreq)),          na.rm=TRUE)
covCore <- mean(recov(names(efreq)[efreq>=2]), na.rm=TRUE)
covCor3 <- mean(recov(names(efreq)[efreq>=3]), na.rm=TRUE)
w(sprintf("  all edges:            %.1f%%", covAll))
w(sprintf("  core (>=2 events):    %.1f%%", covCore))
w(sprintf("  core (>=3 events):    %.1f%%", covCor3))
close(con)

# ---- figure ---------------------------------------------------------
mk <- function(dev){
  dev; par(mfrow=c(1,2), mar=c(4.2,4.2,2.4,1), cex.lab=1.05)
  plot(100*fr, 100*acc/full, type="o", pch=19, col="#d21d13", ylim=c(0,100),
       xlab="Source events (tsunamis) used (%)", ylab="Unique edges recovered (%)",
       main="(a) Sampling sufficiency & completeness")
  abline(h=covCore, lty=2, col="grey40")
  text(42, min(covCore-8, 88), sprintf("core edges (>=2 events):\n%.0f%% recovered out-of-sample", covCore),
       col="grey30", cex=0.8)
  barplot(stab, col="#F2A43A", ylim=c(0,100), ylab="Holds in resamples (%)",
          names.arg=c("C1 source","C2 hubs","C3 NF","C4 cascade"),
          main="(b) Conclusion stability (event bootstrap)")
  abline(h=95, lty=2, col="grey40")
  invisible(dev.off())
}
mk(pdf(file.path(OUT,"fig_robustness.pdf"), width=10, height=4.2))
mk(png(file.path(OUT,"fig_robustness.png"), width=1600, height=680, res=150))
cat("\nWrote results/robustness_summary.txt and fig_robustness.pdf/.png\n")
