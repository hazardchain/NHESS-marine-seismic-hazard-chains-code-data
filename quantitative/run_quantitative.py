"""Quantitative extension of the marine seismic hazard-chain network.

Turns the binary 39-node network into a probabilistic cascade model using only
the event-level information already present in data/:

  1. edge triggering probabilities p_ij estimated from the 58 source events
     (how often j is documented downstream of i, given that i occurred);
  2. the weighted transition matrix A, its spectral radius and the cumulative
     interaction matrix M = (I-A)^-1 - I of the paper's Eqs. (10)-(12);
  3. an independent-cascade Monte Carlo from the seismic seed: per-node
     activation probability, cascade-size distribution, consequence reach;
  4. node state s (activation probability / downstream consequence potential),
     the degree-state product d*s and corr(d,s) of Zheng et al. (2026);
  5. exact (exhaustive) search of the k-node removal set that lowers the
     expected consequence reach the most, k = 1..3;
  6. event-level bootstrap intervals for the key outputs;
  7. probability-weighted check of the three held-out events;
  8. the same p_ij / cascade / corr(d,s) pipeline for marine-15 vs mainland-15.

Run from the repository root:
    python -I quantitative/run_quantitative.py
Outputs go to quantitative/outputs/ and quantitative/figures/.
"""
import csv
import glob
import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from chains_io import load_chains, load_many, event_edge_sets, event_node_sets  # noqa: E402
from taxonomy import LAYER_FINE, COARSE_MAP, LAYER15, ABBR15, FULL_NAME  # noqa: E402

OUT = os.path.join(HERE, "outputs")
FIG = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

MARINE = os.path.join(ROOT, "data", "marine", "tsunami_Dohmen_2025_EQ.docx")
EXTRA = os.path.join(ROOT, "data", "marine", "tsunami_extra_events_2004_1960_2018_EQ.docx")
MAINLAND = sorted(glob.glob(os.path.join(ROOT, "data", "mainland_china", "earthquake_table_China_4_*.docx")))

SEED = 20260730
N_SIMS = 20000       # activation probabilities / cascade sizes
N_SIMS_OPT = 4000    # exhaustive removal search (common random numbers)
N_BOOT = 300         # event-level bootstrap replicates
N_SIMS_BOOT = 3000


# --------------------------------------------------------------------------- #
# 1. Estimation of edge triggering probabilities
# --------------------------------------------------------------------------- #
def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def estimate(ev_edges, ev_nodes):
    """p_ij = (#events with i->j documented) / (#events in which i occurs)."""
    n_node, n_edge = {}, {}
    for ev in ev_edges:
        for v in ev_nodes[ev]:
            n_node[v] = n_node.get(v, 0) + 1
        for e in ev_edges[ev]:
            n_edge[e] = n_edge.get(e, 0) + 1
    rows = {}
    for (i, j), k in n_edge.items():
        n = n_node[i]
        lo, hi = wilson(k, n)
        rows[(i, j)] = dict(n_i=n, n_ij=k, p_mle=k / n,
                            p_shrunk=(k + 0.5) / (n + 1.0), ci_lo=lo, ci_hi=hi)
    return rows, n_node


def build_A(nodes, p_rows, key="p_mle"):
    idx = {v: k for k, v in enumerate(nodes)}
    A = np.zeros((len(nodes), len(nodes)))
    for (i, j), r in p_rows.items():
        if i != j:                       # self-loop PP->PP carries no propagation
            A[idx[i], idx[j]] = r[key]
    return A, idx


# --------------------------------------------------------------------------- #
# 3. Independent-cascade Monte Carlo (live-edge sampling, vectorised)
# --------------------------------------------------------------------------- #
def sample_live(A, n_sims, rng):
    return (rng.random((n_sims, A.shape[0], A.shape[0])) < A[None, :, :]).astype(np.float32)


def propagate(live, seed_idx, removed=None):
    """reached[s, v] = 1 if v is activated in simulation s (seed always active)."""
    S, n, _ = live.shape
    reached = np.zeros((S, n), dtype=np.float32)
    reached[:, seed_idx] = 1.0
    keep = np.ones(n, dtype=np.float32)
    if removed is not None:
        keep[list(removed)] = 0.0
        keep[seed_idx] = 1.0
    for _ in range(n):
        new = (np.matmul(reached[:, None, :], live)[:, 0, :] > 0).astype(np.float32) * keep
        nxt = np.maximum(reached, new)
        if np.array_equal(nxt, reached):
            break
        reached = nxt
    return reached


def cascade_stats(reached, layer_mask_cons):
    size = reached.sum(axis=1) - 1                      # activated nodes besides the seed
    cons = (reached * layer_mask_cons[None, :]).sum(axis=1)
    return dict(p_active=reached.mean(axis=0), size=size, cons=cons)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return pearson(rx, ry)


def aggregate_events(events, mapping):
    """Apply a node->class mapping to every chain of every event."""
    out = {}
    for ev, chains in events.items():
        out[ev] = [[(mapping.get(a, a), mapping.get(b, b)) for a, b in ch] for ch in chains]
    return out


# --------------------------------------------------------------------------- #
# full pipeline for one network (used for marine-39, marine-15, mainland-15)
# --------------------------------------------------------------------------- #
def analyse(events, layer_of, seed_node, tag, rng, n_sims=N_SIMS, do_M=True):
    ev_edges = event_edge_sets(events)
    ev_nodes = event_node_sets(events)
    p_rows, n_node = estimate(ev_edges, ev_nodes)
    nodes = sorted({v for e in p_rows for v in e})
    A, idx = build_A(nodes, p_rows)
    n = len(nodes)
    cons_mask = np.array([1.0 if layer_of[v] == "Consequence" else 0.0 for v in nodes], dtype=np.float32)

    res = dict(tag=tag, nodes=nodes, n_events=len(ev_edges), n_edges=int((A > 0).sum()),
               p_rows=p_rows, n_node=n_node, A=A, idx=idx)

    # spectral radius and cumulative interaction matrix
    rho = float(np.max(np.abs(np.linalg.eigvals(A))))
    res["rho"] = rho
    if do_M and rho < 1:
        M = np.linalg.inv(np.eye(n) - A) - np.eye(n)
        res["M"] = M

    # Monte Carlo from the seismic seed
    live = sample_live(A, n_sims, rng)
    reached = propagate(live, idx[seed_node])
    st = cascade_stats(reached, cons_mask)
    res["p_active"] = st["p_active"]
    res["size"] = st["size"]
    res["cons"] = st["cons"]
    del live

    # downstream consequence potential of every node (seeded at that node)
    live2 = sample_live(A, 5000, rng)
    s_cons = np.zeros(n)
    for v in nodes:
        r = propagate(live2, idx[v])
        s_cons[idx[v]] = (r * cons_mask[None, :]).sum(axis=1).mean()
    res["s_cons"] = s_cons
    del live2

    # degrees
    res["d_out_w"] = A.sum(axis=1)
    res["d_in_w"] = A.sum(axis=0)
    res["d_out"] = (A > 0).sum(axis=1).astype(float)
    res["d_in"] = (A > 0).sum(axis=0).astype(float)
    return res


def corr_table(res, seed_node):
    nodes, idx = res["nodes"], res["idx"]
    s = res["p_active"]; sc = res["s_cons"]
    keep = np.array([v != seed_node for v in nodes])   # exclude the seed (s=1 by construction)
    rows = []
    for dname in ["d_out_w", "d_in_w", "d_tot_w", "d_out", "d_in", "d_tot"]:
        if dname == "d_tot_w":
            d = res["d_out_w"] + res["d_in_w"]
        elif dname == "d_tot":
            d = res["d_out"] + res["d_in"]
        else:
            d = res[dname]
        rows.append((dname, pearson(d, s), pearson(d[keep], s[keep]), pearson(d, sc), pearson(d[keep], sc[keep])))
    return rows


# --------------------------------------------------------------------------- #
def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    summary = {}

    # ---------------- marine 39-node network ----------------
    ev_marine = load_chains(MARINE)
    R = analyse(ev_marine, LAYER_FINE, "SEQ", "marine39", rng)
    nodes, idx, A = R["nodes"], R["idx"], R["A"]
    n = len(nodes)
    print(f"[marine39] events={R['n_events']} nodes={n} edges(no self-loop)={R['n_edges']} rho(A)={R['rho']:.3f}")
    summary["marine39"] = dict(n_events=R["n_events"], n_nodes=n, n_edges=R["n_edges"], rho=R["rho"],
                               row_sum_max=float(A.sum(axis=1).max()),
                               row_sum_max_node=nodes[int(A.sum(axis=1).argmax())])

    # 1. edge probabilities
    rows = []
    for (i, j), r in sorted(R["p_rows"].items(), key=lambda kv: -kv[1]["p_mle"]):
        rows.append([i, j, LAYER_FINE[i], LAYER_FINE[j], r["n_i"], r["n_ij"], f"{r['p_mle']:.3f}",
                     f"{r['p_shrunk']:.3f}", f"{r['ci_lo']:.3f}", f"{r['ci_hi']:.3f}"])
    write_csv(os.path.join(OUT, "marine39_edge_probabilities.csv"),
              ["from", "to", "layer_from", "layer_to", "n_events_from", "n_events_edge",
               "p_mle", "p_jeffreys", "wilson95_lo", "wilson95_hi"], rows)
    pm = np.array([r["p_mle"] for r in R["p_rows"].values()])
    summary["marine39"]["p_quantiles"] = {q: float(np.quantile(pm, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)}
    summary["marine39"]["edges_single_event"] = int(sum(1 for r in R["p_rows"].values() if r["n_ij"] == 1))

    # 2. M matrix: cumulative influence from the two seismic sources
    if "M" in R:
        M = R["M"]
        np.savetxt(os.path.join(OUT, "marine39_A_weighted.csv"), A, delimiter=",", fmt="%.4f", header=",".join(nodes), comments="")
        np.savetxt(os.path.join(OUT, "marine39_M_cumulative.csv"), M, delimiter=",", fmt="%.4f", header=",".join(nodes), comments="")
        mrows = []
        for src in ["SEQ", "EQ"]:
            for v in nodes:
                if v != src:
                    mrows.append([src, v, LAYER_FINE[v], f"{A[idx[src], idx[v]]:.4f}", f"{M[idx[src], idx[v]]:.4f}"])
        mrows.sort(key=lambda r: (r[0], -float(r[4])))
        write_csv(os.path.join(OUT, "marine39_M_from_sources.csv"),
                  ["source", "target", "layer", "direct_A", "cumulative_M"], mrows)
        # uniform-weight comparison used in the first manuscript (0.1 on every edge)
        Au = (A > 0) * 0.1
        summary["marine39"]["rho_uniform_0.1"] = float(np.max(np.abs(np.linalg.eigvals(Au))))

    # 3. activation probabilities and cascade size from SEQ
    pa = R["p_active"]
    arows = sorted([[v, FULL_NAME[v], LAYER_FINE[v], f"{pa[idx[v]]:.4f}"] for v in nodes], key=lambda r: -float(r[3]))
    write_csv(os.path.join(OUT, "marine39_activation_from_SEQ.csv"), ["node", "name", "layer", "p_activated_given_SEQ"], arows)
    size, cons = R["size"], R["cons"]
    summary["marine39"]["cascade_from_SEQ"] = dict(
        mean_size=float(size.mean()), median_size=float(np.median(size)),
        size_q=[float(np.quantile(size, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)],
        mean_consequence_nodes=float(cons.mean()),
        p_any_consequence=float((cons > 0).mean()),
        p_PP=float(pa[idx["PP"]]), p_BI=float(pa[idx["BI"]]), p_EC=float(pa[idx["EC"]]),
        p_CF=float(pa[idx["CF"]]), p_NF=float(pa[idx["NF"]]), p_TS=float(pa[idx["TS"]]),
        p_NFA=float(pa[idx["NFA"]]))
    # size histogram
    hist = np.bincount(size.astype(int), minlength=n)
    write_csv(os.path.join(OUT, "marine39_cascade_size_distribution.csv"), ["activated_nodes_besides_seed", "probability"],
              [[k, f"{hist[k] / size.size:.5f}"] for k in range(len(hist))])

    # 4. node state, d*s, corr(d,s)
    s_act, s_cons = R["p_active"], R["s_cons"]
    dsrows = []
    for v in nodes:
        k = idx[v]
        dsrows.append([v, LAYER_FINE[v], int(R["d_out"][k]), int(R["d_in"][k]), f"{R['d_out_w'][k]:.3f}", f"{R['d_in_w'][k]:.3f}",
                       f"{s_act[k]:.4f}", f"{s_cons[k]:.3f}",
                       f"{R['d_out_w'][k] * s_act[k]:.4f}", f"{R['d_in_w'][k] * s_act[k]:.4f}",
                       f"{R['d_out_w'][k] * s_cons[k]:.4f}", f"{(R['d_out_w'][k] + R['d_in_w'][k]) * s_act[k]:.4f}"])
    dsrows.sort(key=lambda r: -float(r[11]))
    write_csv(os.path.join(OUT, "marine39_node_state_ds.csv"),
              ["node", "layer", "d_out", "d_in", "d_out_w", "d_in_w", "s_act(P|SEQ)", "s_cons(E#cons from node)",
               "dout_w*s_act", "din_w*s_act", "dout_w*s_cons", "dtot_w*s_act"], dsrows)
    crows = corr_table(R, "SEQ")
    write_csv(os.path.join(OUT, "marine39_corr_d_s.csv"),
              ["degree", "corr(d,s_act)", "corr(d,s_act) excl. seed", "corr(d,s_cons)", "corr(d,s_cons) excl. seed"],
              [[r[0]] + [f"{x:.3f}" for x in r[1:]] for r in crows])
    summary["marine39"]["corr_d_s"] = {r[0]: dict(all=r[1], excl_seed=r[2], s_cons=r[3], s_cons_excl=r[4]) for r in crows}

    # 5. exhaustive removal search (seed SEQ; sources cannot be removed)
    cons_mask = np.array([1.0 if LAYER_FINE[v] == "Consequence" else 0.0 for v in nodes], dtype=np.float32)
    live = sample_live(A, N_SIMS_OPT, rng)
    base = propagate(live, idx["SEQ"])
    base_cons = float((base * cons_mask[None, :]).sum(axis=1).mean())
    base_size = float(base.sum(axis=1).mean() - 1)
    # Intervention candidates are the process-layer nodes only. Triggers cannot
    # be prevented, and removing a consequence node would delete the very
    # outcome being counted (the automatic-reduction bias raised by RC1).
    candidates = [v for v in nodes if LAYER_FINE[v] in ("Natural", "Engineering")]
    opt = {}
    all_rows = []
    for k in (1, 2, 3):
        best = []
        for combo in itertools.combinations(candidates, k):
            r = propagate(live, idx["SEQ"], removed=[idx[v] for v in combo])
            ec = float((r * cons_mask[None, :]).sum(axis=1).mean())
            es = float(r.sum(axis=1).mean() - 1)
            pany = float(((r * cons_mask[None, :]).sum(axis=1) > 0).mean())
            best.append((ec, es, pany, combo))
        best.sort(key=lambda t: t[0])
        opt[k] = best[:15]
        for ec, es, pany, combo in best[:15]:
            all_rows.append([k, "+".join(combo), f"{ec:.3f}", f"{100 * (base_cons - ec) / base_cons:.1f}",
                             f"{es:.3f}", f"{100 * (base_size - es) / base_size:.1f}", f"{pany:.3f}"])
        print(f"[opt k={k}] best {best[0][3]} E[#cons] {base_cons:.3f}->{best[0][0]:.3f}  ({len(best)} sets)")
    # reference sets from the manuscript
    for combo in [("CF", "FI", "TS"), ("CF", "NF", "TS"), ("CF",), ("TS",), ("NF",), ("FI",)]:
        r = propagate(live, idx["SEQ"], removed=[idx[v] for v in combo])
        ec = float((r * cons_mask[None, :]).sum(axis=1).mean()); es = float(r.sum(axis=1).mean() - 1)
        pany = float(((r * cons_mask[None, :]).sum(axis=1) > 0).mean())
        all_rows.append(["ref", "+".join(combo), f"{ec:.3f}", f"{100 * (base_cons - ec) / base_cons:.1f}",
                         f"{es:.3f}", f"{100 * (base_size - es) / base_size:.1f}", f"{pany:.3f}"])
    write_csv(os.path.join(OUT, "marine39_removal_search.csv"),
              ["k", "removed_set", "E[#consequence nodes]", "reduction_%", "E[cascade size]", "size_reduction_%",
               "P[any consequence]"], all_rows)
    summary["marine39"]["removal"] = dict(
        baseline_E_cons=base_cons, baseline_E_size=base_size,
        best={k: dict(set=list(opt[k][0][3]), E_cons=opt[k][0][0], reduction_pct=100 * (base_cons - opt[k][0][0]) / base_cons)
              for k in opt})
    del live

    # 6. event-level bootstrap of key outputs
    events_list = list(ev_marine.keys())
    boot = dict(p_PP=[], p_CF=[], p_NF=[], p_TS=[], p_NFA=[], E_cons=[], rho=[], corr_out_w=[], corr_in_w=[],
                top1_removal=[], ds_top=[])
    for b in range(N_BOOT):
        pick = rng.choice(len(events_list), size=len(events_list), replace=True)
        ev_b = {}
        for t, ii in enumerate(pick):
            ev_b[f"{events_list[ii]}#{t}"] = ev_marine[events_list[ii]]
        pe, pn = estimate(event_edge_sets(ev_b), event_node_sets(ev_b))
        nodes_b = sorted({v for e in pe for v in e})
        if "SEQ" not in nodes_b:
            continue
        Ab, ib = build_A(nodes_b, pe)
        boot["rho"].append(float(np.max(np.abs(np.linalg.eigvals(Ab)))))
        cm = np.array([1.0 if LAYER_FINE[v] == "Consequence" else 0.0 for v in nodes_b], dtype=np.float32)
        lv = sample_live(Ab, N_SIMS_BOOT, rng)
        rb = propagate(lv, ib["SEQ"])
        pab = rb.mean(axis=0)
        for key in ("PP", "CF", "NF", "TS", "NFA"):
            boot["p_" + key].append(float(pab[ib[key]]) if key in ib else 0.0)
        boot["E_cons"].append(float((rb * cm[None, :]).sum(axis=1).mean()))
        keep = np.array([v != "SEQ" for v in nodes_b])
        boot["corr_out_w"].append(pearson(Ab.sum(axis=1)[keep], pab[keep]))
        boot["corr_in_w"].append(pearson(Ab.sum(axis=0)[keep], pab[keep]))
        dtot = (Ab.sum(axis=1) + Ab.sum(axis=0)) * pab
        dtot[ib["SEQ"]] = -1
        boot["ds_top"].append(nodes_b[int(dtot.argmax())])
        # single-node removal winner
        bestv, bestval = None, None
        for v in nodes_b:
            if LAYER_FINE[v] not in ("Natural", "Engineering"):
                continue
            rr = propagate(lv, ib["SEQ"], removed=[ib[v]])
            val = float((rr * cm[None, :]).sum(axis=1).mean())
            if bestval is None or val < bestval:
                bestval, bestv = val, v
        boot["top1_removal"].append(bestv)
        del lv
    def ci(x):
        x = np.asarray(x, float)
        return [float(np.quantile(x, 0.025)), float(np.quantile(x, 0.5)), float(np.quantile(x, 0.975))]
    bsum = {k: ci(v) for k, v in boot.items() if k not in ("top1_removal", "ds_top")}
    for key in ("top1_removal", "ds_top"):
        vals, cnt = np.unique(boot[key], return_counts=True)
        bsum[key] = {str(v): int(c) for v, c in sorted(zip(vals, cnt), key=lambda t: -t[1])}
    bsum["n_replicates"] = len(boot["rho"])
    summary["marine39"]["bootstrap"] = bsum
    print("[bootstrap]", json.dumps({k: v for k, v in bsum.items() if k in ("p_PP", "E_cons", "top1_removal", "ds_top")}))

    # 7. held-out events with probabilities
    ev_extra = load_chains(EXTRA)
    ex_edges = event_edge_sets(ev_extra)
    hrows = []
    allp = np.array([r["p_mle"] for r in R["p_rows"].values()])
    for ev, es in ex_edges.items():
        for (i, j) in sorted(es):
            r = R["p_rows"].get((i, j))
            if r is None:
                hrows.append([ev, i, j, "novel", "", "", ""])
            else:
                pct = float((allp <= r["p_mle"]).mean())
                hrows.append([ev, i, j, "in network", r["n_ij"], f"{r['p_mle']:.3f}", f"{pct:.2f}"])
    write_csv(os.path.join(OUT, "heldout_events_edge_probabilities.csv"),
              ["event", "from", "to", "status", "n_events_edge", "p_mle", "percentile_among_network_edges"], hrows)
    covered_p = [float(r[5]) for r in hrows if r[3] == "in network"]
    union_links = set().union(*ex_edges.values())
    summary["heldout"] = dict(
        links_union=len(union_links),
        covered=sum(1 for e in union_links if e in R["p_rows"]),
        mean_p_covered=float(np.mean(covered_p)), median_p_covered=float(np.median(covered_p)),
        mean_p_all_edges=float(allp.mean()), median_p_all_edges=float(np.median(allp)),
        frac_covered_above_network_median=float(np.mean([p > np.median(allp) for p in covered_p])),
        frac_covered_single_event_edges=float(np.mean([r[4] == 1 for r in hrows if r[3] == "in network"])))

    # 8. marine-15 vs mainland-15
    ev_m15 = aggregate_events(ev_marine, COARSE_MAP)
    ev_main = load_many(MAINLAND)
    R15 = analyse(ev_m15, LAYER15, "EQ or AS", "marine15", rng, n_sims=10000)
    RM = analyse(ev_main, LAYER15, "EQ or AS", "mainland15", rng, n_sims=10000)
    comp = []
    for RR in (R15, RM):
        cr = {r[0]: r for r in corr_table(RR, "EQ or AS")}
        comp.append(dict(tag=RR["tag"], n_events=RR["n_events"], n_nodes=len(RR["nodes"]), n_edges=RR["n_edges"], rho=RR["rho"],
                         mean_size=float(RR["size"].mean()), mean_cons=float(RR["cons"].mean()),
                         corr_dout_w_s=cr["d_out_w"][2], corr_din_w_s=cr["d_in_w"][2], corr_dtot_w_s=cr["d_tot_w"][2],
                         corr_dtot_s=cr["d_tot"][2]))
        prow = []
        for v in RR["nodes"]:
            k = RR["idx"][v]
            prow.append([ABBR15.get(v, v), v, LAYER15[v], int(RR["d_out"][k]), int(RR["d_in"][k]),
                         f"{RR['d_out_w'][k]:.3f}", f"{RR['d_in_w'][k]:.3f}", f"{RR['p_active'][k]:.4f}", f"{RR['s_cons'][k]:.3f}",
                         f"{(RR['d_out_w'][k] + RR['d_in_w'][k]) * RR['p_active'][k]:.4f}"])
        write_csv(os.path.join(OUT, f"{RR['tag']}_node_state_ds.csv"),
                  ["abbr", "class", "layer", "d_out", "d_in", "d_out_w", "d_in_w", "s_act", "s_cons", "dtot_w*s_act"], prow)
        erow = [[i, j, r["n_i"], r["n_ij"], f"{r['p_mle']:.3f}", f"{r['ci_lo']:.3f}", f"{r['ci_hi']:.3f}"]
                for (i, j), r in sorted(RR["p_rows"].items(), key=lambda kv: -kv[1]["p_mle"])]
        write_csv(os.path.join(OUT, f"{RR['tag']}_edge_probabilities.csv"),
                  ["from", "to", "n_events_from", "n_events_edge", "p_mle", "wilson95_lo", "wilson95_hi"], erow)
    summary["comparison15"] = comp
    # bootstrap the corr(d,s) difference at 15 classes
    def boot_corr(events, nb):
        keys = list(events.keys()); out = []
        for b in range(nb):
            pick = rng.choice(len(keys), size=len(keys), replace=True)
            ev_b = {f"{keys[i]}#{t}": events[keys[i]] for t, i in enumerate(pick)}
            pe, pn = estimate(event_edge_sets(ev_b), event_node_sets(ev_b))
            nb_nodes = sorted({v for e in pe for v in e})
            Ab, ib = build_A(nb_nodes, pe)
            lv = sample_live(Ab, 2000, rng)
            pab = propagate(lv, ib["EQ or AS"]).mean(axis=0)
            keep = np.array([v != "EQ or AS" for v in nb_nodes])
            out.append(pearson((Ab.sum(axis=1) + Ab.sum(axis=0))[keep], pab[keep]))
        return out
    bm = boot_corr(ev_m15, 200); bl = boot_corr(ev_main, 200)
    diff = np.array(bm) - np.array(bl)
    summary["comparison15_bootstrap"] = dict(marine_corr_ci=ci(bm), mainland_corr_ci=ci(bl), diff_ci=ci(diff),
                                             p_marine_lower=float((diff < 0).mean()))
    print("[15-class] corr(d_tot_w, s) marine", f"{comp[0]['corr_dtot_w_s']:.3f}", "mainland", f"{comp[1]['corr_dtot_w_s']:.3f}",
          "diff CI", [round(x, 3) for x in summary["comparison15_bootstrap"]["diff_ci"]])

    # ---------------- figures ----------------
    try:
        make_figures(R, R15, RM, summary)
    except Exception as e:  # figures are optional
        print("figure generation failed:", repr(e))

    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"done in {time.time() - t0:.0f}s")


# --------------------------------------------------------------------------- #
def make_figures(R, R15, RM, summary):
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nodes, idx, A = R["nodes"], R["idx"], R["A"]
    lay_col = {"Trigger": "#e4572e", "Natural": "#2a9d8f", "Engineering": "#e9a220", "Consequence": "#4a6fa5"}

    # Q1: edge probability distribution + top edges with Wilson intervals
    items = sorted(R["p_rows"].items(), key=lambda kv: -kv[1]["p_mle"])
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    pm = [r["p_mle"] for _, r in items]
    ax[0].hist(pm, bins=20, color="#4a6fa5", edgecolor="white")
    ax[0].set_xlabel("estimated triggering probability $p_{ij}$"); ax[0].set_ylabel("edges")
    ax[0].set_title("(a) distribution over the 151 documented couplings")
    top = items[:25][::-1]
    y = np.arange(len(top))
    ax[1].barh(y, [r["p_mle"] for _, r in top], color="#d21d13", alpha=.8)
    ax[1].errorbar([r["p_mle"] for _, r in top], y,
                   xerr=[[r["p_mle"] - r["ci_lo"] for _, r in top], [r["ci_hi"] - r["p_mle"] for _, r in top]],
                   fmt="none", ecolor="grey", capsize=2)
    ax[1].set_yticks(y); ax[1].set_yticklabels([f"{i}→{j}  (n={r['n_ij']}/{r['n_i']})" for (i, j), r in top], fontsize=8)
    ax[1].set_xlabel("$p_{ij}$ with Wilson 95% interval"); ax[1].set_title("(b) 25 most probable couplings")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "figQ1_edge_probabilities.png"), dpi=200); plt.close(fig)

    # Q2: activation probability from SEQ, by layer
    pa = R["p_active"]
    order = sorted(nodes, key=lambda v: (["Trigger", "Natural", "Engineering", "Consequence"].index(LAYER_FINE[v]), -pa[idx[v]]))
    fig, ax = plt.subplots(figsize=(13, 4.5))
    ax.bar(range(len(order)), [pa[idx[v]] for v in order], color=[lay_col[LAYER_FINE[v]] for v in order])
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=8)
    ax.set_ylabel("P(activated | submarine earthquake)")
    ax.set_title("independent-cascade Monte Carlo from SEQ (20 000 runs, event-frequency $p_{ij}$)")
    for L, c in lay_col.items():
        ax.bar([0], [0], color=c, label=L)
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "figQ2_activation_from_SEQ.png"), dpi=200); plt.close(fig)

    # Q3: cascade-size distribution
    size = R["size"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].hist(size, bins=np.arange(-0.5, size.max() + 1.5), color="#2a9d8f", edgecolor="white", density=True)
    ax[0].set_xlabel("activated nodes besides SEQ"); ax[0].set_ylabel("probability"); ax[0].set_title("(a) cascade size")
    cons = R["cons"]
    ax[1].hist(cons, bins=np.arange(-0.5, cons.max() + 1.5), color="#4a6fa5", edgecolor="white", density=True)
    ax[1].set_xlabel("consequence-layer nodes reached"); ax[1].set_title("(b) consequence reach")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "figQ3_cascade_size.png"), dpi=200); plt.close(fig)

    # Q4: d vs s scatter (marine 39) and 15-class comparison
    d = R["d_out_w"] + R["d_in_w"]; s = pa
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    for v in nodes:
        k = idx[v]
        ax[0].scatter(d[k], s[k], color=lay_col[LAYER_FINE[v]], s=30 + 400 * d[k] * s[k] / max(1e-9, (d * s).max()), alpha=.8)
        ax[0].annotate(v, (d[k], s[k]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    keep = np.array([v != "SEQ" for v in nodes])
    ax[0].set_xlabel("weighted total degree $d$"); ax[0].set_ylabel("state $s$ = P(activated | SEQ)")
    ax[0].set_title(f"(a) marine 39 nodes, corr(d,s) = {pearson(d[keep], s[keep]):.2f} (excl. seed); marker ∝ d·s")
    for a, RR, ttl in ((ax[1], R15, "(b) marine, 15 classes"), (ax[2], RM, "(c) mainland, 15 classes")):
        dd = RR["d_out_w"] + RR["d_in_w"]; ss = RR["p_active"]
        kk = np.array([v != "EQ or AS" for v in RR["nodes"]])
        for v in RR["nodes"]:
            k = RR["idx"][v]
            a.scatter(dd[k], ss[k], color=lay_col[LAYER15[v]], s=40, alpha=.85)
            a.annotate(ABBR15[v], (dd[k], ss[k]), fontsize=8, xytext=(3, 3), textcoords="offset points")
        a.set_xlabel("weighted total degree $d$"); a.set_ylabel("state $s$ = P(activated | EQ)")
        a.set_title(f"{ttl}: corr(d,s) = {pearson(dd[kk], ss[kk]):.2f} (excl. seed)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "figQ4_degree_state.png"), dpi=200); plt.close(fig)

    # Q5: removal search
    rem = summary["marine39"]["removal"]
    base = rem["baseline_E_cons"]
    labels, vals = [], []
    with open(os.path.join(OUT, "marine39_removal_search.csv"), encoding="utf-8") as f:
        rd = list(csv.DictReader(f))
    for k in ("1", "2", "3"):
        for r in [x for x in rd if x["k"] == k][:5]:
            labels.append(f"k={k}: {r['removed_set']}"); vals.append(float(r["reduction_%"]))
    for r in [x for x in rd if x["k"] == "ref"]:
        labels.append(f"ref: {r['removed_set']}"); vals.append(float(r["reduction_%"]))
    fig, ax = plt.subplots(figsize=(9, 7))
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, color=["#d21d13" if l.startswith("k=") else "grey" for l in labels])
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(f"reduction of E[# consequence nodes reached | SEQ] (baseline {base:.2f})  [%]")
    ax.set_title("exhaustive k-node removal over process-layer nodes\n(4 000 common-random-number cascades from SEQ)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "figQ5_removal_search.png"), dpi=200); plt.close(fig)


if __name__ == "__main__":
    main()
