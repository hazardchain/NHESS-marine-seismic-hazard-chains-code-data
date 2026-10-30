# Quantitative extension: from a binary coupling network to a probabilistic cascade model

This folder turns the 39-node marine hazard-chain network into a probabilistic
model **using only the event-level information already contained in `data/`**.
Nothing is assumed about edge weights: every number below is estimated from the
58 documented source events (366 chains) or derived from those estimates.

It answers the central criticism of RC1 on egusphere-2026-4167 ("no transition
rate is estimated, no probability is calculated, no model is calibrated") and
adds the degree–state resilience measures of Zheng et al. (2026, *Nat. Commun.*,
"Discovering network resilience formulas with symbolized reinforcement learning").

```
python -m venv .venv && .venv/bin/pip install numpy matplotlib
.venv/bin/python -I quantitative/run_quantitative.py        # ~4 min, seed fixed
```

Everything is regenerated into `outputs/` (CSV, JSON) and `figures/` (PNG).

## 1. Edge triggering probabilities

For every documented coupling i → j

> p_ij = (number of source events in which i → j is documented) / (number of source events in which i occurs)

with a Wilson 95 % interval and a Jeffreys-shrunk alternative
(`outputs/marine39_edge_probabilities.csv`, Fig. Q1). The self-loop PP → PP is
kept in the record but carries no propagation.

| quantity | value |
|---|---|
| propagating edges | 151 |
| edges documented in **one** event only | 80 (53 %) |
| median p_ij (IQR) | 0.15 (0.05 – 0.33) |
| most probable couplings | MTD→BI 1.00 (10/10), SLS→TS 0.98 (39/40), CF→BI 0.91 (42/46), BI→EC 0.89 (42/47), CF→PP 0.83 (38/46), SEQ→CF 0.73 (36/49), CFL→CF 0.71 (17/24), SEQ→SLS 0.67 (33/49) |

Interpretation and caveats. p_ij is the frequency with which j is *documented
downstream of* i, conditional on i being documented at all. It is a documented
co-occurrence frequency, not a mechanistic triggering probability, and inherits
the reporting bias of the record. Edges seen once in a node that occurs once
(e.g. NFA→WP, 1/1) get p = 1 with an interval of 0.21 – 1.00; the interval, not
the point estimate, is the honest statement for such edges. The coastal
earthquake EQ appears in only 4 of the 58 tsunami events, so its outgoing
probabilities are the least certain in the network.

## 2. Weighted transition matrix and cumulative interaction matrix

With these probabilities the paper's Eqs. (10)–(12) become computable:

| quantity | value |
|---|---|
| spectral radius ρ(A) | 0.48 (event bootstrap 95 %: 0.23 – 0.75), series converges |
| ρ(A) with the uniform 0.1 weight of the first manuscript | 0.27 |
| largest row sum | 4.75 (EQ); rows are *not* probabilities summing to ≤ 1, because one trigger can document several downstream hazards in the same event; A is a branching-process offspring matrix, not a substochastic Markov matrix, which settles the Eq. (8)–(9) inconsistency raised by RC1 |

`outputs/marine39_M_cumulative.csv` holds M = (I − A)⁻¹ − I. The entry m_ij is the
expected number of documented-type pathways from i to j summed over all
lengths. From the submarine earthquake: BI 3.59, EC 3.24, PP 2.82, CF 2.38,
TS 1.25, NF 0.80 (`outputs/marine39_M_from_sources.csv`). The three consequence
nodes lead not because they are directly triggered (their direct A entries are
≈ 0) but because almost every pathway ends in them: this is the "convergence"
reading of the first manuscript, now as a number rather than an in-degree.

## 3. Independent-cascade Monte Carlo from the submarine earthquake

Each edge fires independently with probability p_ij (live-edge sampling,
20 000 runs). `outputs/marine39_activation_from_SEQ.csv`, Figs. Q2–Q3.

| quantity | value (event-bootstrap 95 % where computed) |
|---|---|
| cascade size, nodes activated besides SEQ | mean 9.3, median 9, 90 % range 4 – 15 |
| expected consequence-layer nodes reached | 3.31 (2.84 – 3.79) |
| P(at least one consequence node) | 0.95 |
| P(critical infrastructure failure) | 0.92 (0.86 – 0.96) |
| P(business interruption) | 0.91 |
| P(physical & psychological trauma) | 0.84 (0.74 – 0.92) |
| P(tsunami) | 0.82 (0.74 – 0.89) |
| P(critical network failure) | 0.46 (0.29 – 0.59) |
| P(nuclear facility accident) | 0.04 (0.00 – 0.12) |

This is the "input → output" model the first manuscript did not have: the input
is the seed hazard, the output is a probability over every downstream hazard
and a distribution of cascade size.

## 4. Node state, d·s and corr(d, s) (Zheng et al. 2026)

Zheng et al. find that a node's contribution to resilience is the product of its
degree d and its dynamical state s, and that networks whose degree and state
are *less* correlated are *more* resilient. Here the state is defined without
any new assumption as

* s_act = P(node activated | SEQ) from the Monte Carlo (closest analogue of their steady state), and
* s_cons = expected number of consequence-layer nodes reached when the cascade is seeded at the node (downstream consequence potential).

`outputs/marine39_node_state_ds.csv`, Fig. Q4a.

| rank by d_tot,w · s_act | node | d_out,w | d_in,w | s_act | d·s |
|---|---|---|---|---|---|
| 1 | CF critical infrastructure failure | 2.41 | 5.26 | 0.92 | 7.02 |
| 2 | BI business interruption | 0.89 | 4.87 | 0.91 | 5.24 |
| 3 | TS tsunami | 1.13 | 3.49 | 0.82 | 3.77 |
| 4 | SEQ submarine earthquake (seed) | 3.55 | 0 | 1 | 3.55 |
| 5 | NF critical network failure | 1.00 | 3.85 | 0.46 | 2.24 |

CF is the top d·s node in 298 of 300 event-bootstrap replicates. Note that
d·s is *not* reducible to degree: SEQ has by far the largest out-degree yet
ranks fourth, and CFL (weighted out-degree 2.21, second only to SEQ and CF)
ranks eighth because its activation probability is 0.43.

corr(d, s) excluding the seed (`outputs/marine39_corr_d_s.csv`):

| degree used | corr(d, s_act) | corr(d, s_cons) |
|---|---|---|
| weighted out-degree | 0.03 | 0.79 |
| weighted in-degree | 0.82 | 0.25 |
| weighted total degree | 0.66 | 0.63 |

The directed structure splits Zheng's single correlation in two: in-degree
predicts *whether a node is reached* (s_act), out-degree predicts *what it
leads to* (s_cons). Neither alone is the full picture, which is the directed
analogue of their d·s argument.

### Marine vs mainland at 15 classes

Both networks are rebuilt from their own event records with the same pipeline
(marine chains aggregated to the 15 classes before counting; mainland 94 events
after merging the split ID rows of the regional tables).

| | marine, 15 classes | mainland, 15 classes |
|---|---|---|
| events / edges | 58 / 57 | 94 / 47 |
| ρ(A) | 0.47 | 0.20 |
| mean cascade size from EQ | 6.7 | 4.4 |
| expected consequence classes reached | 3.0 | 1.8 |
| corr(d_tot,w , s_act), excl. seed | **0.76** (0.63 – 0.85) | **0.91** (0.77 – 0.96) |
| bootstrap difference marine − mainland | −0.12 (−0.26 – +0.01); marine lower in 95 % of replicates | |

Fig. Q4b–c. In Zheng et al.'s terms the mainland chain is the "man-made"
pattern, where degree and state coincide on a few hubs; the marine chain spreads
state across nodes of moderate degree (flood, mass slide, trauma). This is the
dynamical counterpart of the first manuscript's centralisation contrast
(0.83 vs 0.55), with an uncertainty statement attached. The interval touches
zero, so the difference should be reported as suggestive, not established.

## 5. Exact intervention search (replaces "optimal intervention")

Candidates are the 27 process-layer nodes (natural + engineering). Triggers
cannot be prevented, and removing a consequence node would delete the very
outcome being counted (the automatic-reduction bias of RC1). Objective:
expected number of consequence-layer nodes reached from SEQ, evaluated on
4 000 common-random-number cascades. All 27 + 351 + 2 925 sets are enumerated;
no heuristic and no learning is needed at this network size.
`outputs/marine39_removal_search.csv`, Fig. Q5.

| removed | E[# consequence nodes] | reduction | P(any consequence) |
|---|---|---|---|
| none | 3.29 | | 0.95 |
| CF | 1.59 | **−51.7 %** | 0.65 |
| TS | 2.96 | −10.1 % | 0.91 |
| SLS | 3.08 | −6.4 % | 0.92 |
| NF | 3.15 | −4.3 % | 0.94 |
| FI | 3.17 | −3.5 % | 0.95 |
| CF + TS (best pair) | 0.86 | −74.0 % | 0.39 |
| CF + NF + TS (best triple) | 0.55 | **−83.4 %** | 0.28 |
| CF + FI + TS (first manuscript's greedy set) | 0.73 | −77.8 % | 0.35 |

CF is the best single removal in 300 of 300 bootstrap replicates. The optimal
triple is {CF, NF, TS}, not the manuscript's {CF, FI, TS}; the same correction
was obtained earlier from the fixed-node-set reachability analysis, so the
"fire outperforms network failure" statement should be withdrawn. Note also
that NF alone is weak (−4 %) but strong *in combination* with CF, which is
exactly the convergence-node role the manuscript assigns to it.

## 6. Held-out events, probability-weighted

The three events outside the record (2004 Indian Ocean, 1960 Valdivia,
2018 Palu) contribute 54 distinct couplings, 40 of which are in the network
(the 74 % of the manuscript). With probabilities attached the statement
becomes sharper (`outputs/heldout_events_edge_probabilities.csv`):

| | value |
|---|---|
| mean p_ij of the covered held-out couplings | 0.37 |
| mean p_ij over all 151 network edges | 0.25 |
| covered held-out couplings above the network median p_ij | 67 % |
| covered held-out coupling rows that are single-event edges | 21 % (vs 53 % of all edges) |

The held-out events reproduce the probable backbone, not a random 74 % of the
edge list; the 14 novel couplings are the low-probability tail. This is the
defensible version of the "backbone is reproduced" claim.

## What this does and does not fix

Fixed, relative to the first manuscript:
* transition rates are estimated (Sect. 1), the Markov-vs-adjacency ambiguity of A is resolved (Sect. 2), M is computed (Sect. 2);
* frequency information is no longer discarded (every edge carries n_ij / n_i);
* "amplification" is replaced by activation probability and d·s (Sects. 3–4);
* "optimal intervention" is replaced by an exhaustive search with a stated objective and bootstrap stability (Sect. 5);
* the hindcast is probability-weighted (Sect. 6).

Still open, and to be stated as such in a revised manuscript:
* p_ij is a documented co-occurrence frequency under reporting bias; 53 % of edges rest on one event;
* the independent-cascade model assumes edges fire independently and ignores timing;
* the node state is an activation probability, not a calibrated loss. Replacing s by a loss-calibrated state (fatality / economic-loss ratio per hazard at a reference intensity) is the next step and requires exposure inputs that are not in this repository;
* the intervention search has no cost term; it ranks topological-probabilistic leverage, not investment efficiency;
* the taxonomy and the coding remain expert judgements; the bootstrap propagates sampling uncertainty only.
