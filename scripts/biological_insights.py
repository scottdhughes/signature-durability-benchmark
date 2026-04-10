"""
biological_insights.py
======================
Mine per_cohort_effects.csv for novel biological patterns that could make
the signature-durability-benchmark paper publishable at a top venue.

Five analyses:
  1. CROSS-TALK  - signatures showing signal in programs they don't belong to
  2. ANTI-CORRELATION - Hallmark program correlation matrix across cohorts
  3. UNIVERSAL vs CONTEXT-SPECIFIC cohorts (permissive vs restrictive)
  4. SIGNATURE SIMILARITY - cosine similarity of 30-cohort effect vectors;
     do stealth-confounded signatures cluster with their real vs confound program?
  5. COHORT QUALITY INDEX - mean |d| across durable signatures per cohort
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cosine
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────────────
EFFECTS_CSV  = "outputs/canonical_v6/per_cohort_effects.csv"
MANIFEST_TSV = "config/cohort_manifest.tsv"

# ── load data ──────────────────────────────────────────────────────────────────
df   = pd.read_csv(EFFECTS_CSV)
mf   = pd.read_csv(MANIFEST_TSV, sep="\t")

# Build wide matrix: rows = cohorts, cols = signatures, values = Cohen's d
wide = df.pivot(index="cohort_id", columns="signature_id", values="cohens_d")

# Map cohort → biological program
cohort_prog = mf.set_index("cohort_id")["biological_program"].to_dict()
# Map signature → "canonical" program (by naming convention)
def sig_program(sig):
    if "ifng" in sig or "ifna" in sig or "blind_durable_ifn" in sig:
        return "interferon"
    if "inflam" in sig or "tnfa" in sig or "brittle_small_study_inflam" in sig:
        return "inflammation"
    if "hypoxia" in sig or "blind_durable_hypoxia" in sig:
        return "hypoxia"
    if "emt" in sig and "inflam" not in sig:
        return "emt"
    if "e2f" in sig or "prolif" in sig or "mixed_senescence" in sig:
        return "proliferation"
    if "senescence" in sig:
        return "proliferation"
    return "other"

sig_prog = {s: sig_program(s) for s in wide.columns}

PROGRAMS = ["interferon", "inflammation", "hypoxia", "emt", "proliferation"]

HALLMARKS = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "hallmark_inflammatory_response",
    "hallmark_hypoxia",
    "hallmark_e2f_targets",
    "hallmark_emt",
    "hallmark_tnfa_nfkb",
]

DURABLE_SIGS = [
    "hallmark_ifng_response",
    "hallmark_ifna_response",
    "hallmark_inflammatory_response",
    "hallmark_hypoxia",
    "hallmark_e2f_targets",
    "hallmark_emt",
    "hallmark_tnfa_nfkb",
    "curated_senescence",
    "blind_durable_ifn_composite",
    "blind_durable_hypoxia_core",
]

SEP = "=" * 72

# ══════════════════════════════════════════════════════════════════════════════
# 1. CROSS-TALK ANALYSIS
#    For each signature, look at its mean |Cohen's d| in cohorts whose program
#    differs from the signature's own program.  High off-program effect = cross-talk.
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("ANALYSIS 1 — CROSS-TALK: Signatures active outside their native program")
print(SEP)

rows = []
for sig in HALLMARKS:
    my_prog = sig_prog[sig]
    for cohort in wide.index:
        coh_prog = cohort_prog.get(cohort, "unknown")
        d = wide.loc[cohort, sig]
        rows.append({
            "signature": sig,
            "sig_program": my_prog,
            "cohort": cohort,
            "cohort_program": coh_prog,
            "cohens_d": d,
            "abs_d": abs(d),
            "cross_talk": coh_prog != my_prog,
        })

ct = pd.DataFrame(rows)

# Mean |d| in own program vs other programs
print("\n--- Mean |Cohen's d| in own vs foreign programs ---")
summary = (
    ct.groupby(["signature", "cross_talk"])["abs_d"]
    .agg(["mean", "max", "count"])
    .rename(columns={"mean": "mean_abs_d", "max": "max_abs_d", "count": "n_cohorts"})
)
summary.index.names = ["signature", "is_cross_talk"]
print(summary.to_string())

# Top cross-talk: signature × cohort combos with large off-program effects
print("\n--- Top 20 off-program (cross-talk) effects (|d| > 0.5) ---")
top_ct = (
    ct[ct["cross_talk"] & (ct["abs_d"] > 0.5)]
    .sort_values("abs_d", ascending=False)
    .head(20)[["signature", "sig_program", "cohort", "cohort_program", "cohens_d"]]
)
print(top_ct.to_string(index=False))

# Surprising specific finding: EMT signature in inflammation cohorts
print("\n--- EMT signature in inflammation cohorts ---")
emt_in_inflam = ct[(ct["signature"] == "hallmark_emt") & (ct["cohort_program"] == "inflammation")]
print(emt_in_inflam[["cohort", "cohort_program", "cohens_d"]].to_string(index=False))

# IFN signatures in hypoxia cohorts
print("\n--- IFN signatures in hypoxia cohorts ---")
ifn_in_hyp = ct[ct["signature"].isin(["hallmark_ifng_response","hallmark_ifna_response"])
                & (ct["cohort_program"] == "hypoxia")]
print(ifn_in_hyp[["signature","cohort","cohort_program","cohens_d"]].to_string(index=False))

# Proliferation (E2F) in inflammation cohorts
print("\n--- E2F/proliferation in inflammation cohorts ---")
e2f_in_inflam = ct[(ct["signature"] == "hallmark_e2f_targets") & (ct["cohort_program"] == "inflammation")]
print(e2f_in_inflam[["cohort","cohort_program","cohens_d"]].to_string(index=False))

# ══════════════════════════════════════════════════════════════════════════════
# 2. ANTI-CORRELATION: Hallmark program correlation matrix across 30 cohorts
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + SEP)
print("ANALYSIS 2 — ANTI-CORRELATION: Program correlation matrix (30 cohorts)")
print(SEP)

hall_wide = wide[HALLMARKS].dropna()

print(f"\nUsing {len(hall_wide)} cohorts with complete Hallmark data")

# Pearson correlation matrix
corr_mat = hall_wide.corr(method="pearson")
print("\n--- Pearson correlation matrix ---")
print(corr_mat.round(3).to_string())

print("\n--- Strongly ANTI-correlated pairs (r < -0.3) ---")
anti_pairs = []
for a, b in combinations(HALLMARKS, 2):
    r = corr_mat.loc[a, b]
    if r < -0.3:
        anti_pairs.append((a, b, r))
anti_pairs.sort(key=lambda x: x[2])
if anti_pairs:
    for a, b, r in anti_pairs:
        print(f"  {a}  ×  {b}  →  r = {r:.3f}")
else:
    print("  None found at r < -0.3")

print("\n--- Strongly CORRELATED pairs (r > 0.5, same program) ---")
for a, b in combinations(HALLMARKS, 2):
    r = corr_mat.loc[a, b]
    if r > 0.5 and sig_prog[a] == sig_prog[b]:
        print(f"  {a}  ×  {b}  →  r = {r:.3f}  [same program: {sig_prog[a]}]")

print("\n--- Strongly CORRELATED pairs (r > 0.5, DIFFERENT program) [cross-program co-activation] ---")
cross_prog_corr = []
for a, b in combinations(HALLMARKS, 2):
    r = corr_mat.loc[a, b]
    if r > 0.5 and sig_prog[a] != sig_prog[b]:
        cross_prog_corr.append((a, sig_prog[a], b, sig_prog[b], r))
cross_prog_corr.sort(key=lambda x: -x[4])
if cross_prog_corr:
    for row in cross_prog_corr:
        print(f"  {row[0]} [{row[1]}]  ×  {row[2]} [{row[3]}]  →  r = {row[4]:.3f}")
else:
    print("  None found at r > 0.5")

# Spearman as robustness check
spear_mat = hall_wide.corr(method="spearman")
print("\n--- Spearman correlation matrix (robustness check) ---")
print(spear_mat.round(3).to_string())

# Highlight the single most anti-correlated pair
flat = []
for a, b in combinations(HALLMARKS, 2):
    flat.append((a, b, corr_mat.loc[a, b], spear_mat.loc[a, b]))
flat.sort(key=lambda x: x[2])
print("\n--- Most anti-correlated pair overall ---")
a, b, rp, rs = flat[0]
print(f"  {a}  ×  {b}:  Pearson r = {rp:.3f},  Spearman r = {rs:.3f}")

# ══════════════════════════════════════════════════════════════════════════════
# 3. UNIVERSAL vs CONTEXT-SPECIFIC COHORTS
#    "Permissive" cohort: many durable signatures show large effect
#    "Restrictive" cohort: most signatures show small effect
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + SEP)
print("ANALYSIS 3 — UNIVERSAL vs CONTEXT-SPECIFIC COHORTS")
print(SEP)

durable_avail = [s for s in DURABLE_SIGS if s in wide.columns]
dur_wide = wide[durable_avail]

# For each cohort: fraction of durable sigs with |d| > 0.5
THRESH = 0.5
cohort_stats = []
for cohort in dur_wide.index:
    row = dur_wide.loc[cohort]
    n_total = row.notna().sum()
    n_active = (row.abs() > THRESH).sum()
    mean_abs_d = row.abs().mean()
    frac_active = n_active / n_total if n_total > 0 else np.nan
    cohort_stats.append({
        "cohort": cohort,
        "program": cohort_prog.get(cohort, "unknown"),
        "n_durable_sigs": n_total,
        "n_active_sigs": n_active,
        "frac_active": frac_active,
        "mean_abs_d": mean_abs_d,
    })

cs_df = pd.DataFrame(cohort_stats).sort_values("frac_active", ascending=False)
print("\n--- Cohort permissiveness (fraction of durable sigs with |d| > 0.5) ---")
print(cs_df.to_string(index=False))

print("\n--- TOP 5 PERMISSIVE cohorts (most signatures work) ---")
print(cs_df.head(5).to_string(index=False))

print("\n--- TOP 5 RESTRICTIVE cohorts (fewest signatures work) ---")
print(cs_df.tail(5).to_string(index=False))

# Within-program comparison: which program's cohorts are most permissive?
print("\n--- Mean permissiveness by biological program ---")
prog_perm = cs_df.groupby("program")["frac_active"].agg(["mean","std","count"])
print(prog_perm.sort_values("mean", ascending=False).round(3).to_string())

# Identify cohorts that activate signatures from MULTIPLE foreign programs
print("\n--- Cross-program activators: cohorts with |d|>0.5 for sigs from ≥3 different programs ---")
for cohort in wide.index:
    progs_active = set()
    for sig in HALLMARKS:
        if abs(wide.loc[cohort, sig]) > 0.5:
            progs_active.add(sig_prog[sig])
    if len(progs_active) >= 3:
        prog_str = ", ".join(sorted(progs_active))
        print(f"  {cohort} [{cohort_prog.get(cohort,'?')}]: activates {prog_str}")

# ══════════════════════════════════════════════════════════════════════════════
# 4. SIGNATURE SIMILARITY via cosine similarity of 30-cohort effect vectors
#    Focus: do stealth-confounded signatures cluster with their real program
#    or with their confounder?
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + SEP)
print("ANALYSIS 4 — SIGNATURE SIMILARITY (cosine similarity of effect vectors)")
print(SEP)

# Drop cohorts with any NaN
full_wide = wide.dropna(axis=0)  # cohorts complete across all sigs
sigs_complete = wide.dropna(axis=1)  # sigs complete across all cohorts

print(f"Cohorts with complete data (all sigs): {len(full_wide)}")
print(f"Signatures with complete data (all cohorts): {sigs_complete.shape[1]}")

# Use all cohorts, impute missing with 0 for cosine (or use only cohorts present for both sigs)
def cosine_sim(v1, v2):
    mask = ~(np.isnan(v1) | np.isnan(v2))
    if mask.sum() < 5:
        return np.nan
    return 1 - cosine(v1[mask], v2[mask])

sigs = list(wide.columns)
n = len(sigs)
sim_mat = np.full((n, n), np.nan)
for i in range(n):
    for j in range(n):
        v1 = wide.iloc[:, i].values.astype(float)
        v2 = wide.iloc[:, j].values.astype(float)
        sim_mat[i, j] = cosine_sim(v1, v2)

sim_df = pd.DataFrame(sim_mat, index=sigs, columns=sigs)

print("\n--- Top 15 most SIMILAR signature pairs (cosine sim) ---")
pairs = []
for i, a in enumerate(sigs):
    for j, b in enumerate(sigs):
        if j > i:
            s = sim_df.loc[a, b]
            if not np.isnan(s):
                pairs.append((a, sig_prog[a], b, sig_prog[b], s))
pairs.sort(key=lambda x: -x[4])
for row in pairs[:15]:
    same = "SAME" if row[1] == row[3] else "DIFF"
    print(f"  [{same}] {row[0]} ({row[1]}) × {row[2]} ({row[3]}): sim={row[4]:.3f}")

print("\n--- Top 10 most DISSIMILAR signature pairs (cosine sim) ---")
anti = sorted(pairs, key=lambda x: x[4])
for row in anti[:10]:
    same = "SAME" if row[1] == row[3] else "DIFF"
    print(f"  [{same}] {row[0]} ({row[1]}) × {row[2]} ({row[3]}): sim={row[4]:.3f}")

# KEY QUESTION: stealth-confounded signatures
STEALTH_SIGS = {
    "stealth_confounded_hypoxia_prolif": {"true": "hypoxia", "confound": "proliferation"},
    "stealth_confounded_inflam_immune":  {"true": "inflammation", "confound": "interferon"},
    "blind_confounded_stealth_emt_immune": {"true": "emt", "confound": "interferon"},
    "blind_confounded_prolif_inflam":    {"true": "proliferation", "confound": "inflammation"},
}

print("\n--- STEALTH CONFOUND CLUSTERING ANALYSIS ---")
print("Do stealth-confounded sigs cluster with their real program or their confounder?\n")

# For each stealth sig, compute its average cosine similarity to sigs of each program
for stealth_sig, meta in STEALTH_SIGS.items():
    if stealth_sig not in sigs:
        print(f"  {stealth_sig}: NOT IN DATA")
        continue
    true_prog = meta["true"]
    confound_prog = meta["confound"]

    prog_sims = {}
    for prog in PROGRAMS:
        prog_sig_list = [s for s in sigs if sig_prog[s] == prog and s != stealth_sig]
        if prog_sig_list:
            sim_vals = [sim_df.loc[stealth_sig, s] for s in prog_sig_list if not np.isnan(sim_df.loc[stealth_sig, s])]
            prog_sims[prog] = np.mean(sim_vals) if sim_vals else np.nan
        else:
            prog_sims[prog] = np.nan

    nearest = max(prog_sims, key=lambda p: prog_sims[p] if not np.isnan(prog_sims.get(p, np.nan)) else -999)

    print(f"  {stealth_sig}")
    print(f"    True program: {true_prog}  |  Confound program: {confound_prog}")
    for prog, s in sorted(prog_sims.items(), key=lambda x: -(x[1] if not np.isnan(x[1]) else -999)):
        marker = " <-- NEAREST" if prog == nearest else ""
        if prog == true_prog:
            marker += " [TRUE]"
        if prog == confound_prog:
            marker += " [CONFOUND]"
        print(f"    {prog}: {s:.3f}{marker}")

    if nearest == true_prog:
        verdict = "CLUSTERS WITH TRUE PROGRAM (well-designed)"
    elif nearest == confound_prog:
        verdict = "*** CLUSTERS WITH CONFOUNDER *** (biologically misleading)"
    else:
        verdict = f"Clusters with {nearest} (unexpected)"
    print(f"    Verdict: {verdict}\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5. COHORT QUALITY INDEX
#    For each cohort, mean |d| across all durable signatures
#    High quality = many durable sigs show strong consistent effects
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + SEP)
print("ANALYSIS 5 — COHORT QUALITY INDEX")
print(SEP)

# Quality = mean |Cohen's d| across durable signatures
cqi = {}
for cohort in wide.index:
    vals = wide.loc[cohort, durable_avail].dropna()
    cqi[cohort] = {
        "program": cohort_prog.get(cohort, "unknown"),
        "mean_abs_d": vals.abs().mean(),
        "max_abs_d": vals.abs().max(),
        "n_sigs_gt1": (vals.abs() > 1.0).sum(),
        "n_sigs_gt0_5": (vals.abs() > 0.5).sum(),
        "median_abs_d": vals.abs().median(),
    }

cqi_df = pd.DataFrame(cqi).T.sort_values("mean_abs_d", ascending=False)
cqi_df = cqi_df.astype({c: float for c in cqi_df.columns if c != "program"})

print("\n--- Cohort Quality Index (ranked by mean |d| across durable signatures) ---")
print(cqi_df.round(3).to_string())

print("\n--- HIGH QUALITY cohorts (mean |d| > 0.7) ---")
hq = cqi_df[cqi_df["mean_abs_d"] > 0.7]
print(hq[["program","mean_abs_d","n_sigs_gt1","n_sigs_gt0_5"]].to_string())

print("\n--- LOW QUALITY / NOISY cohorts (mean |d| < 0.3) ---")
lq = cqi_df[cqi_df["mean_abs_d"] < 0.3]
print(lq[["program","mean_abs_d","n_sigs_gt1","n_sigs_gt0_5"]].to_string())

# Compare platform quality
print("\n--- Quality by platform ---")
plat_map = mf.set_index("cohort_id")["platform"].to_dict()
cqi_df["platform"] = cqi_df.index.map(plat_map)
plat_q = cqi_df.groupby("platform")["mean_abs_d"].agg(["mean","std","count"])
print(plat_q.sort_values("mean", ascending=False).round(3).to_string())

# Compare tissue / biological program quality
print("\n--- Quality by biological program ---")
prog_q = cqi_df.groupby("program")["mean_abs_d"].agg(["mean","std","count"])
print(prog_q.sort_values("mean", ascending=False).round(3).to_string())

# ══════════════════════════════════════════════════════════════════════════════
# BONUS: Outlier signature — hallmark_tnfa_nfkb in EMT cohorts
#         and other non-obvious cross-program observations
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + SEP)
print("BONUS — SURPRISING INDIVIDUAL OBSERVATIONS")
print(SEP)

print("\n--- TNF-α/NF-κB (inflammation) in EMT cohorts ---")
tnfa_emt = df[(df["signature_id"] == "hallmark_tnfa_nfkb") &
              (df["cohort_id"].isin(mf[mf["biological_program"] == "emt"]["cohort_id"]))]
print(tnfa_emt[["cohort_id","cohens_d","direction_consistent"]].to_string(index=False))

print("\n--- Hypoxia signature (hallmark_hypoxia) in INTERFERON cohorts ---")
hyp_ifn = df[(df["signature_id"] == "hallmark_hypoxia") &
             (df["cohort_id"].isin(mf[mf["biological_program"] == "interferon"]["cohort_id"]))]
print(hyp_ifn[["cohort_id","cohens_d","direction_consistent"]].to_string(index=False))

print("\n--- IFN-γ response in HYPOXIA cohorts (HIF1a/IFN cross-talk?) ---")
ifng_hyp = df[(df["signature_id"] == "hallmark_ifng_response") &
              (df["cohort_id"].isin(mf[mf["biological_program"] == "hypoxia"]["cohort_id"]))]
print(ifng_hyp[["cohort_id","cohens_d","direction_consistent"]].to_string(index=False))

print("\n--- Cohorts where IFN and PROLIFERATION are simultaneously UP (d>0.5 both) ---")
for cohort in wide.index:
    ifn = wide.loc[cohort, "hallmark_ifng_response"]
    e2f = wide.loc[cohort, "hallmark_e2f_targets"]
    if ifn > 0.5 and e2f > 0.5:
        print(f"  {cohort} [{cohort_prog.get(cohort,'?')}]: IFNg_d={ifn:.2f}, E2F_d={e2f:.2f}")

print("\n--- Cohorts where EMT and INFLAMMATION are simultaneously UP (d>0.5 both) ---")
for cohort in wide.index:
    emt_d = wide.loc[cohort, "hallmark_emt"]
    inf_d = wide.loc[cohort, "hallmark_inflammatory_response"]
    if emt_d > 0.5 and inf_d > 0.5:
        print(f"  {cohort} [{cohort_prog.get(cohort,'?')}]: EMT_d={emt_d:.2f}, Inflam_d={inf_d:.2f}")

# ── Summary of surprising findings ────────────────────────────────────────────
print("\n" + SEP)
print("SUMMARY OF POTENTIALLY NOVEL / PUBLISHABLE FINDINGS")
print(SEP)

# 1) Most extreme cross-talk event
top_ct_row = ct[ct["cross_talk"]].sort_values("abs_d", ascending=False).iloc[0]
print(f"\n[1] EXTREME CROSS-TALK: {top_ct_row['signature']} (a {top_ct_row['sig_program']} sig)"
      f" scores d={top_ct_row['cohens_d']:.2f} in {top_ct_row['cohort']} (a {top_ct_row['cohort_program']} cohort)")

# 2) Most anti-correlated Hallmark pair
print(f"\n[2] STRONGEST ANTI-CORRELATION: {a} × {b} (Pearson r={rp:.3f}, Spearman r={rs:.3f})")

# 3) Most permissive vs most restrictive cohort
print(f"\n[3] MOST PERMISSIVE cohort: {cs_df.iloc[0]['cohort']} ({cs_df.iloc[0]['program']}, "
      f"{cs_df.iloc[0]['frac_active']:.0%} of durable sigs active)")
print(f"    MOST RESTRICTIVE cohort: {cs_df.iloc[-1]['cohort']} ({cs_df.iloc[-1]['program']}, "
      f"{cs_df.iloc[-1]['frac_active']:.0%} of durable sigs active)")

# 4) Stealth confound clustering
print("\n[4] STEALTH CONFOUND: see detailed analysis above for clustering verdicts")

# 5) Highest quality cohort
print(f"\n[5] HIGHEST QUALITY cohort: {cqi_df.index[0]} ({cqi_df.iloc[0]['program']}, "
      f"mean |d|={cqi_df.iloc[0]['mean_abs_d']:.3f})")
print(f"    LOWEST QUALITY cohort:  {cqi_df.index[-1]} ({cqi_df.iloc[-1]['program']}, "
      f"mean |d|={cqi_df.iloc[-1]['mean_abs_d']:.3f})")

print("\nDone.\n")
