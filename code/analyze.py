"""Preregistered analysis (sections 4 and 6 of preregistration.md). Reads results/<set>/*.json.

Usage:  python code/analyze.py --set main [--results-dir results] [--out-dir results/analysis]

Per term and model (MiniLM primary, bge confirmatory, reported separately):
  L = AUROC(offtarget) - AUROC(target)        (H1)
  R = mean AUROC(rand0, rand1, rand2) - AUROC(target)   (H2)
  D = AUROC(full) - AUROC(target)
For each: median over terms, 95% percentile bootstrap CI of the median (10,000 resamples of terms,
numpy default_rng(0), a fresh generator for every quantity so all quantities use the same resamples),
and a two-sided Wilcoxon signed-rank p (scipy defaults). The MiniLM H1 and H2 p-values are
Holm-corrected together. Material leakage = MiniLM L CI lower bound > 0 AND median >= 0.02.
Secondary and exploratory summaries are labelled as such.
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design_lib as D  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_BOOT = 10000
MODELS = ("minilm", "bge")


def load_results(results_dir, set_name):
    res = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(results_dir, set_name, "*.json")))]
    if any(r["set"] != set_name for r in res):
        raise RuntimeError("result file from another set found")
    return sorted(res, key=lambda r: r["term"])


def au(r, key):
    return r["metrics"][key]["auroc"]


def lrd(r, model):
    t = au(r, f"{model}:target")
    rand = np.mean([au(r, f"{model}:rand{k}") for k in range(3)])
    return au(r, f"{model}:offtarget") - t, float(rand) - t, au(r, f"{model}:full") - t


def bootstrap_ci(x, n_boot=N_BOOT):
    x = np.asarray(x, dtype=float)
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    meds = np.median(x[idx], axis=1)
    lo, hi = np.percentile(meds, [2.5, 97.5])
    return float(lo), float(hi)


def wilcoxon_p(x):
    x = np.asarray(x, dtype=float)
    if np.all(x == 0):
        return float("nan")
    return float(stats.wilcoxon(x).pvalue)


def summarize(x):
    x = np.asarray(x, dtype=float)
    lo, hi = bootstrap_ci(x)
    return {"n": int(len(x)), "median": float(np.median(x)), "ci_low": lo, "ci_high": hi,
            "p_wilcoxon": wilcoxon_p(x), "n_zero": int(np.sum(x == 0)), "mean": float(np.mean(x))}


def holm(pvals):
    """{name: p} -> {name: Holm-adjusted p}."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m, out, running = len(items), {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m - i) * p))
        out[k] = running
    return out


def name_contained(results, rows, text):
    """MiniLM-full percentiles of positives whose summary body contains every content token of the term
    name vs. those that do not (pooled gene-term pairs), with a two-sided Mann-Whitney p."""
    names = {r["term"]: r["name"] for r in rows}
    yes, no, skipped = [], [], []
    for r in results:
        need = set(D.content_tokens(names[r["term"]]))
        if not need:
            skipped.append(r["term"])
            continue
        for p in r["positives"]:
            body = {D.norm(t) for t in D.TOKEN.findall(D.split_prefix(text[p["gene"]])[1])}
            (yes if need <= body else no).append(p["pct_minilm_full"])
    out = {"n_contains": len(yes), "n_not": len(no), "terms_skipped_empty_name": skipped,
           "median_pct_contains": float(np.median(yes)) if yes else None,
           "median_pct_not": float(np.median(no)) if no else None, "p_mannwhitney": None}
    if yes and no:
        out["p_mannwhitney"] = float(stats.mannwhitneyu(yes, no).pvalue)
    return out


def analyze(results, rows, text):
    S = {"n_terms": len(results), "terms": [r["term"] for r in results]}
    per_term, prim = [], {}
    for model in MODELS:
        L, R, Dd = zip(*(lrd(r, model) for r in results))
        prim[model] = {"L": summarize(L), "R": summarize(R), "D": summarize(Dd)}
        for i, r in enumerate(results):
            if model == MODELS[0]:
                per_term.append({"term": r["term"], "name": r["name"],
                                 "target_del_pos": r["deletions"]["target"]["tokens_deleted_pos"]})
            per_term[i].update({f"{model}_L": L[i], f"{model}_R": R[i], f"{model}_D": Dd[i]})
    S["primary"] = prim
    S["holm_minilm"] = holm({"H1_L": prim["minilm"]["L"]["p_wilcoxon"], "H2_R": prim["minilm"]["R"]["p_wilcoxon"]})
    Lm = prim["minilm"]["L"]
    S["material"] = {"verdict": "MATERIAL" if (Lm["ci_low"] > 0 and Lm["median"] >= 0.02) else "NOT MATERIAL",
                     "rule": "MiniLM L CI lower bound > 0 AND median >= 0.02",
                     "ci_low": Lm["ci_low"], "median": Lm["median"]}
    sub = [r for r in results if r["deletions"]["target"]["tokens_deleted_pos"] >= 1]
    S["secondary_restricted_L_minilm"] = summarize([lrd(r, "minilm")[0] for r in sub]) if sub else {"n": 0}
    ex = {}
    tf = lambda r, c: au(r, f"tfidf:{c}")  # noqa: E731
    ex["tfidf_L"] = summarize([tf(r, "offtarget") - tf(r, "target") for r in results])
    ex["tfidf_D"] = summarize([tf(r, "full") - tf(r, "target") for r in results])
    ex["phrase_drop"] = {m: summarize([au(r, f"{m}:full") - au(r, f"{m}:phrase") for r in results]) for m in MODELS}
    keys = sorted(results[0]["metrics"])
    med = {k: float(np.median([au(r, k) for r in results])) for k in keys}
    for m in MODELS:
        med[f"{m}:rand_mean"] = float(np.median([np.mean([au(r, f"{m}:rand{k}") for k in range(3)]) for r in results]))
    ex["median_auroc"] = med
    ex["name_contained_minilm_full"] = name_contained(results, rows, text)
    S["exploratory"] = ex
    return S, per_term


def fmt(x, nd=4):
    return "nan" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


def fmt_p(p):
    return "nan" if p is None or np.isnan(p) else f"{p:.3g}"


def to_text(S):
    L = [f"ANALYSIS ({S['n_terms']} terms)", ""]
    for m in MODELS:
        role = "primary" if m == "minilm" else "confirmatory"
        L += [f"== {m} ({role}) ==", "| quantity | n | median | 95% CI | Wilcoxon p | Holm p | zeros |", "|---|---|---|---|---|---|---|"]
        for q, h in (("L", "H1_L"), ("R", "H2_R"), ("D", None)):
            s = S["primary"][m][q]
            hp = fmt_p(S["holm_minilm"][h]) if (m == "minilm" and h) else "-"
            L.append(f"| {q} | {s['n']} | {fmt(s['median'])} | [{fmt(s['ci_low'])}, {fmt(s['ci_high'])}] | "
                     f"{fmt_p(s['p_wilcoxon'])} | {hp} | {s['n_zero']} |")
        L.append("")
    mat = S["material"]
    L += [f"MATERIAL LEAKAGE: {mat['verdict']}  ({mat['rule']}; CI lower={fmt(mat['ci_low'])}, median={fmt(mat['median'])})", ""]
    s = S["secondary_restricted_L_minilm"]
    if s.get("n"):
        L += ["== SECONDARY: MiniLM L, terms with >=1 target deletion in positives ==",
              f"n={s['n']} median={fmt(s['median'])} CI=[{fmt(s['ci_low'])}, {fmt(s['ci_high'])}] "
              f"Wilcoxon p={fmt_p(s['p_wilcoxon'])} zeros={s['n_zero']}", ""]
    ex = S["exploratory"]
    L += ["== EXPLORATORY (not preregistered hypotheses) ==", "| quantity | n | median | 95% CI | Wilcoxon p |", "|---|---|---|---|---|"]
    for label, s in (("TF-IDF L (offtarget - target)", ex["tfidf_L"]), ("TF-IDF D (full - target)", ex["tfidf_D"]),
                     ("MiniLM phrase drop (full - phrase)", ex["phrase_drop"]["minilm"]),
                     ("bge phrase drop (full - phrase)", ex["phrase_drop"]["bge"])):
        L.append(f"| {label} | {s['n']} | {fmt(s['median'])} | [{fmt(s['ci_low'])}, {fmt(s['ci_high'])}] | {fmt_p(s['p_wilcoxon'])} |")
    L += ["", "Median AUROC per feature set (exploratory):", "| feature set | median AUROC |", "|---|---|"]
    for k, v in ex["median_auroc"].items():
        L.append(f"| {k} | {fmt(v)} |")
    nc = ex["name_contained_minilm_full"]
    L += ["", "MiniLM full: percentile of positives among negatives, by whether the summary contains all content "
          "tokens of the term name (exploratory, pooled gene-term pairs):",
          f"contains all: n={nc['n_contains']} median percentile={fmt(nc['median_pct_contains'], 2)}; "
          f"does not: n={nc['n_not']} median percentile={fmt(nc['median_pct_not'], 2)}; "
          f"Mann-Whitney p={fmt_p(nc['p_mannwhitney'])}; terms skipped (no name content tokens)={len(nc['terms_skipped_empty_name'])}"]
    return "\n".join(L)


def load_text(frozen):
    keys = {}
    with open(os.path.join(frozen, "universe_keys.tsv"), encoding="utf-8") as f:
        f.readline()
        for line in f:
            g, sk, _, _ = line.rstrip("\n").split("\t")
            keys[g] = sk
    jp, _ = D.pick_genept_files(os.path.join(ROOT, "data", "raw", "genept"))
    summ = json.load(open(jp, encoding="utf-8"))
    return {g: summ[k] for g, k in keys.items()}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="main")
    ap.add_argument("--results-dir", default=os.path.join(ROOT, "results"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "results", "analysis"))
    args = ap.parse_args(argv)
    frozen = os.path.join(ROOT, "frozen")
    rows = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == args.set]
    results = load_results(args.results_dir, args.set)
    S, per_term = analyze(results, rows, load_text(frozen))
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, f"analysis_{args.set}.json"), "w") as f:
        json.dump(S, f, indent=1)
    with open(os.path.join(args.out_dir, f"per_term_{args.set}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_term[0]))
        w.writeheader()
        w.writerows(per_term)
    txt = to_text(S)
    with open(os.path.join(args.out_dir, f"analysis_{args.set}.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
