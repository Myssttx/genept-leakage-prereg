"""Exploratory temporal holdout (preregistration Amendment 1) for the main terms marked temporal-eligible
in frozen/design.json.

Usage:  python code/temporal.py [--out-dir results/analysis]
Train on baseline pos (1) + neg (0) with the term's new_pos genes removed from the negatives; test on
new_pos (1) + tneg (0). Same model as CV: StandardScaler + LogisticRegression(C=1, class_weight="balanced",
max_iter=5000), no scaler for TF-IDF. Feature sets: minilm full, minilm target, ada, gene2vec, tfidf full.
Reports test AUROC per term, medians, and the ada - gene2vec gap in CV (results/main) vs temporal.
"""
import argparse
import csv
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design_lib as D  # noqa: E402
import evaluate as E  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURES = ["minilm:full", "minilm:target", "ada:full", "gene2vec", "tfidf:full"]


def split_sets(d):
    """(train genes, train y, test genes, test y) for one term's design entry."""
    newp = set(d["new_pos"])
    neg = [g for g in d["neg"] if g not in newp]
    train = list(d["pos"]) + neg
    test = list(d["new_pos"]) + list(d["tneg"])
    return (train, [1] * len(d["pos"]) + [0] * len(neg), test, [1] * len(d["new_pos"]) + [0] * len(d["tneg"]),
            len(d["neg"]) - len(neg))


def fit_predict_auroc(Xtr, ytr, Xte, yte, scale=True):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
    model = make_pipeline(StandardScaler(), lr) if scale else lr
    model.fit(Xtr, np.asarray(ytr))
    return float(roc_auc_score(np.asarray(yte), model.predict_proba(Xte)[:, 1]))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "results", "analysis"))
    args = ap.parse_args(argv)
    from sklearn.feature_extraction.text import TfidfVectorizer
    frozen = os.path.join(ROOT, "frozen")
    rows = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == "main"]
    design = json.load(open(os.path.join(frozen, "design.json")))
    elig = [r for r in rows if design[r["term"]]["eligible"]]
    keys = {}
    with open(os.path.join(frozen, "universe_keys.tsv"), encoding="utf-8") as f:
        f.readline()
        for line in f:
            g, sk, ak, gk = line.rstrip("\n").split("\t")
            keys[g] = (sk, ak, gk)
    universe = sorted(keys)
    jp, ap_ = D.pick_genept_files(os.path.join(ROOT, "data", "raw", "genept"))
    summ = json.load(open(jp, encoding="utf-8"))
    text = {g: summ[keys[g][0]] for g in universe}

    need = sorted({g for r in elig for k in ("pos", "neg", "new_pos", "tneg") for g in design[r["term"]][k]})
    with open(ap_, "rb") as f:
        ada_all = pickle.load(f)
    ada = {g: np.asarray(ada_all[keys[g][1]], dtype=np.float32).ravel() for g in need}
    del ada_all
    g2v_need = {keys[g][2]: g for g in need}
    g2v = {}
    with open(os.path.join(ROOT, "data", "raw", "gene2vec_dim_200_iter_9.txt"), encoding="utf-8") as f:
        for line in f:
            p = line.split()
            if p and p[0] in g2v_need:
                g2v[g2v_need[p[0]]] = np.asarray(p[1:], dtype=np.float32)
    tfidf = TfidfVectorizer(max_features=20000, min_df=2, sublinear_tf=True)
    tfidf.fit([text[g] for g in universe])
    cache = E.EmbCache(os.path.join(ROOT, "data", "cache", "emb.sqlite"))
    enc = E.SentenceEncoder("minilm")

    out = []
    for r in elig:
        d = design[r["term"]]
        tr, ytr, te, yte, n_dropped = split_sets(d)
        mask = frozenset(r["mask"])
        feats = {}
        for name, genes in (("tr", tr), ("te", te)):
            full = [text[g] for g in genes]
            feats[("minilm:full", name)] = (cache.embed(enc.name, full, enc), True)
            feats[("minilm:target", name)] = (cache.embed(enc.name, [D.apply_mask(t, mask)[0] for t in full], enc), True)
            feats[("ada:full", name)] = (np.stack([ada[g] for g in genes]), True)
            feats[("gene2vec", name)] = (np.stack([g2v[g] for g in genes]), True)
            feats[("tfidf:full", name)] = (tfidf.transform(full), False)
        row = {"term": r["term"], "name": r["name"], "n_train_pos": int(sum(ytr)), "n_train_neg": len(ytr) - int(sum(ytr)),
               "new_pos_removed_from_neg": n_dropped, "n_test_pos": int(sum(yte)), "n_test_neg": len(yte) - int(sum(yte))}
        for fs in FEATURES:
            Xtr, sc = feats[(fs, "tr")]
            Xte, _ = feats[(fs, "te")]
            row[fs] = fit_predict_auroc(Xtr, ytr, Xte, yte, sc)
        out.append(row)
        print(f"temporal {r['term']} done", flush=True)
    cache.close()

    cv = {}
    for r in elig:
        p = os.path.join(ROOT, "results", "main", r["term"].replace(":", "_") + ".json")
        if os.path.exists(p):
            m = json.load(open(p))["metrics"]
            cv[r["term"]] = m["ada:full"]["auroc"] - m["gene2vec"]["auroc"]
    S = {"n_terms": len(out), "median_auroc": {fs: float(np.median([o[fs] for o in out])) for fs in FEATURES},
         "median_gap_ada_minus_gene2vec": {
             "cv": float(np.median([cv[o["term"]] for o in out])) if len(cv) == len(out) else None,
             "temporal": float(np.median([o["ada:full"] - o["gene2vec"] for o in out]))},
         "n_test_pos": {"min": min(o["n_test_pos"] for o in out), "median": float(np.median([o["n_test_pos"] for o in out])),
                        "max": max(o["n_test_pos"] for o in out)},
         "per_term": out}
    os.makedirs(args.out_dir, exist_ok=True)
    json.dump(S, open(os.path.join(args.out_dir, "temporal.json"), "w"), indent=1)
    with open(os.path.join(args.out_dir, "temporal_per_term.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    lines = [f"TEMPORAL HOLDOUT (exploratory): {S['n_terms']} temporal-eligible main terms",
             f"test positives per term: min={S['n_test_pos']['min']} median={S['n_test_pos']['median']} max={S['n_test_pos']['max']}",
             "| feature set | median temporal AUROC |", "|---|---|"]
    lines += [f"| {fs} | {v:.4f} |" for fs, v in S["median_auroc"].items()]
    g = S["median_gap_ada_minus_gene2vec"]
    lines.append(f"median per-term gap ada - gene2vec: CV={g['cv']:.4f} temporal={g['temporal']:.4f} (same {S['n_terms']} terms)"
                 if g["cv"] is not None else f"median per-term gap ada - gene2vec: temporal={g['temporal']:.4f} (CV results missing)")
    txt = "\n".join(lines)
    open(os.path.join(args.out_dir, "temporal.txt"), "w").write(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
