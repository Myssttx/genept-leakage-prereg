"""Evaluate the pilot or main GO terms (preregistration Amendment 1).

Usage (from the project folder):  python code/run_study.py --set pilot|main [--jobs N] [--results-dir DIR]
--set main is allowed only when frozen/FREEZE.txt exists and is non-empty. Embeddings are computed in the
main process and stored in the sqlite cache; each term's feature matrices are then built from the cache
one term at a time (so memory stays bounded) and cross-validated in worker processes started with
'spawn'. Terms whose <results-dir>/<set>/<term>.json already exists are skipped, so <set>_summary.csv gets
one row per term.
"""
import argparse
import csv
import json
import multiprocessing as mp
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_EXPECTED = {"pilot": 10, "main": 150}


def allowed(set_name, frozen_dir):
    """(ok, reason). pilot is always allowed; main only after the public freeze is recorded."""
    if set_name == "pilot":
        return True, "pilot"
    if set_name == "main":
        p = os.path.join(frozen_dir, "FREEZE.txt")
        if os.path.exists(p) and open(p).read().strip():
            return True, "main (FREEZE.txt present)"
        return False, "--set main requires a non-empty frozen/FREEZE.txt"
    return False, f"--set {set_name} is not a valid set (pilot or main)"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--results-dir", default=os.path.join(ROOT, "results"))
    args = ap.parse_args(argv)
    frozen = os.path.join(ROOT, "frozen")
    ok, why = allowed(args.set, frozen)
    if not ok:
        print(f"REFUSED: {why}")
        return 2

    import numpy as np
    from sklearn.feature_extraction.text import TfidfVectorizer
    import design_lib as D
    import evaluate as E

    t_start = time.time()
    rows = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == args.set]
    if len(rows) != N_EXPECTED[args.set]:
        raise RuntimeError(f"expected {N_EXPECTED[args.set]} {args.set} terms, found {len(rows)}")
    design = json.load(open(os.path.join(frozen, "design.json")))
    out_dir = os.path.join(args.results_dir, args.set)
    todo = [r for r in rows if not os.path.exists(os.path.join(out_dir, r["term"].replace(":", "_") + ".json"))]
    print(f"{args.set}: {len(rows)} terms, {len(rows) - len(todo)} already done, {len(todo)} to run", flush=True)
    if not todo:
        return 0

    keys = {}
    with open(os.path.join(frozen, "universe_keys.tsv"), encoding="utf-8") as f:
        f.readline()
        for line in f:
            g, sk, ak, gk = line.rstrip("\n").split("\t")
            keys[g] = (sk, ak, gk)
    universe = sorted(keys)
    json_path, ada_path = D.pick_genept_files(os.path.join(ROOT, "data", "raw", "genept"))
    summaries = json.load(open(json_path, encoding="utf-8"))
    text = {g: summaries[keys[g][0]] for g in universe}
    obo = D.read_obo(os.path.join(ROOT, "data", "raw", "go", "baseline", "go-basic.obo"))

    def build_texts(r):
        d = design[r["term"]]
        if set(d["pos"]) & set(d["neg"]):
            raise RuntimeError(f"pos/neg overlap in {r['term']}")
        genes = list(d["pos"]) + list(d["neg"])
        phrases = [obo[r["term"]]["name"]] + list(obo[r["term"]]["exact"])
        texts, dels = E.term_texts(r, genes, text, phrases)
        return genes, [1] * len(d["pos"]) + [0] * len(d["neg"]), texts, dels

    t_prep = time.time()

    # 1) embeddings into the cache (main process); nothing is kept in memory
    cache_path = os.path.join(ROOT, "data", "cache", "emb.sqlite")
    cache = E.EmbCache(cache_path)
    t_emb = {}
    for short in E.MODELS:
        t0 = time.time()
        enc = E.SentenceEncoder(short)
        before = cache.n_encoded
        for i, r in enumerate(todo):
            _, _, texts, _ = build_texts(r)
            for c in E.TEXT_CONDS:
                cache.embed(enc.name, texts[c], enc)
            if (i + 1) % 25 == 0:
                print(f"  {short}: {i + 1}/{len(todo)} terms embedded ({cache.n_encoded - before} new texts)", flush=True)
        t_emb[short] = (time.time() - t0, cache.n_encoded - before, enc.device)
        print(f"embedded {short} on {enc.device}: {cache.n_encoded - before} new texts in {t_emb[short][0]:.1f}s", flush=True)
        del enc
    cache.close()

    # 2) non-text features
    t0 = time.time()
    with open(ada_path, "rb") as f:
        ada_all = pickle.load(f)
    need = sorted({g for r in todo for g in design[r["term"]]["pos"] + design[r["term"]]["neg"]})
    ada = {g: np.asarray(ada_all[keys[g][1]], dtype=np.float32).ravel() for g in need}
    del ada_all
    g2v_need = {keys[g][2]: g for g in need}
    g2v = {}
    with open(os.path.join(ROOT, "data", "raw", "gene2vec_dim_200_iter_9.txt"), encoding="utf-8") as f:
        for line in f:
            p = line.split()
            if p and p[0] in g2v_need:
                g2v[g2v_need[p[0]]] = np.asarray(p[1:], dtype=np.float32)
    pub, n_no_entrez = E.pubmed_counts(os.path.join(ROOT, "data", "raw", "gene2pubmed_human.tsv"),
                                       os.path.join(ROOT, "data", "raw", "hgnc_complete_set.txt"), need)
    tfidf = TfidfVectorizer(max_features=20000, min_df=2, sublinear_tf=True)
    tfidf.fit([text[g] for g in universe])
    t_base = time.time() - t0
    print(f"baselines ready in {t_base:.1f}s (tfidf vocabulary={len(tfidf.vocabulary_)}; "
          f"needed genes without HGNC entrez_id={n_no_entrez})", flush=True)

    # 3) lazily built jobs -> spawn workers
    def jobs():
        from scipy import sparse
        c = E.EmbCache(cache_path)  # opened in the thread that consumes this generator
        try:
            for r in todo:
                genes, y, texts, dels = build_texts(r)
                feats = {}
                for short, (name, _) in E.MODELS.items():
                    for cond in E.TEXT_CONDS:
                        feats[f"{short}:{cond}"] = (c.embed(name, texts[cond], _no_encoder), True)
                feats["ada:full"] = (np.stack([ada[g] for g in genes]), True)
                feats["gene2vec"] = (np.stack([g2v[g] for g in genes]), True)
                for cond in E.TFIDF_CONDS:
                    feats[f"tfidf:{cond}"] = (sparse.csr_matrix(tfidf.transform(texts[cond])), False)
                feats["pubmed"] = (np.array([[pub[g]] for g in genes]), True)
                feats["random"] = (np.stack([E.random_vector(g) for g in genes]), True)
                yield {"term": r["term"], "name": r["name"], "set": r["set"], "off_term": r["off_term"],
                       "genes": genes, "y": y, "dels": dels, "features": feats}
        finally:
            c.close()

    n_jobs = args.jobs or min(len(todo), max(1, (os.cpu_count() or 2) - 2))
    t0 = time.time()
    summary_path = os.path.join(args.results_dir, f"{args.set}_summary.csv")
    with mp.get_context("spawn").Pool(n_jobs) as pool:
        for k, res in enumerate(pool.imap_unordered(E.evaluate_term, jobs()), 1):
            res["timing"] = {"embedding_seconds_per_model": {m: v[0] for m, v in t_emb.items()},
                             "embedding_new_texts": {m: v[1] for m, v in t_emb.items()},
                             "embedding_device": {m: v[2] for m, v in t_emb.items()}}
            E.write_result(res, out_dir)
            append_summary(summary_path, res)
            print(f"done {res['term']} in {res['seconds_eval']:.1f}s ({k}/{len(todo)})", flush=True)
    t_cv = time.time() - t0
    print(f"TIMING prep={t_prep - t_start:.1f}s embedding={sum(v[0] for v in t_emb.values()):.1f}s "
          f"baselines={t_base:.1f}s cv_wall={t_cv:.1f}s jobs={n_jobs} total={time.time() - t_start:.1f}s "
          f"terms={len(todo)}", flush=True)
    print(f"DONE {args.set}")
    return 0


def _no_encoder(texts):
    raise RuntimeError(f"{len(texts)} texts missing from the embedding cache at feature-building time")


def append_summary(path, res):
    m = res["metrics"]
    row = {"term": res["term"], "name": res["name"], "set": res["set"], "n_pos": res["n_pos"], "n_neg": res["n_neg"],
           "target_del_pos": res["deletions"]["target"]["tokens_deleted_pos"],
           "target_del_all": res["deletions"]["target"]["tokens_deleted"],
           "offtarget_del_all": res["deletions"]["offtarget"]["tokens_deleted"],
           "seconds_eval": round(res["seconds_eval"], 2)}
    for k in sorted(m):
        row[f"auroc[{k}]"] = m[k]["auroc"]
    for short in ("minilm", "bge"):
        row[f"auroc[{short}:rand_mean]"] = sum(m[f"{short}:{c}"]["auroc"] for c in ("rand0", "rand1", "rand2")) / 3
    new = not os.path.exists(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        if new:
            w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    sys.exit(main())
