"""Stage 3 runner: evaluate the 10 PILOT terms only.

Usage (from the project folder):  python code/run_study.py --set pilot [--jobs N]
Any other --set value is refused in this stage. Embedding runs in the main process (with the
sqlite cache); cross-validation runs in worker processes started with 'spawn'. Terms whose
results/pilot/<term>.json already exists are skipped, so pilot_summary.csv gets one row per term.
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
ALLOWED_SETS = {"pilot"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True)
    ap.add_argument("--jobs", type=int, default=0)
    args = ap.parse_args(argv)
    if args.set not in ALLOWED_SETS:
        print(f"REFUSED: --set {args.set} is not allowed in Stage 3 (only: pilot)")
        return 2

    import numpy as np
    from scipy import sparse
    from sklearn.feature_extraction.text import TfidfVectorizer
    import design_lib as D
    import evaluate as E

    t_start = time.time()
    frozen = os.path.join(ROOT, "frozen")
    rows = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == args.set]
    if args.set == "pilot" and len(rows) != 10:
        raise RuntimeError(f"expected 10 pilot terms, found {len(rows)}")
    assert all(r["set"] == "pilot" for r in rows)
    design = json.load(open(os.path.join(frozen, "design.json")))
    out_dir = os.path.join(ROOT, "results", args.set)
    todo = [r for r in rows if not os.path.exists(os.path.join(out_dir, r["term"].replace(":", "_") + ".json"))]
    print(f"{args.set}: {len(rows)} terms, {len(rows) - len(todo)} already done, {len(todo)} to run", flush=True)
    if not todo:
        return 0

    # universe, texts and phrases
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

    jobs_text, overlap = {}, []
    for r in todo:
        d = design[r["term"]]
        if set(d["pos"]) & set(d["neg"]):
            overlap.append(r["term"])
        genes = list(d["pos"]) + list(d["neg"])
        phrases = [obo[r["term"]]["name"]] + list(obo[r["term"]]["exact"])
        texts, dels = E.term_texts(r, genes, text, phrases)
        jobs_text[r["term"]] = (genes, [1] * len(d["pos"]) + [0] * len(d["neg"]), texts, dels)
    if overlap:
        raise RuntimeError(f"pos/neg overlap in {overlap}")
    t_prep = time.time()

    # embeddings (main process, cached)
    cache = E.EmbCache(os.path.join(ROOT, "data", "cache", "emb.sqlite"))
    emb, t_emb = {}, {}
    for short in E.MODELS:
        t0 = time.time()
        enc = E.SentenceEncoder(short)
        before = cache.n_encoded
        for term, (genes, y, texts, dels) in jobs_text.items():
            for c in E.TEXT_CONDS:
                emb[(short, term, c)] = cache.embed(enc.name, texts[c], enc)
        t_emb[short] = (time.time() - t0, cache.n_encoded - before, enc.device)
        print(f"embedded {short} on {enc.device}: {cache.n_encoded - before} new texts in {t_emb[short][0]:.1f}s", flush=True)
        del enc
    cache.close()

    # non-text features
    t0 = time.time()
    with open(ada_path, "rb") as f:
        ada_all = pickle.load(f)
    need = sorted({g for v in jobs_text.values() for g in v[0]})
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
    print(f"baselines ready in {time.time() - t0:.1f}s (tfidf vocabulary={len(tfidf.vocabulary_)}; "
          f"needed genes without HGNC entrez_id={n_no_entrez})", flush=True)
    t_base = time.time() - t0

    jobs = []
    for r in todo:
        genes, y, texts, dels = jobs_text[r["term"]]
        feats = {}
        for short in E.MODELS:
            for c in E.TEXT_CONDS:
                feats[f"{short}:{c}"] = (emb[(short, r["term"], c)], True)
        feats["ada:full"] = (np.stack([ada[g] for g in genes]), True)
        feats["gene2vec"] = (np.stack([g2v[g] for g in genes]), True)
        for c in E.TFIDF_CONDS:
            feats[f"tfidf:{c}"] = (sparse.csr_matrix(tfidf.transform(texts[c])), False)
        feats["pubmed"] = (np.array([[pub[g]] for g in genes]), True)
        feats["random"] = (np.stack([E.random_vector(g) for g in genes]), True)
        jobs.append({"term": r["term"], "name": r["name"], "set": r["set"], "off_term": r["off_term"],
                     "genes": genes, "y": y, "dels": dels, "features": feats})

    n_jobs = args.jobs or min(len(jobs), max(1, (os.cpu_count() or 2) - 2))
    t0 = time.time()
    ctx = mp.get_context("spawn")
    summary_path = os.path.join(ROOT, "results", f"{args.set}_summary.csv")
    with ctx.Pool(n_jobs) as pool:
        for res in pool.imap_unordered(E.evaluate_term, jobs):
            res["timing"] = {"embedding_seconds_per_model": {k: v[0] for k, v in t_emb.items()},
                             "embedding_new_texts": {k: v[1] for k, v in t_emb.items()},
                             "embedding_device": {k: v[2] for k, v in t_emb.items()}}
            E.write_result(res, out_dir)
            append_summary(summary_path, res)
            print(f"done {res['term']} in {res['seconds_eval']:.1f}s", flush=True)
    t_cv = time.time() - t0
    print(f"TIMING prep={t_prep - t_start:.1f}s embedding={sum(v[0] for v in t_emb.values()):.1f}s "
          f"baselines={t_base:.1f}s cv_wall={t_cv:.1f}s jobs={n_jobs} total={time.time() - t_start:.1f}s "
          f"terms={len(jobs)}", flush=True)
    print(f"DONE {args.set}")
    return 0


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
