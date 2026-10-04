"""Stage 3 evaluation of one GO term (Amendment 1 of preregistration.md).

Text conditions, the embedding cache, non-text baselines, shared CV splits and the per-term
metrics. Masking uses design_lib unchanged. Model loading is lazy (only `SentenceEncoder` imports
sentence-transformers), so spawned worker processes stay light.
"""
import hashlib
import json
import os
import random
import re
import sqlite3
import time
import warnings

import numpy as np

import design_lib as D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODELS = {
    "minilm": ("sentence-transformers/all-MiniLM-L6-v2", "1110a243fdf4706b3f48f1d95db1a4f5529b4d41"),
    "bge": ("BAAI/bge-small-en-v1.5", "5c38ec7c405ec4b44b94cc5a9bb96e735b38267a"),
}
TEXT_CONDS = ["full", "target", "offtarget", "rand0", "rand1", "rand2", "phrase"]
RAND_CONDS = ["rand0", "rand1", "rand2"]
TFIDF_CONDS = ["full", "target", "offtarget"]
N_SPLITS, N_REPEATS, CV_SEED = 5, 3, 0


# ---------------------------------------------------------------- text conditions

def rand_mask(text, n, seed_str):
    """Delete n randomly chosen content tokens from the body (prefix untouched). The generator is
    random.Random(int(sha256(seed_str))); a deleted token takes one following space with it."""
    if n <= 0:
        return text, 0
    prefix, body = D.split_prefix(text)
    matches = list(D.TOKEN_WITH_SPACE.finditer(body))
    elig = [i for i, m in enumerate(matches) if D.is_content(D.norm(m.group(1)))]
    if n > len(elig):
        raise ValueError(f"asked to delete {n} tokens but only {len(elig)} content tokens exist")
    rng = random.Random(int(hashlib.sha256(seed_str.encode("utf-8")).hexdigest(), 16))
    chosen = set(rng.sample(elig, n))
    out, last = [], 0
    for i, m in enumerate(matches):
        if i in chosen:
            out.append(body[last:m.start()])
            last = m.end()
    out.append(body[last:])
    return prefix + "".join(out), n


def phrase_mask(text, phrases):
    """Delete exact, case-insensitive, full occurrences of each phrase from the body (prefix untouched).
    Longer phrases first; inner spaces match any whitespace run; a match must not touch a letter or
    digit on either side; one following space goes with it. Returns (text, tokens deleted)."""
    prefix, body = D.split_prefix(text)
    uniq = {}
    for p in phrases:
        if p.strip():
            uniq.setdefault(p.strip().lower(), p.strip())
    n = 0

    def rep(m):
        nonlocal n
        n += len(D.TOKEN.findall(m.group(0)))
        return ""

    for p in sorted(uniq.values(), key=lambda s: (-len(s), s.lower())):
        rx = re.compile(r"(?<![A-Za-z0-9])" + r"\s+".join(re.escape(w) for w in p.split()) + r"(?![A-Za-z0-9]) ?",
                        re.IGNORECASE)
        body = rx.sub(rep, body)
    return prefix + body, n


def term_texts(row, genes, text_by_gene, phrases):
    """{condition: [text per gene]} and {condition: [tokens deleted per gene]} for one term."""
    tmask, omask = frozenset(row["mask"]), frozenset(row["off_mask"])
    texts = {c: [] for c in TEXT_CONDS}
    dels = {c: [] for c in TEXT_CONDS}
    for g in genes:
        full = text_by_gene[g]
        texts["full"].append(full)
        dels["full"].append(0)
        t, n_t = D.apply_mask(full, tmask)
        texts["target"].append(t)
        dels["target"].append(n_t)
        o, n_o = D.apply_mask(full, omask)
        texts["offtarget"].append(o)
        dels["offtarget"].append(n_o)
        for k, c in enumerate(RAND_CONDS):
            r, n_r = rand_mask(full, n_t, row["term"] + g + str(k))
            texts[c].append(r)
            dels[c].append(n_r)
        p, n_p = phrase_mask(full, phrases)
        texts["phrase"].append(p)
        dels["phrase"].append(n_p)
    return texts, dels


# ---------------------------------------------------------------- embedding cache

def cache_key(model_name, text):
    return hashlib.sha256((model_name + text).encode("utf-8")).hexdigest()


class EmbCache:
    """sqlite cache of float32 embeddings keyed by sha256(model name + text)."""

    def __init__(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.con = sqlite3.connect(path)
        self.con.execute("CREATE TABLE IF NOT EXISTS emb (key TEXT PRIMARY KEY, model TEXT, dim INTEGER, vec BLOB)")
        self.con.commit()
        self.n_encoded = 0

    def _get(self, keys):
        found = {}
        keys = list(keys)
        for i in range(0, len(keys), 900):
            chunk = keys[i:i + 900]
            q = "SELECT key, vec FROM emb WHERE key IN (%s)" % ",".join("?" * len(chunk))
            for k, v in self.con.execute(q, chunk):
                found[k] = np.frombuffer(v, dtype=np.float32)
        return found

    def embed(self, model_name, texts, encoder):
        """Return an (n, d) float32 array for `texts`, encoding only texts not yet cached."""
        keys = [cache_key(model_name, t) for t in texts]
        have = self._get(set(keys))
        missing = {}
        for k, t in zip(keys, texts):
            if k not in have and k not in missing:
                missing[k] = t
        if missing:
            mk = list(missing)
            vecs = np.asarray(encoder([missing[k] for k in mk]), dtype=np.float32)
            self.con.executemany("INSERT OR REPLACE INTO emb VALUES (?,?,?,?)",
                                 [(k, model_name, int(v.shape[0]), v.tobytes()) for k, v in zip(mk, vecs)])
            self.con.commit()
            self.n_encoded += len(mk)
            have.update({k: v for k, v in zip(mk, vecs)})
        return np.stack([have[k] for k in keys])

    def close(self):
        self.con.close()


class SentenceEncoder:
    """sentence-transformers model pinned to the preregistered snapshot; batch 64, normalised, MPS if available."""

    def __init__(self, short):
        os.environ.setdefault("HF_HOME", os.path.join(ROOT, "data", "models"))
        os.environ.setdefault("HF_HUB_OFFLINE", "1")
        import torch
        from sentence_transformers import SentenceTransformer
        name, rev = MODELS[short]
        self.name = name
        self.device = "mps" if torch.backends.mps.is_available() else "cpu"
        self.model = SentenceTransformer(name, revision=rev, device=self.device)

    def __call__(self, texts):
        return self.model.encode(list(texts), batch_size=64, normalize_embeddings=True,
                                 convert_to_numpy=True, show_progress_bar=False)


# ---------------------------------------------------------------- baselines

def random_vector(gene, dim=384):
    seed = int.from_bytes(hashlib.sha256(("random-" + gene).encode("utf-8")).digest()[:8], "big")
    return np.random.default_rng(seed).standard_normal(dim)


def pubmed_counts(gene2pubmed_path, hgnc_path, genes):
    """{gene: log1p(# distinct PubMed IDs for its HGNC entrez_id)}; 0 if no id or no rows. Also returns #genes without an id."""
    with open(hgnc_path, encoding="utf-8") as f:
        idx = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        entrez = {}
        for line in f:
            r = line.rstrip("\n").split("\t")
            e = r[idx["entrez_id"]] if idx["entrez_id"] < len(r) else ""
            if e:
                entrez[r[idx["symbol"]]] = e
    want = {entrez[g] for g in genes if g in entrez}
    pmids = {}
    with open(gene2pubmed_path, encoding="utf-8") as f:
        f.readline()
        for line in f:
            r = line.rstrip("\n").split("\t")
            if len(r) >= 3 and r[0] == "9606" and r[1] in want:
                pmids.setdefault(r[1], set()).add(r[2])
    out = {g: float(np.log1p(len(pmids.get(entrez.get(g), ())))) for g in genes}
    return out, sum(1 for g in genes if g not in entrez)


# ---------------------------------------------------------------- CV and metrics

def make_splits(y):
    from sklearn.model_selection import RepeatedStratifiedKFold
    rskf = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=CV_SEED)
    return [(tr, te) for tr, te in rskf.split(np.zeros(len(y)), y)]


def cv_scores(X, y, splits, scale=True):
    """Fit StandardScaler (if scale) + LogisticRegression on each training fold; pool out-of-fold
    scores per repeat. Returns (metrics dict, oof array of shape (N_REPEATS, n))."""
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    n = len(y)
    oof = np.full((N_REPEATS, n), np.nan)
    n_conv = 0
    for i, (tr, te) in enumerate(splits):
        lr = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
        model = make_pipeline(StandardScaler(), lr) if scale else lr
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", ConvergenceWarning)
            model.fit(X[tr], y[tr])
        n_conv += sum(1 for x in w if issubclass(x.category, ConvergenceWarning))
        oof[i // N_SPLITS, te] = model.predict_proba(X[te])[:, 1]
    if np.isnan(oof).any():
        raise RuntimeError("some out-of-fold scores are missing")
    au = [float(roc_auc_score(y, oof[r])) for r in range(N_REPEATS)]
    ap = [float(average_precision_score(y, oof[r])) for r in range(N_REPEATS)]
    return {"auroc": float(np.mean(au)), "auprc": float(np.mean(ap)), "auroc_per_repeat": au,
            "auprc_per_repeat": ap, "convergence_warnings": n_conv}, oof


def pos_percentiles(oof, y):
    """For each positive: 100 * (#neg lower + 0.5 * #neg equal) / #neg per repeat, averaged over repeats."""
    y = np.asarray(y)
    out = np.zeros(int((y == 1).sum()))
    for r in range(oof.shape[0]):
        neg = np.sort(oof[r][y == 0])
        pos = oof[r][y == 1]
        lo = np.searchsorted(neg, pos, side="left")
        hi = np.searchsorted(neg, pos, side="right")
        out += 100.0 * (lo + 0.5 * (hi - lo)) / len(neg)
    return out / oof.shape[0]


def evaluate_term(job):
    """Worker entry point. `job` holds term info, genes, y, per-gene deletions and the feature matrices
    {name: (X, scale)}. Runs every feature set on the same splits. Returns the result dict."""
    from threadpoolctl import threadpool_limits
    t0 = time.time()
    y = np.asarray(job["y"])
    splits = make_splits(y)
    metrics = {}
    with threadpool_limits(1):
        for name, (X, scale) in job["features"].items():
            m, oof = cv_scores(X, y, splits, scale)
            metrics[name] = m
            if name == "minilm:full":
                pct = pos_percentiles(oof, y)
    n_pos = int(y.sum())
    genes = job["genes"]
    positives = [{"gene": genes[i], "pct_minilm_full": float(pct[i]), "target_deleted": int(job["dels"]["target"][i])}
                 for i in range(n_pos)]
    dels = {}
    for c, d in job["dels"].items():
        d = np.asarray(d)
        dels[c] = {"tokens_deleted": int(d.sum()), "tokens_deleted_pos": int(d[y == 1].sum()),
                   "tokens_deleted_neg": int(d[y == 0].sum()), "genes_changed": int((d > 0).sum())}
    return {"term": job["term"], "name": job["name"], "set": job["set"], "off_term": job["off_term"],
            "n_pos": n_pos, "n_neg": int(len(y) - n_pos), "metrics": metrics, "deletions": dels,
            "per_gene_deletions": {c: [int(x) for x in v] for c, v in job["dels"].items()},
            "genes": genes, "positives": positives, "split_check": split_check(splits, len(y)),
            "seconds_eval": time.time() - t0}


def split_check(splits, n):
    """Each repeat's test folds must cover every gene exactly once."""
    ok = True
    for r in range(N_REPEATS):
        cover = np.concatenate([te for tr, te in splits[r * N_SPLITS:(r + 1) * N_SPLITS]])
        ok &= len(cover) == n and len(set(cover.tolist())) == n
    return bool(ok)


def write_result(res, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, res["term"].replace(":", "_") + ".json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)
    return path
