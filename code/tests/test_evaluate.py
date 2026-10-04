import hashlib
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
sys.path.insert(0, CODE)
import design_lib as D  # noqa: E402
import evaluate as E  # noqa: E402


# ---------------------------------------------------------------- text conditions

TEXT = "Gene Symbol REPAIR The alpha kinase binds beta substrates and gamma ligands of the cell."


def content_count(text):
    return len(D.content_tokens(D.split_prefix(text)[1]))


def test_rand_mask_deletes_exactly_n_content_tokens():
    out, n = E.rand_mask(TEXT, 3, "GO:1" + "REPAIR" + "0")
    assert n == 3
    assert out.startswith("Gene Symbol REPAIR ")
    assert content_count(TEXT) - content_count(out) == 3
    # stopwords and short tokens are never deleted
    body = D.split_prefix(out)[1].lower()
    for w in ("the", "and", "of", "cell"):
        assert w in body.split() or w + "." in body.split()


def test_rand_mask_seeding():
    a = E.rand_mask(TEXT, 2, "GO:1REPAIR0")[0]
    assert a == E.rand_mask(TEXT, 2, "GO:1REPAIR0")[0]            # deterministic
    outs = {E.rand_mask(TEXT, 2, f"GO:1REPAIR{k}")[0] for k in range(6)}
    assert len(outs) > 1                                          # seed number changes the draw
    assert E.rand_mask(TEXT, 0, "x") == (TEXT, 0)
    with pytest.raises(ValueError):
        E.rand_mask("Gene Symbol A1 the of", 1, "x")


def test_rand_matches_target_count_per_gene():
    row = {"term": "GO:1", "mask": ["kinase", "ligand"], "off_mask": ["zzz"]}
    texts, dels = E.term_texts(row, ["REPAIR"], {"REPAIR": TEXT}, ["alpha kinase"])
    assert dels["target"] == [2]
    for c in E.RAND_CONDS:
        assert dels[c] == [2] and texts[c][0].startswith("Gene Symbol REPAIR ")
    assert dels["offtarget"] == [0] and texts["offtarget"][0] == TEXT
    assert texts["full"][0] == TEXT and dels["phrase"] == [2]


def test_phrase_mask_exact_case_insensitive():
    t = "Gene Symbol DNA1 Involved in DNA Repair and dna repairing. DNA-repair too; dna  repair."
    out, n = E.phrase_mask(t, ["DNA repair", "dna repair"])
    assert out == "Gene Symbol DNA1 Involved in and dna repairing. DNA-repair too; ."
    assert n == 4
    # prefix is never touched, longer phrase first
    out2, n2 = E.phrase_mask("Gene Symbol REPAIR repair of DNA repair", ["repair", "DNA repair"])
    assert out2 == "Gene Symbol REPAIR of " and n2 == 3


# ---------------------------------------------------------------- embedding cache

def test_embedding_cache(tmp_path):
    calls = []

    def enc(texts):
        calls.append(list(texts))
        return np.array([[len(t), 1.0] for t in texts], dtype=np.float32)

    c = E.EmbCache(str(tmp_path / "c" / "emb.sqlite"))
    X = c.embed("m", ["aa", "b", "aa"], enc)
    assert X.shape == (3, 2) and calls == [["aa", "b"]]
    X2 = c.embed("m", ["b", "aa"], enc)
    assert len(calls) == 1 and np.allclose(X2, X[[1, 0]])        # nothing re-embedded
    c.embed("other-model", ["b"], enc)
    assert len(calls) == 2                                         # key includes the model name
    assert E.cache_key("m", "aa") == hashlib.sha256(b"maa").hexdigest()
    c.close()


# ---------------------------------------------------------------- CV and metrics

def _toy(n_pos=30, n_neg=200, seed=0):
    rng = np.random.default_rng(seed)
    y = np.array([1] * n_pos + [0] * n_neg)
    informative = rng.standard_normal((len(y), 5)) + 2.0 * y[:, None]
    noise = rng.standard_normal((len(y), 5))
    return y, informative, noise


def test_cv_scores_and_shared_splits():
    y, Xi, Xn = _toy()
    splits = E.make_splits(y)
    assert len(splits) == 15 and E.split_check(splits, len(y))
    assert all(np.array_equal(a[1], b[1]) for a, b in zip(splits, E.make_splits(y)))   # deterministic
    mi, oof = E.cv_scores(Xi, y, splits)
    mn, _ = E.cv_scores(Xn, y, splits)
    assert mi["auroc"] > 0.95 and 0.3 < mn["auroc"] < 0.7
    assert len(mi["auroc_per_repeat"]) == 3 and abs(np.mean(mi["auroc_per_repeat"]) - mi["auroc"]) < 1e-12
    assert oof.shape == (3, len(y)) and not np.isnan(oof).any()


def test_tfidf_path_has_no_scaler():
    from scipy import sparse
    y, Xi, _ = _toy()
    Xs = sparse.csr_matrix(np.abs(Xi))
    m, _ = E.cv_scores(Xs, y, E.make_splits(y), scale=False)
    assert m["auroc"] > 0.9
    with pytest.raises(Exception):
        E.cv_scores(Xs, y, E.make_splits(y), scale=True)          # StandardScaler would reject sparse centring


def test_pos_percentiles():
    y = np.array([1, 1, 0, 0, 0, 0])
    oof = np.array([[0.9, 0.2, 0.1, 0.2, 0.3, 0.4]])
    assert np.allclose(E.pos_percentiles(oof, y), [100.0, 37.5])


def test_evaluate_term_end_to_end():
    y, Xi, Xn = _toy(n_pos=25, n_neg=150)
    genes = [f"G{i}" for i in range(len(y))]
    dels = {c: [1] * len(y) for c in E.TEXT_CONDS}
    job = {"term": "GO:9", "name": "t", "set": "pilot", "off_term": "GO:8", "genes": genes, "y": y.tolist(),
           "dels": dels, "features": {"minilm:full": (Xi, True), "random": (Xn, True)}}
    res = E.evaluate_term(job)
    assert set(res["metrics"]) == {"minilm:full", "random"}
    assert res["split_check"] and res["n_pos"] == 25 and res["n_neg"] == 150
    assert len(res["positives"]) == 25 and all(0 <= p["pct_minilm_full"] <= 100 for p in res["positives"])
    assert res["deletions"]["target"]["tokens_deleted_pos"] == 25


def test_random_vector_seeded_per_gene():
    a, b = E.random_vector("TP53"), E.random_vector("TP53")
    assert a.shape == (384,) and np.array_equal(a, b) and not np.array_equal(a, E.random_vector("EGFR"))


def test_run_script_refuses_unknown_sets():
    # never pass "main" here: on the real project it is allowed and would start the main run
    for s in ("all", "bogus"):
        p = subprocess.run([sys.executable, os.path.join(CODE, "run_study.py"), "--set", s],
                           capture_output=True, text=True)
        assert p.returncode != 0 and "REFUSED" in p.stdout


def test_main_requires_freeze(tmp_path):
    import run_study as RS
    assert RS.allowed("pilot", str(tmp_path))[0]
    assert not RS.allowed("main", str(tmp_path))[0]                    # no FREEZE.txt
    (tmp_path / "FREEZE.txt").write_text("\n")
    assert not RS.allowed("main", str(tmp_path))[0]                    # empty FREEZE.txt
    (tmp_path / "FREEZE.txt").write_text("https://github.com/x/y/commit/abc 2026-10-04\n")
    assert RS.allowed("main", str(tmp_path))[0]
    assert not RS.allowed("other", str(tmp_path))[0]
