import os
import re
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import analyze as A  # noqa: E402
import temporal as T  # noqa: E402
import truncation as TR  # noqa: E402


# ---------------------------------------------------------------- analyze

CONDS = ["full", "target", "offtarget", "rand0", "rand1", "rand2", "phrase"]


def fake_result(i, shift, del_pos=3):
    m = {}
    for model in ("minilm", "bge"):
        base = 0.80
        vals = {"full": base, "target": base - shift, "offtarget": base - 0.001, "rand0": base - 0.002,
                "rand1": base - 0.003, "rand2": base - 0.004, "phrase": base - shift / 2}
        for c in CONDS:
            m[f"{model}:{c}"] = {"auroc": vals[c] + 0.0001 * i}
    for k, v in {"tfidf:full": 0.9, "tfidf:target": 0.85, "tfidf:offtarget": 0.89, "ada:full": 0.9,
                 "gene2vec": 0.7, "pubmed": 0.6, "random": 0.5}.items():
        m[k] = {"auroc": v}
    pos = [{"gene": f"G{i}a", "pct_minilm_full": 90.0, "target_deleted": 1},
           {"gene": f"G{i}b", "pct_minilm_full": 60.0, "target_deleted": 0}]
    return {"term": f"GO:{i:07d}", "name": "alpha signaling", "set": "main", "metrics": m, "positives": pos,
            "deletions": {"target": {"tokens_deleted_pos": del_pos}}}


def test_lrd_and_summary():
    r = fake_result(0, 0.05)
    L, R, Dd = A.lrd(r, "minilm")
    assert L == pytest.approx(0.049) and R == pytest.approx(0.047) and Dd == pytest.approx(0.05)
    s = A.summarize([0.01, 0.02, 0.03, 0.04, 0.05])
    assert s["median"] == 0.03 and s["ci_low"] <= 0.03 <= s["ci_high"]
    assert A.summarize([0.01, 0.02, 0.03, 0.04, 0.05]) == s                      # default_rng(0): reproducible
    assert np.isnan(A.wilcoxon_p([0, 0, 0]))


def test_holm():
    adj = A.holm({"H1_L": 0.01, "H2_R": 0.04})
    assert adj == {"H1_L": 0.02, "H2_R": 0.04}
    adj = A.holm({"H1_L": 0.03, "H2_R": 0.02})
    assert adj == {"H2_R": 0.04, "H1_L": 0.04}                                   # monotone step-down


def test_material_and_restricted():
    text = {f"G{i}{s}": (f"Gene Symbol G{i}{s} Involved in alpha signaling." if s == "a" else f"Gene Symbol G{i}{s} other")
            for i in range(20) for s in "ab"}
    rows = [{"term": f"GO:{i:07d}", "name": "alpha signaling"} for i in range(20)]
    res = [fake_result(i, 0.05, del_pos=(0 if i < 5 else 3)) for i in range(20)]
    S, per_term = A.analyze(res, rows, text)
    assert S["material"]["verdict"] == "MATERIAL" and len(per_term) == 20
    assert S["secondary_restricted_L_minilm"]["n"] == 15
    nc = S["exploratory"]["name_contained_minilm_full"]
    assert nc["n_contains"] == 20 and nc["n_not"] == 20 and nc["median_pct_contains"] == 90.0
    res2 = [fake_result(i, 0.0005) for i in range(20)]            # target ~ offtarget: L tiny -> not material
    S2, _ = A.analyze(res2, rows, text)
    assert S2["material"]["verdict"] == "NOT MATERIAL"
    assert "MATERIAL LEAKAGE" in A.to_text(S2)


# ---------------------------------------------------------------- temporal

def test_temporal_split_and_fit():
    d = {"pos": ["P1", "P2"], "neg": ["N1", "N2", "X1"], "new_pos": ["X1", "X2"], "tneg": ["T1"]}
    tr, ytr, te, yte, dropped = T.split_sets(d)
    assert tr == ["P1", "P2", "N1", "N2"] and ytr == [1, 1, 0, 0] and dropped == 1
    assert te == ["X1", "X2", "T1"] and yte == [1, 1, 0]
    rng = np.random.default_rng(0)
    ytr = np.array([1] * 30 + [0] * 200)
    yte = np.array([1] * 20 + [0] * 100)
    Xtr = rng.standard_normal((230, 5)) + 2 * ytr[:, None]
    Xte = rng.standard_normal((120, 5)) + 2 * yte[:, None]
    assert T.fit_predict_auroc(Xtr, ytr, Xte, yte) > 0.95


# ---------------------------------------------------------------- truncation

class WsTokenizer:
    """Fake tokenizer: one piece per [A-Za-z0-9]+ run; no special tokens."""
    def __call__(self, text, add_special_tokens=False, return_offsets_mapping=False):
        offs = [(m.start(), m.end()) for m in re.finditer(r"[A-Za-z0-9]+", text)]
        return {"input_ids": list(range(len(offs))), "offset_mapping": offs}


def test_truncation_positions():
    text = "Gene Symbol A1 one two kinase three kinases four"
    starts = TR.deletion_char_starts(text, frozenset({"kinase"}))
    assert [text[s:s + 6] for s in starts] == ["kinase", "kinase"]
    tok = WsTokenizer()
    offs = tok(text)["offset_mapping"]
    assert TR.piece_index(offs, starts[0]) == 5 and TR.piece_index(offs, starts[1]) == 7
    assert TR.beyond_count(text, frozenset({"kinase"}), tok, keep=6) == (2, 1)    # piece 7 is cut, piece 5 kept
    assert TR.beyond_count(text, frozenset({"kinase"}), tok, keep=8) == (2, 0)
    assert TR.beyond_count(text, frozenset({"zzz"}), tok, keep=1) == (0, 0)
