"""Truncation diagnostic (Deviation 1). Diagnostic only; no hypothesis, analysis or rule depends on it.

Usage:  python code/truncation.py [--out-dir results/analysis]
For MiniLM and bge, using each model's own tokenizer and max_seq_length (which counts the special
tokens the model adds):
  1. % of universe summaries (prefix + body) whose tokenized length exceeds max_seq_length;
  2. for the 150 main terms, % of target-mask deletions in positives whose first word piece falls at or
     beyond the truncation point, i.e. at content position >= max_seq_length - (number of special tokens).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design_lib as D  # noqa: E402
import evaluate as E  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def deletion_char_starts(text, mask):
    """Character offsets (in the full text) of the tokens the target mask deletes."""
    prefix, body = D.split_prefix(text)
    return [len(prefix) + m.start(1) for m in D.TOKEN_WITH_SPACE.finditer(body) if D.norm(m.group(1)) in mask]


def piece_index(offsets, char_start):
    """Index of the first word piece that starts at or after char_start (skipping empty pieces)."""
    for i, (s, e) in enumerate(offsets):
        if e > s and s >= char_start:
            return i
    return len(offsets)


def beyond_count(text, mask, tokenizer, keep):
    """(n deletions, n at or beyond `keep` content pieces) for one text."""
    starts = deletion_char_starts(text, mask)
    if not starts:
        return 0, 0
    offs = tokenizer(text, add_special_tokens=False, return_offsets_mapping=True)["offset_mapping"]
    return len(starts), sum(1 for s in starts if piece_index(offs, s) >= keep)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "results", "analysis"))
    args = ap.parse_args(argv)
    frozen = os.path.join(ROOT, "frozen")
    rows = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == "main"]
    design = json.load(open(os.path.join(frozen, "design.json")))
    keys = {}
    with open(os.path.join(frozen, "universe_keys.tsv"), encoding="utf-8") as f:
        f.readline()
        for line in f:
            g, sk, _, _ = line.rstrip("\n").split("\t")
            keys[g] = sk
    jp, _ = D.pick_genept_files(os.path.join(ROOT, "data", "raw", "genept"))
    summ = json.load(open(jp, encoding="utf-8"))
    universe = sorted(keys)
    text = {g: summ[keys[g]] for g in universe}

    S = {}
    for short in E.MODELS:
        enc = E.SentenceEncoder(short)
        tok, limit = enc.model.tokenizer, enc.model.max_seq_length
        n_special = len(tok("a")["input_ids"]) - len(tok("a", add_special_tokens=False)["input_ids"])
        keep = limit - n_special
        lens = [len(ids) for ids in tok([text[g] for g in universe], add_special_tokens=True, truncation=False,
                                        verbose=False)["input_ids"]]
        n_long = sum(1 for n in lens if n > limit)
        tot = bey = genes_hit = 0
        for r in rows:
            mask = frozenset(r["mask"])
            for g in design[r["term"]]["pos"]:
                n, b = beyond_count(text[g], mask, tok, keep)
                tot += n
                bey += b
                genes_hit += 1 if b else 0
        S[short] = {"model": enc.name, "max_seq_length": limit, "special_tokens": n_special,
                    "universe_summaries": len(universe), "summaries_longer_than_limit": n_long,
                    "pct_summaries_longer_than_limit": 100.0 * n_long / len(universe),
                    "main_target_deletions_in_positives": tot, "deletions_at_or_beyond_limit": bey,
                    "pct_deletions_at_or_beyond_limit": (100.0 * bey / tot) if tot else None,
                    "positive_term_pairs_with_any_deletion_beyond": genes_hit}
        del enc
    os.makedirs(args.out_dir, exist_ok=True)
    json.dump(S, open(os.path.join(args.out_dir, "truncation.json"), "w"), indent=1)
    lines = ["TRUNCATION DIAGNOSTIC (diagnostic only)", "| model | max_seq_length | summaries > limit | % | "
             "target deletions in positives (150 main terms) | at/beyond limit | % |", "|---|---|---|---|---|---|---|"]
    for short, s in S.items():
        pct = "n/a" if s["pct_deletions_at_or_beyond_limit"] is None else f"{s['pct_deletions_at_or_beyond_limit']:.2f}"
        lines.append(f"| {short} | {s['max_seq_length']} | {s['summaries_longer_than_limit']} of {s['universe_summaries']} | "
                     f"{s['pct_summaries_longer_than_limit']:.2f} | {s['main_target_deletions_in_positives']} | "
                     f"{s['deletions_at_or_beyond_limit']} | {pct} |")
    txt = "\n".join(lines)
    open(os.path.join(args.out_dir, "truncation.txt"), "w").write(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
