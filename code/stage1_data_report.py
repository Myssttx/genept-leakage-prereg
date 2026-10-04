"""Stage 1 descriptive data report. Computes NO model outputs and NO AUROC.

Usage (from the project folder):
    python code/stage1_data_report.py            # real data under data/raw
    python code/stage1_data_report.py --root X   # any folder with the same layout (tests use this)

Expected layout under <root>/data/raw/:
    genept/                       files extracted from GenePT_emebdding_v2.zip (+ zip_listing.txt)
    go/current/goa_human.gaf.gz   go/current/go-basic.obo
    go/baseline/goa_human.gaf.gz  go/baseline/go-basic.obo
    hgnc_complete_set.txt
    Homo_sapiens.gene_info.gz
    gene2pubmed_human.tsv         (tax_id 9606 rows only, header kept)
    gene2vec_dim_200_iter_9.txt

Every section prints a header line "== SECTION ==". A section that cannot run
prints "SECTION FAILED: <reason>" and the script continues.
"""
import argparse
import gzip
import json
import os
import pickle
import random
import re
import sys
from collections import Counter

EXPERIMENTAL_CODES = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}

# Sentence openings used by Alliance of Genome Resources automated descriptions
# (Kishore et al., Database 2020). Matched at the start of a sentence.
ALLIANCE_OPENINGS = [
    "predicted to enable", "predicted to be involved in", "predicted to act upstream of",
    "predicted to be located in", "predicted to be part of", "predicted to be active in",
    "predicted to contribute to", "predicted to colocalize with", "predicted to be integral component of",
    "enables ", "involved in ", "acts upstream of", "located in ", "part of ", "is active in",
    "colocalizes with", "contributes to", "implicated in", "biomarker of", "orthologous to human",
    "used to study", "is integral component of", "exhibits ",
]
SENT_SPLIT = re.compile(r"(?<=[.;])\s+")


def header(name):
    print(f"\n== {name} ==", flush=True)


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def classify_summary(text):
    """Return 'alliance_template', 'mixed' or 'curated' plus number of template sentences."""
    sents = [s.strip().lower() for s in SENT_SPLIT.split(text.strip()) if s.strip()]
    if not sents:
        return "empty", 0, 0
    hits = sum(1 for s in sents if any(s.startswith(o) for o in ALLIANCE_OPENINGS))
    if hits >= 2 and hits / len(sents) >= 0.6:
        return "alliance_template", hits, len(sents)
    if hits >= 1:
        return "mixed", hits, len(sents)
    return "curated", hits, len(sents)


def load_json_summaries(path):
    with open(path, encoding="utf-8") as f:
        obj = json.load(f)
    if isinstance(obj, dict):
        return {str(k): (v if isinstance(v, str) else json.dumps(v)) for k, v in obj.items()}
    raise ValueError(f"top-level JSON type is {type(obj).__name__}, expected dict")


def load_pickle_embeddings(path):
    import numpy as np
    with open(path, "rb") as f:
        obj = pickle.load(f)
    if not isinstance(obj, dict):
        raise ValueError(f"pickle type is {type(obj).__name__}, expected dict")
    dims = Counter()
    for v in list(obj.values())[:2000]:
        dims[len(np.asarray(v).ravel())] += 1
    return obj, dims


def read_hgnc(path):
    pc, all_sym, prev = set(), set(), {}
    with open(path, encoding="utf-8") as f:
        cols = f.readline().rstrip("\n").split("\t")
        i_sym, i_grp = cols.index("symbol"), cols.index("locus_group")
        i_prev = cols.index("prev_symbol") if "prev_symbol" in cols else None
        i_alias = cols.index("alias_symbol") if "alias_symbol" in cols else None
        for line in f:
            r = line.rstrip("\n").split("\t")
            s = r[i_sym]
            all_sym.add(s)
            if r[i_grp] == "protein-coding gene":
                pc.add(s)
            for idx in (i_prev, i_alias):
                if idx is not None and idx < len(r) and r[idx]:
                    for p in r[idx].strip('"').split("|"):
                        prev.setdefault(p, s)
    return pc, all_sym, prev


def read_gaf(path):
    """Return header lines, set of (symbol, GO id) BP pairs with experimental codes, counts."""
    hdr, pairs_exp, n_rows, n_bp, n_not = [], set(), 0, 0, 0
    with opener(path) as f:
        for line in f:
            if line.startswith("!"):
                if any(k in line for k in ("gaf-version", "date-generated", "go-version", "generated-by")):
                    hdr.append(line.strip())
                continue
            r = line.rstrip("\n").split("\t")
            if len(r) < 15:
                continue
            n_rows += 1
            if r[8] != "P":
                continue
            n_bp += 1
            if "NOT" in r[3]:
                n_not += 1
                continue
            if r[6] in EXPERIMENTAL_CODES:
                pairs_exp.add((r[2], r[4]))
    return hdr, pairs_exp, n_rows, n_bp, n_not


def read_obo_meta(path):
    meta, n_terms, n_bp, n_obsolete = [], 0, 0, 0
    in_term = False
    with opener(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(("format-version", "data-version", "date:")) and not in_term:
                meta.append(line)
            if line == "[Term]":
                in_term, n_terms = True, n_terms + 1
            elif line.startswith("namespace: biological_process"):
                n_bp += 1
            elif line == "is_obsolete: true":
                n_obsolete += 1
    return meta, n_terms, n_bp, n_obsolete


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    args = ap.parse_args()
    raw = os.path.join(args.root, "data", "raw")
    gdir = os.path.join(raw, "genept")
    rng = random.Random(0)

    summaries, embeddings = {}, {}

    header("GENEPT FILES")
    try:
        listing = os.path.join(gdir, "zip_listing.txt")
        if os.path.exists(listing):
            print(open(listing, encoding="utf-8").read().rstrip())
        for fn in sorted(os.listdir(gdir)):
            p = os.path.join(gdir, fn)
            if os.path.isdir(p):
                continue
            print(f"{fn}\t{os.path.getsize(p)} bytes")
            if fn.lower().endswith(".json"):
                summaries[fn] = load_json_summaries(p)
                print(f"  json dict entries: {len(summaries[fn])}")
            elif fn.lower().endswith((".pickle", ".pkl")):
                emb, dims = load_pickle_embeddings(p)
                embeddings[fn] = emb
                print(f"  pickle dict entries: {len(emb)}; vector dims (first 2000): {dict(dims)}")
    except Exception as e:
        print(f"GENEPT FILES FAILED: {e!r}")

    header("TEXT / EMBEDDING KEY PAIRING")
    try:
        for sf, s in summaries.items():
            for ef, e in embeddings.items():
                ks, ke = set(s), set(e)
                print(f"{sf} vs {ef}: text={len(ks)} emb={len(ke)} both={len(ks & ke)} "
                      f"text_only={len(ks - ke)} emb_only={len(ke - ks)}")
        if len(summaries) >= 2:
            names = sorted(summaries)
            a, b = summaries[names[0]], summaries[names[1]]
            common = set(a) & set(b)
            same = sum(1 for k in common if a[k] == b[k])
            longer = sum(1 for k in common if len(b[k]) > len(a[k]))
            print(f"{names[0]} vs {names[1]}: common keys={len(common)} identical text={same} "
                  f"{names[1]} longer={longer}")
    except Exception as e:
        print(f"TEXT / EMBEDDING KEY PAIRING FAILED: {e!r}")

    header("SUMMARY TEXT DESCRIPTION")
    try:
        for sf, s in summaries.items():
            lens = sorted(len(v.split()) for v in s.values())
            n = len(lens)
            q = lambda p: lens[min(n - 1, int(p * n))]
            print(f"{sf}: n={n} words min={lens[0]} q25={q(.25)} median={q(.5)} q75={q(.75)} max={lens[-1]}")
            lits = {k: sum(1 for v in s.values() if k.lower() in v.lower())
                    for k in ["Alliance of Genome Resources", "provided by RefSeq", "provided by",
                              "Predicted to", "Gene Ontology", "GO:"]}
            print(f"  literal substring counts: {lits}")
    except Exception as e:
        print(f"SUMMARY TEXT DESCRIPTION FAILED: {e!r}")

    header("ALLIANCE TEMPLATE CLASSIFICATION")
    try:
        for sf, s in summaries.items():
            cls = {}
            for k, v in s.items():
                cls[k] = classify_summary(v)[0]
            c = Counter(cls.values())
            tot = sum(c.values())
            print(f"{sf}: " + ", ".join(f"{k}={v} ({100*v/tot:.1f}%)" for k, v in sorted(c.items())))
            for label in ["alliance_template", "mixed", "curated"]:
                keys = sorted(k for k, v in cls.items() if v == label)
                for k in rng.sample(keys, min(5, len(keys))):
                    txt = s[k].replace("\n", " ")
                    print(f"  [{label}] {k}: {txt[:400]}{'...' if len(txt) > 400 else ''}")
    except Exception as e:
        print(f"ALLIANCE TEMPLATE CLASSIFICATION FAILED: {e!r}")

    header("HGNC AND ID OVERLAP")
    try:
        pc, all_sym, prev = read_hgnc(os.path.join(raw, "hgnc_complete_set.txt"))
        print(f"HGNC symbols={len(all_sym)} protein-coding={len(pc)} prev/alias entries={len(prev)}")
        sources = {f"text:{k}": set(v) for k, v in summaries.items()}
        sources.update({f"emb:{k}": set(v) for k, v in embeddings.items()})
        g2v = os.path.join(raw, "gene2vec_dim_200_iter_9.txt")
        if os.path.exists(g2v):
            with open(g2v, encoding="utf-8") as f:
                first = f.readline().split()
                syms = {first[0]} if len(first) > 3 else set()
                dim = len(first) - 1 if len(first) > 3 else None
                for line in f:
                    syms.add(line.split(" ", 1)[0].split("\t", 1)[0])
            sources["gene2vec"] = syms
            print(f"gene2vec first line tokens={len(first)} (dim={dim}), genes={len(syms)}")
        for name, ks in sources.items():
            direct = len(ks & pc)
            rescued = len({prev[k] for k in ks - all_sym if k in prev} & pc)
            print(f"{name}: keys={len(ks)} protein-coding direct={direct} ({100*direct/len(pc):.1f}% of HGNC pc) "
                  f"rescued via prev/alias={rescued} not_in_HGNC={len(ks - all_sym)}")
        inter = set(pc)
        for ks in sources.values():
            inter &= ks
        print(f"protein-coding genes present in ALL sources (direct match): {len(inter)}")
    except Exception as e:
        print(f"HGNC AND ID OVERLAP FAILED: {e!r}")

    header("GO RELEASES")
    try:
        gaf = {}
        for rel in ["baseline", "current"]:
            hdr, pairs, n_rows, n_bp, n_not = read_gaf(os.path.join(raw, "go", rel, "goa_human.gaf.gz"))
            meta, n_terms, n_bp_terms, n_obs = read_obo_meta(os.path.join(raw, "go", rel, "go-basic.obo"))
            gaf[rel] = pairs
            print(f"{rel}: GAF header {hdr}")
            print(f"  rows={n_rows} BP rows={n_bp} BP NOT-qualified={n_not} "
                  f"experimental BP (gene,term) pairs, unpropagated={len(pairs)} genes={len({g for g, _ in pairs})}")
            print(f"  OBO {meta}; terms={n_terms} BP terms={n_bp_terms} obsolete={n_obs}")
        new = gaf["current"] - gaf["baseline"]
        gone = gaf["baseline"] - gaf["current"]
        print(f"experimental BP pairs in current but not baseline (unpropagated)={len(new)}; "
              f"in baseline but not current={len(gone)}")
        per_term = Counter(t for _, t in new)
        print(f"terms gaining >=10 new genes (unpropagated): {sum(1 for v in per_term.values() if v >= 10)}")
    except Exception as e:
        print(f"GO RELEASES FAILED: {e!r}")

    header("GENE2PUBMED")
    try:
        cnt = Counter()
        with opener(os.path.join(raw, "gene2pubmed_human.tsv")) as f:
            f.readline()
            for line in f:
                r = line.split("\t")
                if r[0] == "9606":
                    cnt[r[1]] += 1
        vals = sorted(cnt.values())
        n = len(vals)
        print(f"human GeneIDs with >=1 PubMed link={n}; median papers={vals[n//2] if n else 'NA'}; max={vals[-1] if n else 'NA'}")
    except Exception as e:
        print(f"GENE2PUBMED FAILED: {e!r}")

    print("\nEND OF REPORT")


if __name__ == "__main__":
    sys.exit(main())
