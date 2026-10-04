"""Stage 2: build the study design from labels and text only and freeze it to frozen/.

Usage (from the project folder):  python code/stage2_freeze.py
Refuses to run if frozen/terms.json exists. If the universe has fewer than 15000 genes or fewer than
160 terms qualify, it prints the numbers, writes nothing and exits with code 2. No model is trained
and no score is computed.
"""
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import design_lib as D  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
FROZEN = os.path.join(ROOT, "frozen")
PATHS = {
    "hgnc": os.path.join(RAW, "hgnc_complete_set.txt"),
    "genept_dir": os.path.join(RAW, "genept"),
    "gene2vec": os.path.join(RAW, "gene2vec_dim_200_iter_9.txt"),
    "base_gaf": os.path.join(RAW, "go", "baseline", "goa_human.gaf.gz"),
    "base_obo": os.path.join(RAW, "go", "baseline", "go-basic.obo"),
    "cur_gaf": os.path.join(RAW, "go", "current", "goa_human.gaf.gz"),
    "cur_obo": os.path.join(RAW, "go", "current", "go-basic.obo"),
}
N_MAIN, N_PILOT, MIN_UNIVERSE = 150, 10, 15000


def mmm(xs):
    return f"min={min(xs)} median={statistics.median(xs)} max={max(xs)}" if xs else "n/a"


def main():
    if os.path.exists(os.path.join(FROZEN, "terms.json")):
        print("REFUSED: frozen/terms.json already exists")
        return 1
    R = D.build_design(PATHS, n_main=N_MAIN, n_pilot=N_PILOT)
    S = R["stats"]
    print("== FILES ==")
    print(f"summaries: {os.path.relpath(R['files']['summaries'], ROOT)}")
    print(f"ada: {os.path.relpath(R['files']['ada'], ROOT)}")
    print("== HGNC AND SOURCES ==")
    print(f"HGNC approved={S['hgnc_approved']} protein-coding={S['hgnc_protein_coding']} "
          f"unique prev/alias={S['hgnc_unique_aliases']} ambiguous prev/alias dropped={S['hgnc_ambiguous_aliases_dropped']}")
    for name in ("summary", "ada", "gene2vec"):
        print(f"{name}: keys={S[name + '_keys']} mapping={S[name + '_map']} protein-coding={S[name + '_protein_coding']}")
    print(f"UNIVERSE: {S['universe']} genes (universe genes whose summary is only the prefix: {S['universe_empty_body']})")
    print("== GO LABELS ==")
    for rel, label in (("base", "baseline"), ("cur", "current")):
        print(f"{label}: obo data-version={S[rel + '_obo_version']} gaf={S[rel + '_gaf']}")
        print(f"{label}: universe genes with >=1 experimental BP annotation = {S[rel + '_annotated_genes']}")
    print("== TERM SELECTION ==")
    print(f"candidates (BP, not obsolete, not root, {20}-{200} positives, non-empty mask) = {S['candidates']}")
    print(f"selected = {S['selected']} (needed {N_MAIN + N_PILOT})")

    if S["universe"] < MIN_UNIVERSE or S["selected"] < N_MAIN + N_PILOT:
        print(f"STOP: universe={S['universe']} (need >= {MIN_UNIVERSE}), selected={S['selected']} "
              f"(need {N_MAIN + N_PILOT}). Nothing written to frozen/.")
        return 2

    rows, design = R["rows"], R["design"]
    main_rows = [r for r in rows if r["set"] == "main"]
    pilot_rows = [r for r in rows if r["set"] == "pilot"]
    print(f"main={len(main_rows)} pilot={len(pilot_rows)}")
    elig = [r for r in main_rows if design[r["term"]]["eligible"]]
    print(f"main terms temporal-eligible (>=10 new_pos, live in current) = {len(elig)}")
    st = {}
    for r in rows:
        k = (r["set"], design[r["term"]]["temporal_status"])
        st[k] = st.get(k, 0) + 1
    print(f"temporal status counts: {dict(sorted(st.items()))}")
    for label, rs in (("main", main_rows), ("all 160", rows)):
        print(f"{label}: target deletions {mmm([r['deletions'] for r in rs])}; "
              f"off-target deletions {mmm([r['off_deletions'] for r in rs])}")
        print(f"{label}: n_pos {mmm([r['n_pos'] for r in rs])}; n_neg {mmm([len(design[r['term']]['neg']) for r in rs])}; "
              f"new_pos {mmm([len(design[r['term']]['new_pos']) for r in rs])}; "
              f"tneg {mmm([len(design[r['term']]['tneg']) for r in rs])}")
    ov = [len(set(design[r["term"]]["new_pos"]) & set(design[r["term"]]["neg"])) for r in rows]
    print(f"new_pos genes that are also in the term's baseline negatives: total={sum(ov)} in {sum(1 for x in ov if x)} terms")
    sel_terms = {r["term"] for r in rows}
    offs = [r["off_term"] for r in rows]
    print(f"off-target terms: distinct={len(set(offs))}; also a selected term={sum(1 for o in offs if o in sel_terms)}")
    zero = [r["term"] for r in rows if r["deletions"] == 0]
    print(f"selected terms whose target mask deletes 0 tokens over annotated summaries: {len(zero)}")

    print("== FROZEN TERMS ==")
    print("set\tterm\tname\tn_pos\tmask\tdeletions\toff_term\toff_name\toff_mask\toff_deletions\tn_neg\tn_new_pos\teligible")
    for r in rows:
        d = design[r["term"]]
        print(f"{r['set']}\t{r['term']}\t{r['name']}\t{r['n_pos']}\t{' '.join(r['mask'])}\t{r['deletions']}\t"
              f"{r['off_term']}\t{r['off_name']}\t{' '.join(r['off_mask'])}\t{r['off_deletions']}\t"
              f"{len(d['neg'])}\t{len(d['new_pos'])}\t{d['eligible']}")

    print("== MASKING EXAMPLES (first 2 main terms; first positive gene, by symbol, with >=1 target deletion) ==")
    for r in main_rows[:2]:
        m = frozenset(r["mask"])
        for g in design[r["term"]]["pos"]:
            t, n = D.apply_mask(R["text"][g], m)
            if n:
                print(f"TERM {r['term']} {r['name']} | mask={r['mask']} | gene {g} | deleted {n}")
                print(f"  BEFORE: {R['text'][g]}")
                print(f"  AFTER : {t}")
                break

    sums = D.write_frozen(FROZEN, R)
    print("== frozen/SHA256SUMS ==")
    print("".join(sums).rstrip())
    print("END OF DESIGN REPORT")
    return 0


if __name__ == "__main__":
    sys.exit(main())
