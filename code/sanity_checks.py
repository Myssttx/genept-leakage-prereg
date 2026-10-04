"""Sanity checks for a finished run (Stage 4). Prints PASS/FAIL per check; exit code 1 if any fail.

Usage:  python code/sanity_checks.py --set main [--out-dir results/analysis]
"""
import argparse
import glob
import hashlib
import json
import os
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_EXPECTED = {"pilot": 10, "main": 150}


def frozen_ok(frozen):
    ok = True
    for line in open(os.path.join(frozen, "SHA256SUMS")):
        h, name = line.split()
        ok &= hashlib.sha256(open(os.path.join(frozen, name), "rb").read()).hexdigest() == h
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", default="main")
    ap.add_argument("--results-dir", default=os.path.join(ROOT, "results"))
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "results", "analysis"))
    args = ap.parse_args(argv)
    frozen = os.path.join(ROOT, "frozen")
    rows = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == args.set]
    design = json.load(open(os.path.join(frozen, "design.json")))
    res = [json.load(open(p)) for p in glob.glob(os.path.join(args.results_dir, args.set, "*.json"))]
    med = lambda k: statistics.median(r["metrics"][k]["auroc"] for r in res)  # noqa: E731
    ov_frozen = sum(len(set(design[r["term"]]["pos"]) & set(design[r["term"]]["neg"])) for r in rows)
    ov_eval = sum(len(set(r["genes"][:r["n_pos"]]) & set(r["genes"][r["n_pos"]:])) for r in res)
    checks = [
        ("random AUROC median between 0.40 and 0.60", 0.40 <= med("random") <= 0.60, f"{med('random'):.4f}"),
        ("minilm full AUROC median > 0.60", med("minilm:full") > 0.60, f"{med('minilm:full'):.4f}"),
        (f"{N_EXPECTED[args.set]} result files", len(res) == N_EXPECTED[args.set]
         and {r['term'] for r in res} == {r['term'] for r in rows}, f"{len(res)} files"),
        ("no pos/neg overlap", ov_frozen == 0 and ov_eval == 0, f"frozen={ov_frozen} evaluated={ov_eval}"),
        ("frozen/ checksums still OK", frozen_ok(frozen), "frozen/SHA256SUMS"),
    ]
    extra = [("CV splits cover every gene once per repeat", all(r["split_check"] for r in res),
              f"{sum(r['split_check'] for r in res)}/{len(res)}"),
             ("logistic-regression convergence warnings", True,
              str(sum(m["convergence_warnings"] for r in res for m in r["metrics"].values())))]
    lines = [f"SANITY CHECKS ({args.set})"]
    lines += [f"{'PASS' if ok else 'FAIL'}  {n}  ({d})" for n, ok, d in checks]
    lines += [f"INFO  {n}: {d}" for n, _, d in extra]
    txt = "\n".join(lines)
    os.makedirs(args.out_dir, exist_ok=True)
    open(os.path.join(args.out_dir, f"sanity_{args.set}.txt"), "w").write(txt + "\n")
    print(txt)
    return 0 if all(ok for _, ok, _ in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
