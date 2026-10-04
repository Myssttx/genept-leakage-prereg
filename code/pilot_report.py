"""Stage 3 pilot report: table, timing projection and sanity checks. Reads results/pilot/*.json and the
pilot rows of frozen/ only; never touches main terms."""
import glob
import json
import os
import re
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    frozen = os.path.join(ROOT, "frozen")
    pilot = [r for r in json.load(open(os.path.join(frozen, "terms.json"))) if r["set"] == "pilot"]
    design = json.load(open(os.path.join(frozen, "design.json")))
    order = {r["term"]: i for i, r in enumerate(pilot)}
    res = [json.load(open(p)) for p in glob.glob(os.path.join(ROOT, "results", "pilot", "*.json"))]
    assert all(r["set"] == "pilot" for r in res)
    res.sort(key=lambda r: order[r["term"]])
    print(f"pilot result files: {len(res)} of {len(pilot)}")

    cols = [("minilm full", "minilm:full"), ("minilm target", "minilm:target"), ("minilm offtarget", "minilm:offtarget"),
            ("minilm mean rand", None), ("bge full", "bge:full"), ("ada", "ada:full"), ("gene2vec", "gene2vec"),
            ("tfidf full", "tfidf:full"), ("pubmed", "pubmed"), ("random", "random")]
    print("| term | name | n_pos | target del. in pos | " + " | ".join(c for c, _ in cols) + " |")
    print("|" + "---|" * (4 + len(cols)))
    vals = {c: [] for c, _ in cols}
    for r in res:
        m = r["metrics"]
        cells = []
        for c, k in cols:
            v = (sum(m[f"minilm:{x}"]["auroc"] for x in ("rand0", "rand1", "rand2")) / 3) if k is None else m[k]["auroc"]
            vals[c].append(v)
            cells.append(f"{v:.3f}")
        print(f"| {r['term']} | {r['name']} | {r['n_pos']} | {r['deletions']['target']['tokens_deleted_pos']} | "
              + " | ".join(cells) + " |")
    print("| median | | | | " + " | ".join(f"{statistics.median(vals[c]):.3f}" for c, _ in cols) + " |")

    # timing
    log = open(os.path.join(ROOT, "handback", "stage3_pilot_run.log")).read()
    tm = re.search(r"TIMING prep=([\d.]+)s embedding=([\d.]+)s baselines=([\d.]+)s cv_wall=([\d.]+)s jobs=(\d+) "
                   r"total=([\d.]+)s terms=(\d+)", log)
    secs = [r["seconds_eval"] for r in res]
    print("\n== TIMING ==")
    print(f"per-term CV seconds (one worker, single-threaded): min={min(secs):.1f} median={statistics.median(secs):.1f} "
          f"max={max(secs):.1f} mean={statistics.mean(secs):.1f}")
    if tm:
        prep, emb, base, cvw, jobs, total, n = (float(x) for x in tm.groups())
        print(f"pilot wall: prep={prep:.1f}s embedding={emb:.1f}s baselines={base:.1f}s cv_wall={cvw:.1f}s "
              f"(jobs={int(jobs)}) total={total:.1f}s for {int(n)} terms -> {total / n:.1f}s wall per term")
        proj_serial = 150 * (statistics.mean(secs) + (emb + prep) / n) + base
        proj_par = 150 * statistics.mean(secs) / jobs + 150 * (emb + prep) / n + base
        print(f"projection for 150 terms: one worker = {proj_serial:.0f}s ({proj_serial / 3600:.2f} h); "
              f"{int(jobs)} workers = {proj_par:.0f}s ({proj_par / 3600:.2f} h)")
        print("(projection = 150 x mean per-term CV time [/ workers] + 150 x pilot embedding+prep time per term "
              "+ one-off baseline setup; the embedding part is an upper bound because cached texts are reused)")

    # sanity checks
    print("\n== SANITY CHECKS ==")
    med = {c: statistics.median(v) for c, v in vals.items()}
    checks = [
        ("random AUROC median between 0.40 and 0.60", 0.40 <= med["random"] <= 0.60, f"{med['random']:.4f}"),
        ("minilm full AUROC median > 0.60", med["minilm full"] > 0.60, f"{med['minilm full']:.4f}"),
        ("ada AUROC median > 0.60", med["ada"] > 0.60, f"{med['ada']:.4f}"),
        ("gene2vec AUROC median > 0.55", med["gene2vec"] > 0.55, f"{med['gene2vec']:.4f}"),
    ]
    bad_genes, n_genes = 0, 0
    for r in res:
        t = r["per_gene_deletions"]["target"]
        for c in ("rand0", "rand1", "rand2"):
            rr = r["per_gene_deletions"][c]
            bad_genes += sum(1 for a, b in zip(t, rr) if a < b)
        n_genes += len(t)
    checks.append(("target deleted >= each rand condition in every gene", bad_genes == 0,
                   f"violations={bad_genes} over {n_genes} genes x 3 rand conditions"))
    ov = sum(len(set(design[r["term"]]["pos"]) & set(design[r["term"]]["neg"])) for r in pilot)
    gl = sum(len(set(r["genes"][:r["n_pos"]]) & set(r["genes"][r["n_pos"]:])) for r in res)
    checks.append(("no pos/neg overlap", ov == 0 and gl == 0, f"frozen overlap={ov}; evaluated overlap={gl}"))
    checks.append(("CV splits cover every gene once per repeat", all(r["split_check"] for r in res),
                   f"{sum(r['split_check'] for r in res)}/{len(res)} terms"))
    for name, ok, detail in checks:
        print(f"{'PASS' if ok else 'FAIL'}  {name}  ({detail})")
    conv = sum(m["convergence_warnings"] for r in res for m in r["metrics"].values())
    print(f"logistic-regression convergence warnings across all pilot fits: {conv}")
    return 0 if all(ok for _, ok, _ in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
