# Build plan: GenePT leakage paper, in audited stages

Written Oct 4, 2026. Applies the fixes in `review.md` to the plan in `plan.md`. Where they disagree, **this file and `review.md` win**. `preregistration.md` overrides all of them.

## Current status (updated Oct 4, 2026)

| Stage | Status |
|---|---|
| 1. Go/no-go checks, environment and data | **Done.** See `handback/stage1.md`. GenePT Appendix B.3 does not mask words, so the idea stands. |
| 2. Design, preregistration, freeze | **Done.** Run from a self-contained spec instead of the 4-part prompt. Design in `frozen/` (17507-gene universe, 150 main + 10 pilot GO BP terms). Preregistration pushed publicly; freeze commit `67395a9` (`frozen/FREEZE.txt`). |
| 3. Pilot | **Done.** Amendment 1 (commit `446754f`, before any outcome), `code/evaluate.py`, and the pilot run on the 10 pilot terms; all sanity checks passed. The project moved out of iCloud to `~/genept-leakage`. |
| 4. Full runs + analysis | **Current.** Deviation 1 (truncation diagnostic), the 150 main terms, then `code/analyze.py`, `code/temporal.py` and `code/truncation.py`. |
| 5–7. Analysis, paper, submission | As below. |

The stage descriptions below are the original plan. Where they differ from the status table or `preregistration.md`, those win.

## How the loop works

1. You paste `STAGEn_PROMPT.md` into Claude Code on your laptop. Each prompt carries the full text of every new or changed file with SHA-256 checks, so nothing needs downloading. Once the public GitHub repo exists (Stage 2), later prompts can just say `git pull`.
2. Claude Code does that stage only, then writes `handback/stageN.md` (8 fixed sections, defined in `CLAUDE.md`) and stops.
3. You upload `handback/stageN.md` to this thread.
4. I audit it against the checks listed for that stage below, write down what changes, and write `STAGE(n+1)_PROMPT.md`, plus any code or patch it needs.

Claude Code only executes. Every design decision is made in the stage prompt or in code I hand over. Anything not covered goes under "Questions for the reviewer" and Claude Code stops.

Only Stage 1's prompt exists now. Each later prompt is written after the previous audit, because what it says depends on what the audit finds.

---

## Stage 1: Go/no-go checks, environment and data
**Goal:** know whether the paper is still worth building, and have every input file downloaded, checksummed and described, without computing any result.

Claude Code:
- Reads GenePT's leakage appendix (B.3) and the Alliance automated-summary paper, and quotes them.
- Checks free disk space, then sets up the environment.
- Downloads GenePT v2, two GO releases (one before the summary snapshot, plus the current one), the GO ontology, HGNC, Gene2vec, NCBI gene2pubmed (human rows only) and the two local embedding models.
- Runs `code/stage1_data_report.py`, a descriptive report only: gene counts and ID overlap, summary lengths, the share of summaries matching the Alliance GO-derived template, which embedding file pairs with which text, and GO release dates.

You, in parallel:
- Confirm IDASB's deadline, page limit, presentation format and review type from its page.
- Email the BIBM registration chair to ask whether a solo author gets the student rate.

**My audit checks:** B.3 does not already do term masking. Disk and files are complete, with checksums. GenePT text/embedding pairing is unambiguous. ID overlap is at least 80%. A pre-snapshot GO release exists. Alliance template share is reported. No outcome (AUROC or classifier output) was computed.

**Gate to Stage 2:**
- the venue is open (or a named fallback venue is chosen),
- B.3 leaves the masking question open,
- data checks pass.

**If B.3 already masks terms:** stop, or reframe as an independent replication with stronger controls. You decide.

## Stage 2: Pipeline code, preregistration, freeze
**Goal:** all analysis code written and unit-tested on synthetic data, and the full analysis plan frozen publicly before any real outcome is seen.

- I write most of the code after the Stage 1 audit. Claude Code installs it, runs the tests and fixes only spec-vs-code bugs.
- Code modules:
  - labels: GO propagation, experimental evidence codes, positives and negatives;
  - term sampler: 20–200 genes, Jaccard ≤ 0.5 between terms, seed 0, 150 terms plus 10 pilot terms from outside the set;
  - masker: content-word level, fixed GO-generic stoplist, three masks (target, off-target matched term, random words), phrase and parent-term ablations;
  - embedder with cache: MiniLM and bge-small, plus the released ada vectors;
  - classifier: L2 logistic regression with C = 1, standardised inputs, 5-fold stratified CV × 3 repeats;
  - baselines B1–B5;
  - temporal split from the GO release difference;
  - stats as preregistered.
- `preregistration.md` covers hypotheses, primary metric (AUROC(off-target) − AUROC(target)), threshold for "material" leakage, tests, multiple-comparison correction, exclusions and what counts as exploratory. It is committed to a **public GitHub repo** with the term list and stoplist, and the commit URL and timestamp go in the file.

**My audit checks:**
- tests cover each masker edge case (plurals, hyphens, case, overlapping synonyms);
- the term list, stoplist and seeds are frozen in the commit;
- no real-data outcome appears in any log.

**Gate:** preregistration publicly committed.

## Stage 3: Pilot on 10 non-frozen terms
**Goal:** prove the pipeline works end to end on real data without touching the frozen terms.

Sanity checks with set thresholds (exact numbers written into the Stage 3 prompt):
- random embeddings B5 score about 0.5;
- unmasked TF-IDF scores well above B5;
- after target masking, residual term-mention rate is about 0;
- the off-target mask removes a similar number of words as the target mask;
- the cache reproduces the same vectors;
- runtime per term is measured, giving a full-run time estimate.

**My audit checks:** all sanity checks pass, and the timing fits the remaining calendar.
**Gate:** if a check fails, the fix is a code bug fix logged as a deviation, never a spec change. Any spec change goes in the preregistration's Deviations section.

## Stage 4: Full runs
**Goal:** every condition × 150 frozen terms × 2 local models (+ ada unmasked) × baselines, plus the temporal holdout, with results written to `results/` as JSON along with run logs and checksums.

**My audit checks:**
- completeness (no missing term × condition cells, or each missing cell explained);
- no NaNs;
- repeat-seed spread is reasonable;
- the code commit used matches the preregistration commit, plus logged fixes.

## Stage 5: Preregistered analysis
**Goal:** run exactly the frozen analysis and produce tables and figures with numbers read from `results/` by script.

Outputs:
- primary leakage estimate with bootstrap CI and Wilcoxon p (Holm-corrected);
- RQ1 adjusted mention coefficient;
- RQ1 split by curated vs Alliance-template summaries;
- temporal-gap comparison;
- all ablations, labelled exploratory where they were not preregistered.

**My audit checks:**
- every number is traceable to a file;
- result direction is reported honestly against the threshold;
- figures read correctly.

## Stage 6: Paper
**Goal:** IEEE 2-column paper within the venue's page limit.

- I write the outline and claims list from the Stage 5 audit. Claude Code drafts the LaTeX and pulls numbers by script.
- Citations come only from the verified list in `plan.md` section 3 (plus Kishore et al. 2020 once opened). Each is re-opened before it goes in.
- The paper includes an AI-assistance disclosure per IEEE policy.

**My audit:** a full read as a reviewer; overclaiming check (no "first"); match between numbers and results files; anonymisation if review is double-blind.

## Stage 7: Submission
**Goal:** re-check the official venue page the day of submission (deadline, template, page limit, portal, presentation rule), confirm the registration plan fits the $800 budget, build the final PDF and upload. Optionally post an arXiv/bioRxiv preprint if the venue allows it.

---

## Calendar

| If IDASB's deadline is Oct 10 | If later or a fallback venue |
|---|---|
| Oct 4–5: Stage 1 | about 2 days per stage |
| Oct 6: Stage 2 (code is pre-written, so mostly tests and the GitHub freeze) | |
| Oct 7: Stages 3 and 4 | |
| Oct 8: Stage 5 | |
| Oct 9: Stage 6 | |
| Oct 10: Stage 7 | |

The Oct 10 column has no slack and competes with the TMLR study. If Stage 1 slips past Oct 5, I recommend moving to a later venue rather than compressing Stages 5–6.
