# GenePT leakage study: do gene summaries give away the function labels?

Read `plan.md` (the original plan), `review.md` (the outside review and fixes) and `STAGES.md` (the staged build) before doing anything. Where they disagree, `STAGES.md` and `review.md` win over `plan.md`.

The author is a high-school student submitting solo, first target IDASB 2026 (a workshop at IEEE BIBM 2026). This project is separate from the TMLR study; do not touch that folder.

## Hard rules
- **Never fabricate results, numbers or citations.** Every number in the paper must come from a file in `results/`. Every citation must be opened and checked.
- **Do not train any classifier or compute any AUROC/AUPRC on real data until the preregistration is frozen** in a public GitHub repo (Stage 2). Unit tests on synthetic data are always fine.
- Do not look at per-term outcomes for the frozen term list before Stage 4. Stage 3's pilot uses only the 10 pilot terms.
- After the freeze, any change to term selection, masking, the stoplist, models, metrics or hypotheses goes in the preregistration's "Deviations" section with a date and a reason, and is reported in the paper.
- Do not claim "first". Say "we found no prior study that ...".

## Your role: execution only
Every design decision (code architecture, data handling, thresholds, what to do when a check fails) is made in the stage prompt or already in the code. Do not make new design decisions. If you hit something the prompt does not cover, do not improvise: write it under "Questions for the reviewer" in the handback with your recommended option, finish what you can, and stop. Fix a bug only when the code disagrees with its written specification. Never change a specification to fit the data.

## Staged workflow with audits
Each stage has its own prompt (`STAGEn_PROMPT.md`), written by the reviewing Claude in the project thread after it audits the previous handback. **Do only the stage you were given, then stop.** Never start the next stage on your own. Stage scopes are in `STAGES.md`.

### Handback report (required at the end of every stage)
Write `handback/stageN.md` with exactly these sections. The reviewer has no access to this computer, so the report must stand alone.
1. **Summary:** 3–5 sentences on what was done and whether the stage's goal was met.
2. **Commands run:** the exact commands, in order.
3. **Outputs:** key outputs copied verbatim (test results, counts, report text). Never round or paraphrase numbers.
4. **Code changes:** the full unified diff of every file changed in this stage (`diff -ru handback/stageN_code_before code`, or `git diff` once the repo exists). "none" if empty.
5. **Deviations:** anything done differently from the stage prompt, and why.
6. **Problems and uncertainties:** bugs found, anything odd, checks that could not be done.
7. **Questions for the reviewer:** decisions you need, each with your recommended option.
8. **Files:** new or changed files with sizes, plus SHA-256 for every data file.

Before you change any code in a stage, copy `code/` to `handback/stageN_code_before/`. At the end, tell the user exactly: "Stage N is done. Upload handback/stageN.md to the project thread for review." and stop.
