# GenePT leakage study: do gene summaries give away the function labels?

Read `preregistration.md` (the frozen study plan and its amendments), `STAGES.md` (the staged build), `review.md` (the outside review) and `plan.md` (the original plan) before doing anything. Where they disagree, `preregistration.md` wins, then `STAGES.md`, then `review.md`, then `plan.md`.

First target venue: IDASB 2026 (a workshop at IEEE BIBM 2026). This project is separate from the TMLR study; do not touch that folder.

## Current state (Oct 4, 2026)
- **Stage 1 done:** data downloaded and checksummed (`data/README.md`, `handback/stage1.md`).
- **Stage 2 done:** design built and frozen in `frozen/`. The preregistration was pushed publicly to github.com/Myssttx/genept-leakage-prereg; the freeze commit is `67395a9` (see `frozen/FREEZE.txt`). **Never edit, delete or regenerate anything in `frozen/`.**
- **Stage 3 (current):** evaluation code plus a run on the 10 **pilot** terms only (`code/run_study.py --set pilot`). Nothing may be run on the 150 main terms.
- **Stage 4 (next):** full run on the 150 main terms, then the preregistered analysis.

## Hard rules
- **Never fabricate results, numbers or citations.** Every number in the paper must come from a file in `results/`. Every citation must be opened and checked.
- **No classifier is trained and no AUROC/AUPRC is computed for the 150 main terms before Stage 4.** Stage 3 uses only the 10 pilot terms. Unit tests on synthetic data are always fine.
- Changes to the plan made before any outcome is seen go in a dated amendment in `preregistration.md`. After outcomes exist, any change to term selection, masking, the stoplist, models, metrics or hypotheses goes in its "Deviations" section with a date and a reason, and is reported in the paper.
- Do not claim "first". Say "we found no prior study that ...".

## How stages run
- **Prompts:** the user pastes each stage's prompt, which is the specification. Implement exactly what it says.
- **Decisions:** record every decision the prompt did not settle in the final report.
- **Bug fixes:** fix a bug only when the code disagrees with its written specification, and never change a specification to make a check pass.
- **Scope:** do only the stage you were given, then stop. Never start the next stage on your own.
- **Final report:** print the report the stage prompt asks for as your last message, since the user pastes terminal output into the project thread. It must stand alone, and numbers are copied verbatim, never rounded or paraphrased.
- **Before changing code:** copy `code/` to `handback/stageN_code_before/` first.
