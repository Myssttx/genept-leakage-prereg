# Outside review of the IDASB plan (plan.md, Idea 1)

Written Oct 4, 2026, reading `plan.md` as a reviewer who did not write it. Nothing here is a result.

## Verdict

**Worth building, with changes.** The core question is real, small and checkable on a laptop. But the plan has one venue problem it cannot fix by itself, three design flaws that a careful BIBM reviewer would catch, and one unchecked fact that could turn the paper from "modest" into "clear finding" or could sink it.

Rating as written: borderline workshop paper. Rating with the fixes below: solid workshop paper, possibly journal-extendable (the IJCBDD special issue invite mentioned in the CFP).

---

## 1. Blocking issues (must resolve before any code)

### 1.1 The venue may already be closed
BIBM's default workshop deadline was **Sept 27**. Nobody has seen IDASB's own deadline. If it was not extended, everything below needs a different venue, and that venue's deadline, page limit and fee have not been checked. **Do not start Stage 2 until the venue is known.** Stage 1 (data download) is cheap and reusable either way.

### 1.2 GenePT's own leakage appendix is unread
GenePT's bioRxiv v2 says leakage is addressed "in detail in Appendix B.3". If B.3 already masks function words in the summaries, the core contribution is gone. Two attempts to read it from here failed (the appendix is not in the HTML or the fetched PDF text). **Stage 1 reads it first.**

### 1.3 Disk space
Project memory says the Mac had about 2 GB free on Oct 4. The GenePT zip is 574 MB and expands further; two GO annotation releases, the GO ontology, Gene2vec and two small embedding models add roughly 1 GB more. **Need about 5 GB free.** Stage 1 checks this before downloading.

---

## 2. Possible headline finding the plan missed (verify in Stage 1)

NCBI Gene shows, for many human genes, a summary written automatically by the Alliance of Genome Resources **from GO annotations themselves** ("Predicted to enable ... activity. Involved in ... Located in ..."). The method is described in Kishore et al., *Database* 2020, "Automated generation of gene summaries at the Alliance of Genome Resources" (found by search; PMC7304461 not opened, so treat as UNVERIFIED until Stage 1 opens it).

GenePT says it took "the NCBI gene database's summary section after removing hyperlinks and date information" (quoted from bioRxiv v2). **If** many GenePT summaries are these GO-derived descriptions, then for those genes the input text is partly a rewrite of the label, and leakage is by construction.

Why it matters:
- It gives the paper a concrete, checkable mechanism, not just "summaries mention terms".
- It sets up a natural stratification: **curated RefSeq summaries vs GO-derived Alliance summaries**.
- Risk: a reviewer may say "obviously leaky". The answer is that nobody has reported the share or measured the effect, and the GenePT benchmarks and their successors all evaluate on GO.

This is a hypothesis until Stage 1 counts template phrases in the actual GenePT text.

---

## 3. Design flaws (fixed in the stage plan)

### 3.1 RQ1 is confounded by how well-studied a gene is
Well-studied genes have longer summaries, mention more terms, **and** get better embeddings for other reasons. So "mention-positives are ranked higher" could be study bias, not leakage.
**Fix:** RQ1 becomes a per-term logistic regression of "positive gene ranked above the median negative" on mention flag, log PubMed count, and summary length; report the mention coefficient. Keep the raw mention vs non-mention comparison as descriptive only.

### 3.2 The random-word mask is too weak a control
Deleting random non-stopwords mostly removes filler; deleting a GO term's words removes dense biology vocabulary. Any masking of biology words will hurt more than random words.
**Fix:** add an **off-target term mask**: for each target term, mask the words of a different, randomly chosen GO term of similar name length, in the same number of genes. The leakage estimate becomes **AUROC(off-target mask) − AUROC(target mask)**. Keep the random-word mask as a secondary control.

### 3.3 Masking parent-term words would wreck the text
Parent terms include words like "regulation", "cellular", "process", "response". Masking them across all summaries is not a leakage test. **Fix:** primary masking = target term name + exact synonyms, at content-word level, with a fixed GO-generic stoplist (Stage 2 writes it before any outcome is seen). Phrase-only masking and parent-term masking become ablations.

### 3.4 Smaller fixes
- **ΔAUROC_leak simplifies.** (C0 − C1) − (C0 − C1-rand) is just C1-rand − C1. Report it that way.
- **Temporal holdout:** GAF column 14 dates can be refreshed on re-review. Use the difference between two **GO releases** (one archived release before the GenePT summary snapshot, one current) instead. Cleaner and reviewer-proof.
- **Term overlap:** GO terms share genes, so 150 terms are not 150 independent samples. Filter the sample so no two terms share more than 50% of genes (Jaccard ≤ 0.5), and report it.
- **Classifier:** inner CV for C with 20 positives is noisy and adds code. Fix C = 1.0 on standardised embeddings; ablate C ∈ {0.1, 10}.
- **Which GenePT text:** the Zenodo zip has NCBI-only and NCBI+UniProt summaries. The plan must state which text pairs with which released embedding. Stage 1 pins this.
- **One local model is thin.** Use two (all-MiniLM-L6-v2 and bge-small-en-v1.5). Optional: re-embed with a paid OpenAI model (under $1 at current list prices, unverified) so the masking result is shown on an OpenAI embedding too. Your call in Stage 2.
- **Pilot terms** must come from outside the frozen term set, so the pilot cannot influence the main analysis.
- **AI-use policy:** check BIBM/IEEE rules on disclosing AI assistance before writing.

---

## 4. What a BIBM reviewer would still attack (after the fixes)

| Attack | Answer available |
|---|---|
| "Only GO BP, only human." | Scope stated; the mask script works on any term list. Add Reactome as an ablation only if time allows. |
| "Small local models, not GenePT's ada embeddings." | Unmasked local vs ada AUROC is reported side by side; optional paid re-embed. |
| "Paraphrase leakage is not caught." | Masking is a lower bound; reported residual mention rate. |
| "This is a negative/sanity paper, not a method." | IDASB lists "data science approaches" and "LLMs in biomedical fields"; the deliverable is a corrected benchmark plus released masked text. |
| "Alliance summaries obviously leak; trivial." | The share and the effect size were not reported anywhere we found; the benchmarks people cite still score on GO. |

## 5. Feasibility for a solo high-school author
- Compute: fine on a laptop.
- Time: tight for an Oct 10 deadline while the TMLR study runs (TMLR target Oct 30). Realistic only if Stages 1-4 each take about a day. If IDASB's deadline is Oct 10 or earlier and Stage 1 slips, drop to a later venue rather than rushing Stage 5.
- Overlap: no overlap with the TMLR study (different domain, data and question). Can run in parallel without the TMLR "parallel submission" problem.
