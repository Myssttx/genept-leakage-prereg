# Preregistration: how much of GenePT's GO prediction accuracy depends on the summary naming the term?

Status: frozen with the design in `frozen/`. The commit URL and time are in `frozen/FREEZE.txt`.
Written 2026-10-04, before any model was trained or any score was computed on real data.

## 1. Question

GenePT embeds each gene's NCBI text summary and uses the embedding to predict gene function. Many summaries contain the words of the very Gene Ontology (GO) term being predicted. **How much of the prediction accuracy for a GO Biological Process term depends on the summary containing that term's own words?**

We measure this by deleting the term's words from every summary (the target mask) and comparing against deleting the words of a different, matched GO term (the off-target mask).

## 2. Data

The `data/` folder is not in this repository. The inputs are identified by their SHA-256.

| File | Source | SHA-256 |
|---|---|---|
| `genept/NCBI_summary_of_genes.json` | GenePT v2, Zenodo record 10833191 | `3db721b7cfbc35795428c6a7ab4e8a9c59fb5f5ccb8c64af50dd3c97c2d4de2e` |
| `genept/GenePT_gene_embedding_ada_text.pickle` | GenePT v2, Zenodo record 10833191 | `fd297510ddd3040744033fde0b0f2cf15a40ac8b2fd2fb02f10667295e55c862` |
| `gene2vec_dim_200_iter_9.txt` | Gene2vec (Du et al. 2019), GitHub | `3441d642d4b153e4366dac02d930b3d941b3002d5484a8b3b1b85974d4feeb61` |
| `hgnc_complete_set.txt` | HGNC, downloaded 2026-10-04 | `2b4224ea847df2fc6982f5b2a52804c5d92fbb8810b134afb636f6029452dc03` |
| `go/baseline/goa_human.gaf.gz` | GO release 2024-03-28 | `a561c363d3f15caf173f3cf6101988b6dbed088de645004746b2d9f3f1f3c264` |
| `go/baseline/go-basic.obo` | GO release 2024-03-28 | `036ac170c4cc38e35af1749e208b1066b8a9e42edec4cc3aa5ed8befddb1e672` |
| `go/current/goa_human.gaf.gz` | current.geneontology.org, 2026-10-04 | `db472faff1785878521693af62646546cea6af4a386609dd01c72c0554a46a30` |
| `go/current/go-basic.obo` | current.geneontology.org, 2026-10-04 (data-version releases/2026-07-26) | `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc` |

## 3. Design rules

All rules are implemented in `code/design_lib.py` and run by `code/stage2_freeze.py`.

1. **Universe.** HGNC protein-coding genes that appear in all three of: the GenePT NCBI summary JSON, the GenePT ada embedding pickle, and Gene2vec.
   - **Picking the GenePT files:**
     - The text file is the single `.json` whose name contains "ncbi" and not "uniprot".
     - The embedding file is the single pickle whose name contains "ada" and none of "uniprot", "protein", "model_3".
     - Any ambiguity raises an error.
   - **Mapping source keys to HGNC symbols:**
     - A key that is an approved HGNC symbol maps to itself.
     - Otherwise the key is looked up among HGNC previous and alias symbols. Aliases that point to more than one gene are dropped.
     - If several keys of one source map to the same gene, a key equal to the approved symbol wins. If there are only alias keys and more than one of them, the gene is dropped from that source.
2. **Labels.**
   - **Rows used:** GO Biological Process rows from `goa_human.gaf.gz` with evidence codes EXP, IDA, IPI, IMP, IGI or IEP, excluding rows whose qualifier contains NOT.
   - **Term IDs:**
     - Alternative GO IDs are mapped to their primary ID.
     - Rows pointing to obsolete or unknown terms are dropped.
     - GAF gene symbols are mapped to HGNC as in rule 1.
   - **Propagation:** labels are propagated up `is_a` and `part_of` using the same release's `go-basic.obo`.
   - **Releases:** the baseline release (2024-03-28) defines the design; the current release is used only for the temporal sets.
   - **Annotated genes:** universe genes with at least one such annotation.
3. **Text.**
   - **Prefix:** every summary starts with "Gene Symbol X ". This prefix is split off and never masked.
   - **Tokens:** runs of `[A-Za-z0-9]`, lowercased. A plural "s" is stripped if the token is longer than 3 characters and does not end in "ss", "us" or "is".
   - **Content tokens:** tokens of at least 2 characters that are not stopwords.
   - **Stoplist:** `frozen/stoplist.txt`. Stopwords are matched after the same normalisation, so a listed word and its normalised form are both stopped.
4. **Mask set** of a term: the content tokens of its name plus its EXACT synonyms.
   - Masking deletes every body token whose normalised form is in the mask set; a deleted token takes one following space with it.
   - **Deletions** of a term: the number of tokens its mask deletes, summed over the summaries of all annotated genes.
5. **Terms.**
   - **Candidates:** non-obsolete BP terms other than the root, with 20 to 200 positives among annotated universe genes and a non-empty mask set.
   - **Order:** candidates are sorted by GO ID and shuffled with `random.Random(0)`.
   - **Acceptance:** a term is accepted if its positive-gene Jaccard with every accepted term is at most 0.5 and it has an off-target match.
   - **Stop and split:** selection stops at 160 terms. The first 150 accepted are **main** and the last 10 are **pilot**.
6. **Off-target match** for term T: another candidate term that meets all of the following:
   - It is not an ancestor or descendant of T.
   - It shares no mask token with T.
   - Its positive-gene Jaccard with T is at most 0.05.
   - Its deletion count is closest to T's, within max(5, 25% of T's count).

   Ties are broken with `random.Random("off-" + T)`.
7. **Negatives** per term: annotated universe genes that are not positive for the term and not positive (after propagation) for any of its direct `is_a`/`part_of` parents. Up to 2000 are sampled with `random.Random("neg-" + term)`.
8. **Temporal sets** per term, from the current release:
   - **new_pos:** universe genes positive in current but not in baseline.
   - **Temporal negatives:** current-annotated genes that meet all of the following, sampled up to 2000 with `random.Random("tneg-" + term)`:
     - not positive in current;
     - not positive for any of the term's direct parents in the current ontology;
     - not among the term's baseline positives or sampled negatives.
   - **Eligibility:** a term is temporal-eligible if it has at least 10 new_pos and is not obsolete (or missing) in current.
9. **Outputs:**
   - `frozen/universe.txt`, plus `frozen/universe_keys.tsv`, which gives the source key used for each gene.
   - `frozen/terms.json`
   - `frozen/design.json` (pos, neg, new_pos, tneg and eligible per term)
   - `frozen/stoplist.txt`
   - `frozen/SHA256SUMS`

   The script refuses to overwrite `frozen/terms.json`.

### Frozen design (facts from `code/stage2_freeze.py`, no outcomes)
- **Universe:** 17507 genes. Universe genes with at least one experimental BP annotation: 9691 (baseline) and 10268 (current).
- **Terms:**
  - Candidates: 1933. Selected: 150 main and 10 pilot.
  - Temporal-eligible main terms: 33.
  - One main term (GO:0046683) has 0 target deletions.
- The full design report is `handback/stage2_design_report.txt`.

## 4. Analysis plan

- **Models:**
  - Primary: `sentence-transformers/all-MiniLM-L6-v2` (snapshot `1110a243fdf4706b3f48f1d95db1a4f5529b4d41`).
  - Confirmatory: `BAAI/bge-small-en-v1.5` (snapshot `5c38ec7c405ec4b44b94cc5a9bb96e735b38267a`).
- **Embedded text:** the prefix plus the (masked) body, with the full stored text otherwise unchanged.
- **Conditions per main term**, using the same genes (pos ∪ neg) and the same splits:
  - **target mask:** the term's own mask applied to every summary;
  - **off-target mask:** the matched off-target term's mask applied to every summary.
- **Classifier:**
  - L2 logistic regression with C = 1 and `class_weight="balanced"` (scikit-learn), labels 1 for pos and 0 for neg.
  - Standardisation (`StandardScaler`) is fit on the training folds only.
  - 5-fold stratified CV repeated 3 times: `RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=0)`.
- **AUROC per term and condition:** in each repeat, the out-of-fold predicted probabilities of the 5 folds are pooled and one AUROC is computed; the term's AUROC is the mean over the 3 repeats.
- **Primary hypothesis H1 (MiniLM, 150 main terms):** Δ_t = AUROC(off-target mask) − AUROC(target mask) > 0.
  - **Estimate:** the median of Δ_t over terms.
  - **Interval:** a 95% percentile bootstrap CI, from 10,000 resamples of terms with replacement, `numpy.random.default_rng(0)`.
  - **Test:** a Wilcoxon signed-rank test of Δ_t against 0, two-sided, α = 0.05, `scipy.stats.wilcoxon` defaults.
  - **Materiality:** leakage is called **material** if the CI lower bound is > 0 and the median is ≥ 0.02.
- **Confirmatory:** the same H1 analysis with bge-small. It is reported alongside the primary result and does not replace it.
- **Exclusions:**
  - None. All 150 main terms are analysed, including any term whose masks delete nothing.
  - Pilot terms are used only to check that the pipeline runs and are excluded from H1 and the confirmatory analysis.
- **Secondary:** the temporal sets are frozen here, but no confirmatory temporal test is specified; any temporal analysis will be reported as exploratory. The same goes for any analysis not listed above.
- **Software:** Python 3.12.7, scikit-learn 1.9.1, scipy 1.18.1, numpy 2.5.3, sentence-transformers 6.1.0 (`requirements.txt`).

## 5. Deviations

None yet.
