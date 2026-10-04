# IDASB 2026 (IEEE BIBM workshop): venue check, candidate ideas, and full plan

Written October 4, 2026. **No experiments have been run; nothing below is a result.** Citations marked "verified" were opened on a publisher, PubMed/PMC, bioRxiv or arXiv page on Oct 4, 2026.

> **Update Oct 4, 2026:** an outside review (`review.md`) changed parts of section 2. The build now follows `STAGES.md`; where they disagree, `STAGES.md` and `review.md` win.

---

## 0. Venue facts (read this first)

What I could confirm on official BIBM 2026 pages on Oct 4:

| Item | Status | Source |
|---|---|---|
| Host conference | IEEE BIBM 2026, Dec 1-4, 2026, Dallas, TX | [BIBM 2026 CFP](https://www3.cs.stonybrook.edu/~bibm2026/callforpaper2026.html) |
| Default workshop paper deadline | **Sept 27, 2026 (already passed)** | [BIBM call for workshops](https://www3.cs.stonybrook.edu/~bibm2026/CallforWorkshopProposals.html) |
| Workshop notification / camera-ready | Oct 18 / Nov 8, 2026 | same |
| Extensions | Individual workshops set their own. Example: HP4MoDa extended to **Oct 10**; M³AI-Bio still lists Sept 27 | [HP4MoDa](https://philippe-fournier-viger.com/HP4MODA_2026/submission.html), [M³AI-Bio](https://sites.google.com/view/m3ai-bio-2026/home) |
| Page limit | Other BIBM 2026 workshops: 8 pages full / 4 short, IEEE 2-column, references included; M³AI-Bio lists $100 per extra page | same two workshop pages |
| Proceedings | BIBM workshop proceedings, IEEE Computer Society Press | BIBM call for workshops |
| Author registration | IEEE student **$560**, student non-member **$675**, IEEE member $800, non-member $960 (early, by Nov 8). Extra paper fees apply. Author registration non-refundable. | [BIBM 2026 registration](https://www3.cs.stonybrook.edu/~bibm2026/BIBM2026-registration0.html) |
| In-person required? | **Not confirmed for IDASB.** BIBM lets each workshop choose onsite, online or hybrid. HP4MoDa is online; M³AI-Bio is hybrid. | BIBM call for workshops, workshop pages |

**Not verified (the IDASB site was unreachable from my environment):** IDASB's own deadline, whether it was extended, its page limit, and its presentation format. The IDASB 2024 site was at `combio-lezhang.online/IDASB2024/`, so the 2026 page is probably `combio-lezhang.online/IDASB2026/`. Please check that page and the BIBM submission system (wi-lab.com cyberchair, workshop "ws_submit") for the IDASB entry.

**Budget:**
- Registration fits under $800 **only if you qualify for a student rate** ($560 with IEEE student membership, which costs extra; $675 without). Whether a solo author qualifies for BIBM's student rate is **not stated**; email the BIBM registration chair before paying. At the non-member rate ($960) it is **over budget**.
- If IDASB requires in-person presentation, add travel to Dallas, Dec 1-4. The related BIBM UGHS symposium requires in-person presentation ([UGHS](https://bibm2026-hs.github.io/)); its deadline (Sept 7) has passed.

**Schedule reality:** if IDASB did not extend past Sept 27, this venue is closed. If it extended to about Oct 10, you have about six days, overlapping the AAIML Oct 10 and TMLR Oct 30 work. Only Idea 1 below is small enough for that window.

**Overlap check:** none of the ideas below overlap the TMLR study (attribution vs sensor-fault damage in time-series models) or Paper 1 (C-MAPSS RUL). Different domain, data, models and question.

---

## 1. Candidate ideas, ranked

### Idea 1 (recommended): leakage audit of text-derived gene embeddings
GenePT and its successors embed each gene's NCBI text summary and use the embedding to predict gene function (GO terms, pathways). Those summaries often literally name the function being predicted. Zhong et al. (2025, v1) wrote that GenePT's strong performance "could be heavily driven by data leakage". Their v2 added a temporal GO holdout (annotations after March 2024), so **a date-based holdout alone is not new**. What no paper we found does is **mask the predicted term's words in the summaries** and compare against a random-mask control. That is this paper's core. **Open risk:** GenePT's bioRxiv v2 says leakage is addressed "in detail in Appendix B.3", which we could not open; read it first, before committing (day 1 of the work plan). Small, laptop-scale, public data, fits the CFP's "LLMs in biological and biomedical fields" and "data science approaches" topics, and is publishable whichever way it comes out. Full plan in section 2.

### Idea 2: study-bias baseline for gene-embedding benchmarks
Question: how much of a gene embedding's accuracy on GO/disease-gene tasks is matched by a trivial predictor using only how well-studied a gene is (PubMed paper count from NCBI gene2pubmed)? Study bias and multifunctionality bias are known for interaction networks (Gillis & Pavlidis 2011, 2012), but recent gene-embedding benchmarks do not report a popularity-only baseline as far as I found. Feasible on a laptop. **Weaker than Idea 1** because the effect is partly known and the conclusion is narrower. It is better used as a control inside Idea 1 (it is included there as baseline B4).

### Idea 3: stability of co-expression hub genes under sample subsampling
Question: how many samples does a public GEO/TCGA co-expression network need before its top hub genes stabilise? Feasible, but close to existing WGCNA robustness work and to network-stability studies, so novelty would be incremental at best, and the literature review alone would not fit a six-day window. **Not recommended for this deadline.**

---

## 2. Full plan for Idea 1

### Working title
**How Much Do Gene Summaries Give Away? A Leakage Audit of Text-Derived Gene Embeddings for Function Prediction**

### Target conference
IDASB 2026 at IEEE BIBM 2026 (if still open). Fallback: an open bioinformatics venue with a later deadline, or a 2027 BIBM workshop. Do not submit it to two archival venues at once.

### Research question
When text-derived gene embeddings (GenePT-style: embed the NCBI gene summary) predict whether a gene belongs to a GO Biological Process term, how much of the accuracy comes from the summary explicitly naming that term?

- **RQ1 (lexical stratification).** Is accuracy higher for positive genes whose summary mentions the term's name or synonyms than for positives whose summary does not?
- **RQ2 (term masking).** How much does accuracy drop when the term's words are removed from the summaries, compared with removing the same number of random words?
- **RQ3 (secondary, replication).** Using GAF annotation dates, does the text-vs-Gene2vec gap seen in random CV shrink on annotations added after the summaries were frozen? Zhong et al. (2025, v2) already ran a temporal holdout and found GenePT trailing PPI embeddings there; we report this only as a check that ties RQ1/RQ2 to their finding, not as a contribution.

### Exact novelty claim
A controlled measurement of how much GenePT-style function prediction depends on the summary naming the predicted term: lexical stratification (RQ1) and term masking with a random-mask control (RQ2), with a released masking script and frozen term list. **Not claimed:** first to raise leakage (GenePT's own authors and Zhong et al. raised it), temporal holdout (Zhong et al. v2; CAFA; Tan et al. 2026 for their own descriptions), or a new embedding method. Do not use the word "first" anywhere; the novelty search covered bioRxiv/arXiv/journals but cannot be exhaustive.

Contribution type: **differentiated but modest** — a measurement and benchmark correction, not a new model. Honest framing matters for a workshop reviewer.

### Closest prior work
See section 3 (verified list). Summary of why it does not already answer the question:
- GenePT and its follow-ups (GenePT Revisited, GEbench, "Understanding the LLM-based gene embeddings") report accuracy with random cross-validation splits and do not mask terms or hold out by date.
- Zhong et al. (2025) v2 runs a temporal GO holdout including GenePT, which shows generalisation drops but does not say *why*; it does not separate explicit term mentions from other knowledge.
- Tan et al. (2026, Nature Communications) test "information overlap between curated annotations and gene descriptions" for their own RAG-generated descriptions via newly added MSigDB terms, not for GenePT/NCBI summaries and not by masking.
- GenePT (Chen & Zou) says it addresses leakage in Appendix B.3; content unread (must check).
- CAFA established temporal holdout for protein function prediction but does not study text-summary embeddings.
- Wadi et al. (2016) and Tomczak et al. (2018) study annotation drift, not leakage into embeddings.

### Data (all public, free)
1. **GenePT v2** (Zenodo record 10833191, CC-BY 4.0, 574 MB zip): NCBI gene summaries for human genes plus the original text-embedding-ada-002 embeddings. Summary snapshot date is not stated on the record; record v2 published March 18, 2024, so treat the summaries as frozen **no later than March 2024** and use a cutoff after that.
2. **GO annotations**: current `goa_human.gaf` from geneontology.org. GAF 2.2 column 14 is "Date on which the annotation was made" (verified in GO docs). Caveat: the date may be refreshed on re-review, so it approximates first appearance; say so in the paper. Plus `go-basic.obo` for term names, synonyms and the hierarchy.
3. **Gene2vec** (Du et al. 2019) co-expression embeddings, public on GitHub, as the non-text control.
4. **NCBI gene2pubmed** for the study-bias baseline.
5. **Local embedder** (free, CPU): `sentence-transformers` with a small model such as `all-MiniLM-L6-v2` or `BAAI/bge-small-en-v1.5`, so masked and unmasked text are embedded by the same model.

### Experimental protocol (freeze before running anything final)
1. **Gene universe.** Human protein-coding genes present in GenePT, Gene2vec and the GAF. Record counts.
2. **Term selection (frozen by seed 0).** GO BP terms with 20-200 genes annotated by experimental evidence codes (EXP, IDA, IPI, IMP, IGI, IEP), with true-path propagation. Randomly sample 150 terms; write the list to `terms_frozen.txt` before any model is trained.
3. **Positives / negatives.** Positives: genes annotated to the term (propagated). Negatives: genes with at least one experimental BP annotation but not to this term or its descendants or ancestors.
4. **Lexical mention flag** per (gene, term): case-insensitive match of the term name, exact synonyms and a stemmed version against the gene's summary.
5. **Conditions** (all embedded with the same local model unless noted):
   - C0 full summary (local model) and C0-ada (released GenePT embeddings).
   - C1 term-masked: delete term name, exact/related synonyms and parent term names from every summary; re-embed only genes whose text changed (cache by hash).
   - C1-rand: in the same genes, delete the same number of random non-stopword tokens (5 seeds).
   - C2 temporal: cutoff date T = 2024-06-01 (fixed now). Train on annotations dated before T; test positives are genes whose first annotation to the term is after T; test negatives drawn from genes never annotated to it. Run for C0, C0-ada and Gene2vec.
6. **Classifier.** L2 logistic regression per term, class-balanced, C chosen by inner 3-fold CV on training data only. 5-fold stratified CV repeated 3 times for C0/C1/C1-rand.

### Baselines
- B1 Gene2vec (co-expression, no text).
- B2 TF-IDF bag-of-words on the same summary (does plain word matching already explain it?).
- B3 gene-symbol-only embedding (model knows the gene name only).
- B4 study-bias predictor: log PubMed count from gene2pubmed only.
- B5 random Gaussian embeddings of the same dimension (floor).

### Metrics
Per-term AUROC and AUPRC (AUPRC reported against the positive rate). Main quantity: **ΔAUROC_leak = (C0 − C1) − (C0 − C1-rand)**, i.e. drop from term masking beyond drop from random masking. RQ1: AUROC on mention-positives vs non-mention-positives (same negatives). RQ3: text-minus-Gene2vec AUROC gap in CV vs in temporal holdout.

### Ablations
- Mask only exact term name vs name + synonyms vs + parent terms.
- Local model vs released ada embeddings (C0 only).
- Term size bins (20-50, 50-100, 100-200).
- Evidence codes: experimental only vs all including IEA.

### Statistical tests (fixed now)
- Unit of analysis: GO term (n ≈ 150).
- Paired Wilcoxon signed-rank tests across terms for C0 vs C1, C1 vs C1-rand; bootstrap 95% CIs (10,000 resamples over terms) on median ΔAUROC_leak.
- RQ1: per-term difference (mention vs non-mention) only for terms with ≥5 positives in each group; Wilcoxon plus bootstrap CI.
- Holm correction across the three RQs.
- Pre-stated threshold: call leakage "material" if the 95% CI of median ΔAUROC_leak excludes 0 and the median is ≥ 0.02.

### Expected reviewer attacks and responses
1. *"Masking makes the text ungrammatical, so any drop is from damage, not leakage."* → C1-rand removes the same number of words from the same genes; ΔAUROC_leak subtracts that.
2. *"Paraphrases still leak (e.g. 'programmed cell death' for apoptosis)."* → True; masking gives a lower bound on leakage. We mask synonyms and parents, report residual mention rate, and say so.
3. *"Knowing a function from the literature before GO records it isn't leakage."* → Agreed; the temporal test is framed as prospective accuracy, not leakage. The leakage claim rests on RQ1/RQ2.
4. *"You used a small local model, not GenePT's model."* → Released ada embeddings are used for C0 and the temporal test; masking needs re-embedding, which only the local model allows for free. Optional: re-embed with a paid API for a few dollars if budget allows.
5. *"Zhong et al. already did a temporal holdout."* → Yes, and we cite it as the prior result; our contribution is the masking and mention analysis that explains where text embeddings' advantage comes from.
5b. *"GenePT's Appendix B.3 already handled leakage."* → Must read before submission; if B.3 already masks terms, this idea loses its core and should be dropped or reframed as a replication with controls.
6. *"Only GO BP, only human."* → Scope limitation stated; masking script is general.
7. *"Logistic regression is too simple."* → It matches how GenePT and the benchmarks evaluate embeddings; the question is about the embedding, not the classifier.

### Compute / time
Laptop CPU. Embedding ~20k short summaries with a MiniLM-class model is minutes to tens of minutes; masked re-embedding only touches genes whose text changes. 150 terms × 6 conditions × 15 fits of logistic regression is minutes. Total under a few hours of compute. Main cost is your time.

### Failure conditions
- Fewer than ~30 terms have ≥10 post-cutoff positives → drop RQ3 to a short exploratory section; RQ1/RQ2 still stand.
- Gene ID mismatch between GenePT (symbols) and GAF (UniProt / symbols) loses >20% of genes → map via HGNC; report losses.
- GenePT summaries turn out to mention almost no term names → RQ1 has no power; the result itself ("summaries rarely name the term, leakage is small") is still reportable.

### Minimum publishable result
RQ1 + RQ2 across ~150 frozen GO terms with the random-mask control and baselines B1-B5, whichever direction the effect goes. "Leakage is small" and "leakage is large" are both useful to the people using these embeddings.

### Strong-result scenario
Material ΔAUROC_leak and a clear mention vs non-mention gap that together explain the generalisation drop Zhong et al. v2 observed, plus a released masked benchmark and script others can reuse.

### Work plan (assumes IDASB deadline ≈ Oct 10; shift if different)
| Date | Task |
|---|---|
| Oct 4 | Confirm IDASB deadline, page limit, format and student-rate eligibility. If closed, stop here. Read GenePT bioRxiv v2 Appendix B.3 and Zhong et al. v2; if either already masks terms, stop or reframe. |
| Oct 5 | Download GenePT v2, GAF, OBO, Gene2vec, gene2pubmed; ID mapping; freeze term list and cutoff; write pre-registration note. |
| Oct 6 | Lexical flags, masking script with unit tests; embed C0, C1, C1-rand. |
| Oct 7 | Run all CV conditions and baselines; temporal split. |
| Oct 8 | Stats, figures (ΔAUROC per term; mention vs non-mention). |
| Oct 9 | Write paper (IEEE 2-column, ≤8 pages incl. refs); related work from section 3 only. |
| Oct 10 | Proofread, check official page again, submit. |

---

## 3. Closest prior work (verification status)

Checked Oct 4, 2026 by opening the publisher, bioRxiv, arXiv, PLOS or Crossref record. PMC/PubMed pages were blocked by CAPTCHA, so PMC ids are not used as sources.

| # | Paper | Status | Why it matters here |
|---|---|---|---|
| 1 | Chen & Zou. "GenePT: A Simple But Effective Foundation Model for Genes and Cells Built From ChatGPT." bioRxiv 2023, doi:10.1101/2023.10.16.562533 | Verified. **Appendix B.3 (their leakage check) NOT read** | The method being audited; authors say they address summary leakage in B.3 |
| 2 | Chen & Zou. "Simple and effective embedding model for single-cell biology built from ChatGPT." *Nature Biomedical Engineering* 9:483–493 (2025), doi:10.1038/s41551-024-01284-6 | Verified (main text; supplement not checked) | Journal version of GenePT |
| 3 | Zhong, Li, Dannenfelser, Yao. "Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks." bioRxiv 2025, doi:10.1101/2025.01.29.635607 (v1 and v2) | Verified; no journal version found | Closest work. v1 flags GenePT leakage; v2 temporal holdout after March 2024, mean AUROC 0.79→0.62, GenePT-Model3 0.774 vs Mashup 0.777 |
| 4 | Hedley, Torr, Märtens. "GenePT Revisited: Do Better Text Embeddings Make Better Gene Embeddings?" bioRxiv 2026, doi:10.64898/2026.04.16.718976 | Verified (ICLR 2026 link seen in search only, unverified) | Newer backbones +1–17%; no leakage controls |
| 5 | Huang, Hou, Zhao, et al. "Guidance for high-quality functional gene embeddings from large language models." bioRxiv 2026, doi:10.64898/2026.04.30.721875 | Verified | Quality depends on input text containing explicit functional information, which supports the leakage hypothesis |
| 6 | Cai, Gan, Zhang, Li. "Understanding the LLM-based gene embeddings." bioRxiv 2025, doi:10.64898/2025.12.19.695582 | Verified | Description embeddings recover >93% of Hallmark/C2 pathways; symbol-only >64% |
| 7 | Gan & Li. "Small, Open-Source Text-Embedding Models as Substitutes to OpenAI Models for Gene Analysis." bioRxiv 2025, doi:10.1101/2025.02.15.638462 | Verified preprint; journal version unverified | Justifies using a small local embedder |
| 8 | Tan, Wang, Liang, et al. "An embedding-based framework enables statistical testing of gene-set function hypotheses inferred by large language models." *Nature Communications* 17:9248 (2026), doi:10.1038/s41467-026-75972-z | Verified | Tests annotation/description overlap via newly added MSigDB terms, for their own descriptions |
| 9 | Brechtmann, Bechtler, Londhe, Mertes, Gagneur. "Evaluation of input data modality choices on functional gene embeddings." *NAR Genomics and Bioinformatics* 2023, doi:10.1093/nargab/lqad095 | Verified | Literature-based embeddings win on curated lists through study bias, not on GWAS signal |
| 10 | Jararweh, Macaulay, Arredondo, et al. "LitGene: a transformer-based model that uses contrastive learning to integrate textual information into gene representations." bioRxiv 2024, doi:10.1101/2024.08.07.606674 | Verified | Text gene model trained with GO signal (circularity risk) |
| 11 | Kan-Tor, Danziger, Zohar, Ninio, Shimoni. "Does your model understand genes? A benchmark of gene properties for biological and text models." arXiv:2412.04075 (2024) | Verified v1 | Text models strong on gene-property tasks |
| 12 | Du, Jia, Dai, Tao, Zhao, Zhi. "Gene2vec: distributed representation of genes based on co-expression." *BMC Genomics* 20:82 (2019), doi:10.1186/s12864-018-5370-x | Verified | Non-text control embedding (B1) |
| 13 | Hu, Alkhairy, Lee, et al. "Evaluation of large language models for discovery of gene set function." *Nature Methods* 22:82–91 (2025), doi:10.1038/s41592-024-02525-x | Verified | LLMs name gene-set functions; memorisation question |
| 14 | Gillis & Pavlidis. "The Impact of Multifunctional Genes on 'Guilt by Association' Analysis." *PLoS ONE* 6(2):e17258 (2011); and "'Guilt by Association' Is the Exception Rather Than the Rule in Gene Networks." *PLoS Comput Biol* 8(3):e1002444 (2012) | Verified | Multifunctionality/study bias confound (baseline B4) |
| 15 | Radivojac, Clark, Oron, et al. *Nature Methods* 10:221–227 (2013), doi:10.1038/nmeth.2340; Zhou, Jiang, Bergquist, et al. *Genome Biology* 20:244 (2019), doi:10.1186/s13059-019-1835-8 | Verified | CAFA time-delayed evaluation |
| 16 | Wadi, Meyer, Weiser, Stein, Reimand. *Nature Methods* 13:705–706 (2016), doi:10.1038/nmeth.3963; Tomczak, Mortensen, Winnenburg, et al. *Scientific Reports* 8:5115 (2018), doi:10.1038/s41598-018-23395-2 | Verified | Annotation drift over time |
| 17 | Kapoor & Narayanan. "Leakage and the reproducibility crisis in machine-learning-based science." *Patterns* 4(9):100804 (2023), doi:10.1016/j.patter.2023.100804 | Verified | Leakage taxonomy for framing |

Novelty search result: no paper found that masks GO/pathway terms in NCBI gene summaries before embedding. That is a search result, not proof; do not write "first".
