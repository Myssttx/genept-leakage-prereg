# Stage 1 handback: go/no-go checks, environment and data

Executor: Claude Code (claude-opus-5-5). Folder: `~/Documents/genept-leakage`. Dates: 2026-10-04 (local CDT; UTC in fetch logs).
No classifier was trained, and no AUROC, AUPRC, accuracy or similarity score was computed on real data. No gene summary was embedded; the only text embedded was `'placeholder sentence'` (Step 6).

## 0. Gate answers

- **G1 (IDASB 2026 deadline):** the three sources conflict:
  - IDASB website, Key Dates: "Paper submission deadline: Sept 27, 2026".
  - BIBM workshop list: "Workshop paper submission due: Sep 27th, 2026".
  - IDASB flyer PDF (the workshop *proposal*), linked from the IDASB Call for Papers page: "October 15, 2026: Due date for full workshop papers submission".
  - BIBM CyberChair IDASB page: "Submission deadline: to be announced soon." It shows an active "Submit a New Paper" link. Of the 58 workshops in the CyberChair list, only REGSV is marked "(closed!)".
- **G2 (does GenePT B.3 mask or remove function words from summaries?):** **no**. B.3 has two checks: a "Temporal Split" argument based on release dates, and "Content Overlap", which counts interaction pairs. For Lit-BM it removed test *pairs*, not summary words: "we removed the positive pairs present in the NCBI summary and evaluated the performance again."
- **G3 (template share, copied from the report):**
  - `NCBI_UniProt_summary_of_genes.json: alliance_template=375 (1.1%), curated=26882 (79.8%), mixed=6446 (19.1%)`
  - `NCBI_summary_of_genes.json: alliance_template=2578 (7.6%), curated=28436 (84.4%), mixed=2689 (8.0%)`
  - These shares are probably affected by a "Gene Symbol <SYM>" prefix on every summary (see §6 P5 and §7 Q3).
- **G4 (NCBI spot-check):** **3 of 3** `alliance_template` genes (PLGLB2, RGPD4, BBX) currently carry "[provided by Alliance of Genome Resources, ...]". The 3 `mixed` genes do too. These are *current* NCBI records (attributions dated Jul 2025 / Jun 2026), not the 2024 snapshot (§6 P9).
- **G5 (GO baseline release):** **2024-03-28**.
- **G6 (free disk space at the end):** 20Gi (`/dev/disk3s5 460Gi 402Gi 20Gi 96%`), see §3.10.

## 1. Summary

I created the project from the prompt (all 7 SHA-256 checks `OK`) and built the Python 3.12.7 venv; the tests gave `4 passed`. The first disk check found only 2.9Gi free, so I stopped. After the user freed space (28Gi free), I resumed.

Data downloads:
- All 8 data items were downloaded, checksummed and described in `data/README.md`.
- The GenePT zip MD5 matches Zenodo.
- The GO baseline is release 2024-03-28.
- Both local embedding models load and return `(1, 384)`.

The descriptive report ran to `END OF REPORT` with no `FAILED` line.

Gate findings:
- **GenePT Appendix B.3 does not mask terms**, so the masking question is still open.
- The IDASB deadline is **conflicting**: Sept 27 on the website vs Oct 15 in the proposal flyer, and "to be announced soon" in CyberChair, where the IDASB submit link is still active.
- The NCBI spot-check confirms Alliance attribution for the sampled template genes.

Not completed:
- The GenePT v2 bioRxiv supplementary-material page returned HTTP 429 on every attempt (13 tries, 05:56–06:27Z; §6 P2). B.3 itself is in the main PDF, so no gate answer depends on that page.
- Step 4c (Zhong et al. v2) succeeded on a later retry and is quoted in §3.5.

## 2. Commands run

In order. Working directory is `~/Documents/genept-leakage` unless stated otherwise. Long helper scripts are saved in `handback/logs/`.

```
# Step 0
mkdir -p ~/Documents/genept-leakage/code/tests && cd ~/Documents/genept-leakage
# (7 files written with the Write tool)
cat > /tmp/stage1.sha256 <<'SUMS' ... SUMS        # exactly as in the prompt
shasum -a 256 -c /tmp/stage1.sha256

# Step 1 (first attempt: 2.9Gi free -> STOP; user freed space; resumed)
df -h .
mkdir -p handback && cp -r code handback/stage1_code_before
df -h .                                            # second check: 28Gi

# Step 2
python3 --version                                  # Python 3.12.7 (/opt/anaconda3/bin/python3)
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
pip freeze > handback/stage1_pip_freeze.txt
python -m pytest -q code/tests | tee handback/stage1_pytest.txt

# Step 5 (started in background, in parallel with Steps 3-4)
mkdir -p data/downloads data/raw/genept data/raw/go/current data/raw/go/baseline data/papers handback/logs
./handback/logs/download_main.sh > handback/logs/download_main.log 2>&1 &   # script saved in handback/logs; items 1,2,3,5,6,8,7 via curl -fSL --retry 3
curl -sSL https://release.geneontology.org/ > handback/logs/go_release_listing.html
grep -oE '>[0-9]{4}-[0-9]{2}-[0-9]{2}<' handback/logs/go_release_listing.html | tr -d '<>' | awk '$1>="2023-12-01" && $1<="2024-12-31"'
curl -sSIL https://release.geneontology.org/2024-03-28/annotations/goa_human.gaf.gz    # HTTP/2 200
curl -sSIL https://release.geneontology.org/2024-03-28/ontology/go-basic.obo           # HTTP/2 200
curl -fSL --retry 3 -o data/raw/go/baseline/goa_human.gaf.gz https://release.geneontology.org/2024-03-28/annotations/goa_human.gaf.gz
curl -fSL --retry 3 -o data/raw/go/baseline/go-basic.obo https://release.geneontology.org/2024-03-28/ontology/go-basic.obo
#   inside download_main.sh, gene2pubmed exactly as in the prompt:
#   curl -sSL https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz | gunzip -c | awk -F'\t' 'NR==1 || $1=="9606"' > data/raw/gene2pubmed_human.tsv
#   wc -l data/raw/gene2pubmed_human.tsv

# Step 3 (all pages saved under data/web/)
curl -sSL -o data/web/idasb_index.html 'http://www.combio-lezhang.online/IDASB2026/index.html'
curl -sSL -o data/web/idasb_<f> 'http://www.combio-lezhang.online/IDASB2026/js/<f>'   # app.5643eed8.js, demo.9ef88cf9.js, chunk-vendors.e26e98df.js
python3 (inline) -> extract Vue text nodes _v("...") from demo.9ef88cf9.js -> data/web/idasb_textnodes.txt
python3 (inline) -> whole-bundle regex search of all 3 bundles for: extended|extension|October|Oct d|Sept|deadline|blind|anonym|onsite|on-site|hybrid|in-person|virtual|present|registration|register|page limit|extra page|template|latex|IJCBDD|cyberchair|wi-lab|2026
curl -sSL -o data/web/idasb_IDASB_2026.pdf 'http://www.combio-lezhang.online/IDASB2026/IDASB_2026.pdf'
/opt/anaconda3/bin/python -c "import pypdf ..."       # text of flyer (pypdf 6.16.2 from system Python, not the project venv)
curl -sSL -o data/web/bibm_index.html 'https://www3.cs.stonybrook.edu/~bibm2026/'
curl -sSL -o data/web/bibm_<p> 'https://www3.cs.stonybrook.edu/~bibm2026/<p>'   # menu.html, home.html, CallforWorkshopProposals.html, workshopwebpges.html, BIBM2026-registration0.html, callforpaper2026.html
curl -sSL -o data/web/bibm_ws_submit.html 'https://wi-lab.com/cyberchair/2026/bibm26/scripts/ws_submit.php'
curl -sSL -o data/web/bibm_ws_S44.html 'https://wi-lab.com/cyberchair/2026/bibm26/scripts/submit.php?subarea=S44&undisplay_detail=1&wh=/cyberchair/2026/bibm26/scripts/ws_submit.php'
python3 handback/logs/html2text.py <file>              # plain-text dump used for quoting

# Step 4a
curl -sSL -A '<Chrome UA>' -o data/papers/genept_biorxiv_v2.full.pdf 'https://www.biorxiv.org/content/10.1101/2023.10.16.562533v2.full.pdf'   # 200
/opt/anaconda3/bin/python (pypdf) -> data/papers/genept_biorxiv_v2.full.txt
grep -n -i 'B\.3\|leak\|Appendix B\|information leakage' data/papers/genept_biorxiv_v2.full.txt
sed -n '455,468p' / sed -n '760,860p' data/papers/genept_biorxiv_v2.full.txt
sips -s format png (single-page PDF of p.20) ; PIL crop of Fig. B3 legend   # to read AUC legend values
curl -sSL -A '<Chrome UA>' -o data/papers/genept_biorxiv_v2_supplementary-material.html 'https://www.biorxiv.org/content/10.1101/2023.10.16.562533v2.supplementary-material'   # 429 (x3, then background retry loop, see §3.5)

# Step 4b
curl -sSL -A '<Chrome UA>' 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7304461/'         # 200 but reCAPTCHA page
curl -sSL 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=PMCID:PMC7304461&format=json&resultType=core'
curl -sSL -A '<Chrome UA>' 'https://doi.org/10.1093/database/baaa037'                  # 403 (academic.oup.com "Just a moment...")
curl -sSL -o data/papers/kishore2020_PMC7304461_fullText.xml 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7304461/fullTextXML'   # 200
python3 (inline XML parse) -> data/papers/kishore2020_PMC7304461_fullText.txt ; sentence grep for NCBI|display|distribut|download|evidence|...

# Step 4c
curl -sSL -A '<Chrome UA>' -o data/papers/zhong2025_biorxiv_v2.full.html 'https://www.biorxiv.org/content/10.1101/2025.01.29.635607v2.full'      # 429
curl -sSL -A '<Chrome UA>' -o ... 'https://www.biorxiv.org/content/10.1101/2025.01.29.635607v2.full.pdf'                                         # 429
curl -sSL 'https://api.biorxiv.org/details/biorxiv/10.1101/2025.01.29.635607'          # metadata only
curl -sSL 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22<doi>%22&format=json&resultType=core'   # PPR973307 / PPR745627, no full text
curl -sSL 'https://www.ebi.ac.uk/europepmc/webservices/rest/PPR973307/fullTextXML'     # 500
(Chrome extension tried: "Browser extension is not connected")
background loop: every 120 s, retry both bioRxiv URLs -> handback/logs/biorxiv_retry.log   # 429 at 06:04 and 06:06; stopped
WebFetch of the .full and .full.pdf URLs                                         # both "HTTP 429 Too Many Requests", Retry-After: 19 / 14
./handback/logs/biorxiv_retry2.sh | tee -a handback/logs/biorxiv_retry.log      # honours Retry-After, browser Accept headers; Zhong .full -> 200 at 06:08:31Z
python3 handback/logs/html2text.py data/papers/zhong2025_biorxiv_v2.full.html > data/papers/zhong2025_biorxiv_v2.full.txt
grep -n -i 'temporal\|holdout\|cutoff\|March 2024\|2024\|release' data/papers/zhong2025_biorxiv_v2.full.txt

# Step 5 (continued)
unzip -l data/downloads/GenePT_emebdding_v2.zip > data/raw/genept/zip_listing.txt
unzip -j -o data/downloads/GenePT_emebdding_v2.zip -d data/raw/genept/      # -j = flat; no renames
df -h .
for f in <all data files>; do stat -f %z; shasum -a 256; done > handback/logs/data_checksums.tsv
curl -sSL https://zenodo.org/api/records/10833191 > data/web/zenodo_10833191.json ; md5 -q data/downloads/GenePT_emebdding_v2.zip
python3 (inline) -> data/README.md ; printf 'data/\n.venv/\n' > .gitignore

# Step 6
export HF_HOME="$PWD/data/models"
python -c "from sentence_transformers import SentenceTransformer as S; [print(m, S(m).encode(['placeholder sentence']).shape) for m in ['sentence-transformers/all-MiniLM-L6-v2','BAAI/bge-small-en-v1.5']]" > handback/stage1_step6_models.txt 2> handback/logs/step6_stderr.txt

# Step 7
python code/stage1_data_report.py > handback/stage1_data_report.txt 2>&1

# Step 8 (first pass discarded: zsh parsed $g[sym] as an array subscript; see §6 P12)
curl -g -sSL -A '<Chrome UA>' 'https://www.ncbi.nlm.nih.gov/gene/?term=PLGLB2[sym]+AND+human[orgn]'   # 200 but reCAPTCHA page
for g in PLGLB2 RGPD4 BBX RAB39A NSF ZSCAN18 MIR3925 RNA5SP87 HPYR1; do
  curl -g -sSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${g}%5Bsym%5D+AND+human%5Borgn%5D&retmode=json"
  curl -sSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=${ID}&retmode=xml"     # read <Entrezgene_summary>
done | tee handback/logs/step8_ncbi_eutils.txt

# Extra descriptive check (not requested; see §5 D6)
python (venv) using stage1_data_report.load_pickle_embeddings on GenePT_gene_protein_embedding_model_3_text.pickle. -> handback/logs/extra_model3_pickle_describe.txt

# Wrap-up
diff -ru handback/stage1_code_before code ; diff -ru -x __pycache__ handback/stage1_code_before code
df -h .
```

## 3. Outputs

### 3.1 Step 0 checksums
```
CLAUDE.md: OK
STAGES.md: OK
review.md: OK
plan.md: OK
requirements.txt: OK
code/stage1_data_report.py: OK
code/tests/test_stage1_report.py: OK
```

### 3.2 Step 1 disk checks
First attempt (stopped here as instructed):
```
Filesystem      Size    Used   Avail Capacity iused ifree %iused  Mounted on
/dev/disk3s5   460Gi   421Gi   2.9Gi   100%    3.4M   31M   10%   /System/Volumes/Data
```
After the user freed space (Sun Oct 4 00:52:13 CDT 2026):
```
Filesystem      Size    Used   Avail Capacity iused ifree %iused  Mounted on
/dev/disk3s5   460Gi   396Gi    28Gi    94%    3.1M  299M    1%   /System/Volumes/Data
```

### 3.3 Step 2 environment and pytest
`Python 3.12.7`. `pip install -r requirements.txt` succeeded with no failures. Key versions: numpy==2.5.3, scipy==1.18.1, pandas==3.0.6, scikit-learn==1.9.1, torch==2.14.1, sentence-transformers==6.1.0, transformers==5.18.0, pytest==9.1.1. The full list (47 lines) is in `handback/stage1_pip_freeze.txt`.
```
....                                                                     [100%]
4 passed in 1.28s
```

### 3.4 Step 3: venue facts (verbatim quotes)

**IDASB 2026 site.** It is a Vue single-page app with no server-rendered text, so all visible text comes from the JS bundle `/IDASB2026/js/demo.9ef88cf9.js`. I extracted its text nodes and searched the whole bundle (`app`, `demo`, `chunk-vendors`) for deadline, extension, review-type, format and registration words.
- Fetched 2026-10-04T05:53:30Z (index.html: `<title>IDASB 2026 (IEEE BIBM workshop on Integrative Data Analysis in Systems Biology)</title>`) and 05:53:37Z (bundles).
- Routes in the bundle: `/Home`, `/AboutIDASB`, `/CallForPapers`, `/Committees`, `/Contact`, `/KeyDates`, `/PaperSubmission`, `/Program`, `/Registration`, `/Venue`.

Key Dates table (the same table occurs twice in the bundle):
> Paper submission deadline: / Sept 27, 2026
> Paper acceptance notification: / Oct 18, 2026
> Camera-ready manuscript due: / Nov 8, 2026
> Workshop dates: / Dec 1-4, 2026

Paper Submission:
> Please submit a full-length paper (8 page IEEE 2-column format) through the online submission system. Electronic submissions in pdf format are required.
> Please submit your paper through the BIBM 2026 online submission system [link: https://wi-lab.com/cyberchair/2026/bibm26/index.php].

Registration:
> Online registration information will be available from the BIBM 2026 main conference website [link: https://www3.cs.stonybrook.edu/~bibm2026/].

Call for Papers / Home:
> We invite you to submit papers with unpublished, original research describing recent advances on the areas related to this workshop. All papers will undergo peer review by the conference program committee. All papers accepted will be included in the Workshop Proceedings published by the IEEE Computer Society Press and will be available at the workshops.

Topics (two of them):
> Data science approaches and applications
> Large language models and applications in biological and biomedical fields

Venue:
> The IEEE BIBM 2026 main conference will be held in Dallas, USA, December 1-4, 2026.

Program:
> The IDASB 2026 workshop program will be updated soon.

Contact:
> 1. The workshop website will be alive when the proposal is accepted.
> 2.The selected papers will be published in a special issue in IJCBDD, as we did in the past ten years.

Not found anywhere in the IDASB bundle (whole-bundle search): any extension notice, the review type (single- or double-blind), the presentation format (onsite, online or hybrid), any extra-page fee, and any template detail beyond "8 page IEEE 2-column format".

**IDASB flyer PDF** `http://www.combio-lezhang.online/IDASB2026/IDASB_2026.pdf`, linked from the Call for Papers page as "flyer in pdf format". Fetched 2026-10-04T05:54:22Z. It is 2 pages; page 1 is headed "Proposal for 2026 IEEE BIBM Workshop".
> 4. Important dates
> October 15, 2026: Due date for full workshop papers submission
> November 15, 2026: Notification of paper acceptance to authors
> November. 23, 2026: Camera-ready of accepted papers
> December 02, 2026: Workshop.

> 8. Notes to BIBM workshop chair:
> 1. The workshop website will be alive when the proposal is accepted.
> 2. The selected papers will be published in a special issue in IJCBDD, as we did in the past ten years.

**BIBM 2026 home page** `https://www3.cs.stonybrook.edu/~bibm2026/`. It is a frameset; the content is in `home.html` and the menu in `menu.html`. Fetched 05:54:43Z and 05:54:48Z.
> Important Dates
> Electronic submission of full papers: July 5, 2026
> Notification of paper acceptance: Sept 25, 2026
> Camera-ready of accepted papers:Oct 25, 2026
> Conference: Dec 1-4, 2026

Menu items include "Workshop Paper Submission" → `https://wi-lab.com/cyberchair/2026/bibm26/scripts/ws_submit.php` and "Workshop Webpages" → `workshopwebpges.html`.

**BIBM Call for Workshop Proposals** `https://www3.cs.stonybrook.edu/~bibm2026/CallforWorkshopProposals.html`, fetched 2026-10-04T05:55:03Z.
> All papers accepted for workshops will be included in the main Conference Proceedings published by the IEEE Computer Society Press, made available at the Conference.
> The workshop organizers will also have the discretion of editing selected papers (after their expansion and revision) into books or special journal issues. Workshops should cover at least a single session (8 regular papers). The workshop organizers should ensure registration and presence of authors of accepted papers.
> Workshop style
> Full onsite, full online, or hybrid
> Organizers can select the workshop style
> If the organizers decide to organize full online or hybrid workshop, they need to take care of the virtual meeting platform for the workshop
> Important dates (Please do not change the final camera-ready submission date, feel free to adjust other dates)
> Sept 27, 2026: Due date for full workshop papers submission
> Oct 18, 2026: Notification of paper acceptance to authors
> Nov 8, 2026: Camera-ready of accepted papers
> Dec 1-4, 2026: Workshops

**BIBM Workshop Webpages** `https://www3.cs.stonybrook.edu/~bibm2026/workshopwebpges.html`, fetched 2026-10-04T05:55:03Z. The IDASB entry, in full:
> WS#44 | The 17th Integrative Data Analysis in Systems Biology and Data Science (IDASB 2026) (IDASB2026)
> Main organizer: Huiru Zheng
> Main organizer contact: h.zheng@ulster.ac.uk
> All organizers: Huiru Zheng, Ming Xiao, Zhongming Zhao
> Workshop website: http://www.combio-lezhang.online/IDASB2026/index.html#/Home
> Workshop paper submission due: Sep 27th, 2026
> Notification date: Oct 18th, 2026
> Camera-ready deadline: Nov 8th, 2026

**BIBM workshop submission system** (CyberChair, linked from the BIBM menu; on wi-lab.com, not the stonybrook site).
- `https://wi-lab.com/cyberchair/2026/bibm26/scripts/ws_submit.php`, fetched 2026-10-04T05:55:03Z. There are 58 workshop entries. Exactly one line contains "closed!":
  > Workshop on Regulatory Structural Variation in Disease (REGSV) (closed!)

  The IDASB line is a live link:
  > 44 | The 17th Integrative Data Analysis in Systems Biology and Data Science (IDASB 2026)
- The IDASB entry page, `.../scripts/submit.php?subarea=S44&...`, fetched 2026-10-04T05:55:20Z:
  > Submission deadline: to be announced soon.
  > Submit a New Paper
  > Submissions must follow IEEE Computer Society Proceedings Manuscript Formatting Guidelines ( https://www.ieee.org/conferences/publishing/templates.html), allowing up to 8 pages for full papers and up to 4 pages for short papers.
  > Note: If submitting a short paper, please double-check that your specific workshop accepts them.
  > If your paper was submitted but not accepted in the main conference, you cannot directly submit it to a workshop but it must be transferred into the workshop. You can change this option (workshop) at the conference paper submission site or (after the deadline) request the organizer of the workshop to which you want your paper to be transferred.
  > The camera ready submission site is closed.

**BIBM main Call for Papers** `https://www3.cs.stonybrook.edu/~bibm2026/callforpaper2026.html`, fetched 2026-10-04T05:55:03Z. This page is for the main conference.
> Please submit a full-length paper (up to 8 page IEEE 2-column format) through the online submission system (you can download the format instruction here: http://www.ieee.org/conferences_events/conferences/publishing/templates.html).The 8-page limit includes all content, including the main paper, references, and any appendices. Authors do not need to include an Ethics Statement or an Acknowledgment section in the submission. The paper submission is double blind. Please don’t include any authors and/or affiliations in the paper.
> Electronic submissions (in PDF or Postscript format) are required. Selected participants will be asked to submit their revised papers in a format to be specified at the time of acceptance.
> BIBM 2026 will offer as many student travel awards (around 40 awards) to student authors as possible based on the sponsorship fees (including post-doc). All student authors are eligible to apply for the student travel awards.

**BIBM 2026 registration** `https://www3.cs.stonybrook.edu/~bibm2026/BIBM2026-registration0.html`, fetched 2026-10-04T05:55:03Z.
> Join us for the 2026 IEEE International Conference on Bioinformatics and Biomedicine, December 1–4, 2026, at the Hyatt Regency Dallas in Dallas, Texas.
> Early registration deadline: November 8, 2026, at 11:59 PM Central Time.
> All fees are in U.S. dollars (USD).

> Conference Registration Fees

| Registration Category | Author Registration | Early Registration Through November 8 | Standard Registration After November 8 |
|---|---|---|---|
| IEEE Member | $800 | $800 | $900 |
| Non-Member | $960 | $960 | $1,080 |
| IEEE Student Member | $560 | $560 | $630 |
| Student Non-Member | $675 | $675 | $760 |
| IEEE Life Member | $440 | $440 | $495 |

> For authors: Select Conference Registration – Author in the registration system. Please review Important Author Information in the registration system for paper registration requirements and author deadlines. Author registrations are non-refundable.
> For other attendees: Select Conference Registration. Register by November 8 to receive the early rate.

> Additional Papers — Authors registering additional papers should select the appropriate option during registration.

| Registration Category | 1 Additional Paper | 2 Additional Papers | 3 Additional Papers* | 4 Additional Papers* |
|---|---|---|---|---|
| IEEE Member | $400 | $1,000 | $1,800 | $2,600 |
| Non-Member | $400 | $1,000 | $1,960 | $2,920 |
| IEEE Student Member | $400 | $1,000 | $1,560 | $2,120 |
| Student Non-Member | $400 | $1,000 | $1,675 | $2,350 |
| IEEE Life Member | $400 | $1,000 | $1,440 | $1,880 |

> *For details on registering three or four additional papers, please visit the registration page.
> Additional Pages — Additional page fee: $100 per page.
> Select additional pages separately for each paper during registration. Please check the page limits and additional-page allowances for your paper category before submitting payment.
> Cancellation Policy — Refund and cancellation requests must be emailed to ieeecs-reg+BIBM@computer.org by November 8, 2026, at 11:59 PM Central Time.
> A $100 administrative fee applies to cancelled registrations.
> Author registrations are non-refundable.
> Questions? For assistance with registration, fees, or payment, please contact ieeecs-reg+BIBM@computer.org.

Student eligibility or proof of student status: **not stated** on this page. I checked the raw HTML: "student" appears only in the table row labels.

### 3.5 Step 4: literature (verbatim quotes)

#### 4a. GenePT Appendix B.3
Source: `https://www.biorxiv.org/content/10.1101/2023.10.16.562533v2.full.pdf`, fetched 2026-10-04T05:55:40Z. HTTP 200, 15,940,010 bytes. Page footer: "this version posted March 5, 2024". B.3 is in the main PDF (pages 18–20), so the Nature BME version was not needed.

Main text, §4.2 end:
> Finally, it’s crucial to confirm that the promising results in Sections 4.1 and 4.2 are not simply the result of information leakage, such as test set data being included in the original NCBI gene summaries used as input for GenePT. We address these concerns in detail in Appendix B.3.

**Appendix B.3, in full:**
> **B.3 Addressing potential information leakage**
> It’s important to ensure that the promising results in Section 4 are not merely due to the test set being represented in the original NCBI gene summaries used as input for GenePT gene embeddings. We address the concerns of information leakage by using a temporal split for gene functionality prediction and quantifying the extent of information leakage for gene-gene and protein-protein interaction predictions.
>
> 1. Temporal Split: In terms of the temporal split, we discovered that the gene functionality data was released by Theodoris et al. [1], in June 2023 at https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/example input files. Conversely, the NCBI gene information does not contain information on these specific data, and the embedding model was released in December 2022. Therefore, the information leakage for the gene property prediction should be minimal. Regarding gene functionality, we intended this as a sanity check to ensure that our embeddings contain first-order information, and we have adjusted our original statement accordingly.
>
> 2. Content Overlap: We have also quantified the extent to which explicit gene-gene interactions and protein-protein interactions are mentioned in the NCBI gene summaries (see Table B3). This ensures that we are not merely memorizing training data. In this context, we quantified, for gene-gene interaction and protein-protein interaction datasets, the number of explicitly mentioned interaction pairs. We note that the number of present pairs is quite small across these datasets. In particular, except for the Lit-BM PPI dataset , the small percentage (less than 1%) of leakage in the other three datasets had a minimal impact on the classification result.
>
> Additionally, to quantify the impact of data leakage in the Lit-BM PPI dataset (approximately 4%), we removed the positive pairs present in the NCBI summary and evaluated the performance again. We observed that there is no meaningful difference between the performance with and without the pairs that were mentioned in the NCBI summary. Therefore, this establishes that our GenePT approach was able to encode more information than mere memorization.
>
> | Dataset | Positive pairs | Pairs present in NCBI summary |
> |---|---|---|
> | GGI [18] | 141,130 | 310 |
> | HuRI PPI [33] | 53,548 | 583 |
> | Lit-BM PPI [34] | 12,234 | 524 |
> | Heart tissue PPI [35] | 303,932 | 469 |
>
> **Table B3**: Quantifying potential information leakage in gene-gene interaction and protein-protein interaction datasets.
>
> [Fig. B3 panel titles: "ROC for PPI prediction on Lit-BM dataset (with leakage)" and "ROC for PPI prediction on Lit-BM dataset"]
> (a) Performance on full test set — legend: GenePT (AUC = 0.67); GenePT (Name only) (AUC = 0.62); Gene2vec (AUC = 0.60); Geneformer (AUC = 0.63); scGPT (AUC = 0.58); Permuted (AUC = 0.49)
> (b) Performance after removing overlap — legend: GenePT (AUC = 0.68); GenePT (Name only) (AUC = 0.63); Gene2vec (AUC = 0.60); Geneformer (AUC = 0.61); scGPT (AUC = 0.55); Permuted (AUC = 0.46)
>
> **Fig. B3**: Test set prediction performance using GenePT gene embeddings on the Lit-BM protein-protein interaction datasets [34] before (panel (a)) and after (panel (b)) removing the mentioned pairs of interacting proteins in the NCBI gene summary database.

I read the Fig. B3 legend values from a high-resolution render of page 20 because pypdf text extraction does not capture figure text.

For context, GenePT's own task description (§3):
> – Gene Functionality Class Prediction: This is a multi-class prediction challenge based on the 15 most common functional gene classes. Labels for these classes were curated as part of the Geneformer paper.
> – Gene-Gene Interaction Prediction: We utilized a benchmark for gene-gene interaction (GGI) based on shared gene ontology annotations published by Du et al. [18].

GenePT methods, on the summary source:
> For each gene, we extracted its information from the NCBI gene database’s summary section after removing hyperlinks and date information.

Answers:
1. **Removed, masked or replaced words in the summaries?** Not stated; no word removal or masking is described. What they removed was test-set *pairs*: "we removed the positive pairs present in the NCBI summary and evaluated the performance again."
2. **Tasks and labels:**
   - Temporal-split argument: "gene functionality prediction" and "gene property prediction" (Theodoris et al. data).
   - Content overlap: "GGI [18]", "HuRI PPI [33]", "Lit-BM PPI [34]", "Heart tissue PPI [35]".
   - GO terms as labels in B.3: not stated. §3 describes the GGI benchmark as "based on shared gene ontology annotations published by Du et al. [18]".
3. **Dates or time splits?** Only a release-date argument, no date-based split of labels: "the gene functionality data was released by Theodoris et al. [1], in June 2023 ... the embedding model was released in December 2022. Therefore, the information leakage for the gene property prediction should be minimal."
4. **Numbers before and after:**
   - Table B3 counts (above), plus "less than 1%" for three datasets and "approximately 4%" for Lit-BM.
   - Lit-BM ROC AUC before → after removing overlap: GenePT 0.67 → 0.68; GenePT (Name only) 0.62 → 0.63; Gene2vec 0.60 → 0.60; Geneformer 0.63 → 0.61; scGPT 0.58 → 0.55; Permuted 0.49 → 0.46.
5. **Conclusion, in their words:** "We observed that there is no meaningful difference between the performance with and without the pairs that were mentioned in the NCBI summary. Therefore, this establishes that our GenePT approach was able to encode more information than mere memorization." Also: "Therefore, the information leakage for the gene property prediction should be minimal."

**Supplementary-material page:** NOT fetched. Every request to `https://www.biorxiv.org/content/10.1101/2023.10.16.562533v2.supplementary-material` returned HTTP 429 (17-byte body `error code: 1015`): at 05:56:19Z, 05:57:04Z and 06:01:07Z, then 10 more times in the two background retry loops (`handback/logs/biorxiv_retry.log`, copied in §3.5c). Since B.3 was found in the main PDF, its absence does not affect the B.3 quote.

#### 4b. Kishore et al. 2020 (Alliance automated summaries)
- `https://pmc.ncbi.nlm.nih.gov/articles/PMC7304461/` (2026-10-04T05:56:25Z) returned "Checking your browser - reCAPTCHA".
- `https://doi.org/10.1093/database/baaa037` redirected to `https://academic.oup.com/database/article/doi/10.1093/database/baaa037/5859735` with HTTP 403 "Just a moment... Enable JavaScript and cookies to continue" (05:56:35Z).
- So I quoted from the Europe PMC full-text XML of the same PMC record: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7304461/fullTextXML`, fetched 2026-10-04T05:56:45Z (HTTP 200, 121,333 bytes).

Bibliographic data (from the XML `article-meta`):
- **Title:** Automated generation of gene summaries at the Alliance of Genome Resources.
- **Authors, in full:** Ranjana Kishore; Valerio Arnaboldi; Ceri E Van Slyke; Juancarlos Chan; Robert S Nash; Jose M Urbano; Mary E Dolan; Stacia R Engel; Mary Shimoyama; Paul W Sternberg; the Alliance of Genome Resources (collab).
- **Journal:** Database: the journal of biological databases and curation. Volume 2020; elocation-id (article number) baaa037; epub 2020-06-19.
- **DOI:** 10.1093/database/baaa037 (PMID 32559296, PMCID PMC7304461).

Where the descriptions are displayed or distributed:
> In addition to being displayed on Alliance gene pages, these summaries are also included on several MOK gene pages.
> They can be downloaded in bulk from the Downloads page on the Alliance website (https://alliancegenome.org/downloads#gene-descriptions), and they are displayed on current Alliance gene pages.
> In addition, Alliance summaries are currently imported for display on MOK gene pages that had no gene summaries in the past (ZFIN), supplement existing summaries at the other MOKs (MGD, WB, FB and RGD), and are used to enhance existing curated summaries (SGD).
> Our solution to this problem has resulted in automated, modular and standardized gene summaries for seven species that are updated with each new release of the Alliance, displayed on Alliance gene pages and available for download.
> The gene summaries software is used at the Alliance and is integrated into the build process of the web portal to generate a summary on each gene page, and species-specific summary files for download in different file formats, including the custom gene summary JSON format defined at the Alliance (https://github.com/alliance-genome/agr_schemas/tree/master/ingest/genedescription).

**NCBI Gene as a display or distribution site: not stated.** The only mention of NCBI in the full text:
> The automated gene summaries generated by our method include gene functional data for seven species, six MOK species, in addition to human data provided by RGD, which maintains a full set of human gene records and annotations imported and integrated from diverse sources such as NCBI (19), Ensembl (20), UniProt-GOA (21) annotations mapped to HGNC IDs (HUGO Gene Nomenclature Committee; 22) and Online Mendelian Inheritance in Man (23).

Which data the summaries are generated from, and which evidence types:
> The gene data used for this purpose include curated associations (annotations) to ontology terms from the Gene Ontology, Disease Ontology, model organism knowledgebase (MOK)-specific anatomy ontologies and Alliance orthology data.
> We generated automated gene summaries using gene annotations to ontology terms such as GO, DO, etc. and additional gene-related data such as orthologs.
> Note that all data used in the gene summaries are obtained from the Alliance File Management System API (https://fms.alliancegenome.org) and from the Alliance database at build ime.
> All gene annotations to ontology terms that MOKs submit to the Alliance for a gene of interest were considered for generating gene summaries. However, in order to maintain readability of the final description, when a gene had annotations supported by experimental evidences (according to ECO; 11), these annotations were preferred over others. When a gene had no experimental annotations, annotations supported by other types of evidences were included such as the following (in order of preference): high-throughput experimental evidence codes, phylogenetic evidence codes, curator and author statement evidence codes, computational evidence codes and electronic evidence codes (8, 11, 28).
> For annotations based on non-experimental evidence codes, such as those derived from phylogenetic, sequence based or computational analyses [e.g. GO annotations that are a result of projects such as PAINT and INTERPRO2GO (28, 30, 31)] the phrase ‘Predicted to’ is prepended to the verb phrase of the template, e.g. Predicted to exhibit <list of GO molecular function terms>.
> Templates were built for the generation of sentences for each data category by using the actual ontology terms that the gene of interest is annotated to, together with a verb phrase such as ‘exhibits’, ‘involved in’, ‘is expressed in’, etc.
> Note that, even though each data category can only have either experimental or predicted data due to the evidence code prioritization, the final summary may contain a mixture of experimental and predicted annotations across different categories, as is the case for the C. elegans gene cdk-4, shown in Figure 2.

#### 4c. Zhong et al. v2 temporal holdout
Source: `https://www.biorxiv.org/content/10.1101/2025.01.29.635607v2.full`. Earlier attempts returned HTTP 429. Fetched successfully at 2026-10-04T06:08:31Z (HTTP 200, 285,825 bytes) by the retry script `handback/logs/biorxiv_retry2.sh`. Page: "Benchmarking gene embeddings from sequence, expression, network, and text models for functional prediction tasks". Authors: Jeffrey Zhong, Lechuan Li, Ruth Dannenfelser, Vicky Yao; doi: https://doi.org/10.1101/2025.01.29.635607; "Posted November 10, 2025". Quotes are from the HTML full text. My text dump marks hyperlinks as `[link:URL]`; I removed those markers, so citations read as on the page (e.g. "GO [38] terms", "(Figure 1)").

Results (opening section):
> While our benchmarking suite is structured to measure different gene relationships and accounts for the multi-functionality of genes, the diverse data inputs used in constructing gene embeddings still create possible data leakage, particularly in settings where tasks involve the prediction of shared gene pathways or well-known annotations. Data leakage is less of a problem for methods that train primarily on a single data type, such as amino acid sequence and gene expression-based models, or even protein-protein interaction networks, but is potentially concerning in embeddings derived from biomedical literature and multi-input models that may use functional annotations. To mitigate this, we added a dedicated temporal holdout task for the single gene and gene set evaluations, particularly gene ontology (GO) annotations that were released after the latest embedding method we evaluate was trained (March 2024). Furthermore, for predicting pairwise gene relationships, we created an evaluation set of genetic interactions for a less-studied organism Saccharomyces pombe (S. pombe) which, to the best of our knowledge, was not used as input for any embedding method. Methods that perform well on these tasks are likely to be more generalizable.

Results › Gene-level attribute benchmarks:
> To evaluate how well each embedding captures individual gene-level attributes, we performed two primary classification tasks: gene function prediction using GO [38] terms and disease-gene prediction using OMIM [39] annotations. For each GO or OMIM term, we trained a separate SVM classifier for every embedding, using the embedding vectors as input features to predict held-out genes sharing that annotation (Figure 1). To focus on distinct, nonredundant biological categories and reduce overlapping or hierarchical annotations, we curated “slim” sets of GO Biological Process and Disease Ontology [40] terms (see Methods). Additionally, to prevent circularity and minimize data leakage, we included only GO annotations supported by experimental evidence codes and excluded any GO or OMIM term associated with fewer than 20 genes to ensure adequate sample size for training and testing. After these filters, the benchmark comprised 56 GO terms and 103 OMIM disease terms, with an additional 19 GO terms reserved as a temporal holdout to approximate generalizability.

> In the temporal holdout task (“temporal” column, Figure 2C), classifiers were trained on pre-2024 annotations and tested on genes annotated after March 2024. Because none of these updated gene annotations were included in the training data for any embedding method, this evaluation provides a leakage-free measure of generalization. As expected, performance declined in this more challenging setting (mean AUROC across all methods dropped from 0.79 to 0.62). Despite this drop, the relative ranking of models remained consistent, with PPI-based embeddings (Mashup, Node2Vec, PPI-Raw) showing the strongest generalization, performing better than literature-based models (mean AUROC GenePT-Model3: 0.774; Mashup: 0.777; Node2Vec 0.704; PPI-Raw 0.682). Performance varied across individual GO and OMIM terms, with some (e.g., “endodermal cell differentiation”, “congenital disorder of glycosylation”) consistently easier to predict, while others (e.g., “glucose homeostasis”, “spinal disease”, and “motor neuron disease”) being more difficult across all models (Figures S1 and S2).

Figure 2 caption ("Figure 2. Benchmarking results for single gene tasks."):
> Embedding methods were evaluated for their performance as features in predicting two classes of disease attributes: (A) gene function using GO and (B) disease gene identification using OMIM. In both cases embeddings are categorized and sorted by the median performance of their primary input data type with each dot representing one AUROC score for a single embedding type and GO term (A) or disease (B) pair. (C) Summarizes findings by ranking embedding methods according to the mean AUROC across tasks. Columns denoted as “full” use the entire set of genes covered in the released embeddings. The results for the GO temporal holdout is summarized in the temporal column. The time column highlights the differences in average training times for a single SVM classifier using the different embeddings.

Results › Gene set comparison benchmarks:
> As with all of our benchmarking tasks, we aimed to evaluate generalizability by including test cases unlikely to be included as training for the original embedding models. Here, we again designed a temporal holdout of GO terms, using terms with annotation updates before and after March 2024 where a minimum of 5 genes were added. For each such term, we applied the ANDES gene set matching test between its pre– and post-update versions, reasoning that they should remain functionally similar across updates. Although matching performance varied widely across terms (Figure S4), PPI-based embeddings again emerged on top (mean nrank, PPI-raw: 0.634, Node2Vec: 0.634; Figure 4C). The highest performing embeddings for the KEGG-GO and disease-tissue matching task also performed well in this temporal task, with BioConceptVec skip-gram (mean nrank, 0.632) and fast-text mean nrank, 0.629) having comparable performance to the PPI embeddings.

Discussion:
> In this study, we undertook a comprehensive benchmarking effort to evaluate gene embeddings across a wide range of gene-centered prediction tasks. One of the most notable findings is the strong performance of text-based embeddings and the overall weak performance of expression-based models across all tasks. GenePT-Model3 was among the top performers in nearly all evaluations except the disease-tissue gene set matching task, which is the least well-characterized. The impressive performance of text-based embeddings raises concerns about potential data leakage, since these embeddings are trained on curated documents or vast corpora of biomedical literature that may implicitly contain knowledge overlapping with our evaluation benchmarks. While it is challenging to completely separate the collective knowledge of biomedical research from our gold standards, we attempted to mitigate it with carefully designed temporal and organism-specific holdouts. In these benchmarks, we did not observe any major performance drops that would indicate substantial leakage. Interestingly, PPI-based methods performed especially well in these evaluations. We speculate that their explicit encoding of gene-gene relationships may enhance their generalizability to unseen data.

Methods › Gene ontology gold standards:
> The gold standard for individual test cases was created by taking different subsets of the GO BP hierarchy. For the single gene analyses, we used the 2024-09-08 release of GO BP and retained GO terms with more than 20 annotated genes. To create a nonredundant but representative evaluation, we then filtered these terms to a “slim” set of 56 terms using the Alliance of Genome Resources (AGR) slim annotations by randomly sampling three GO terms per AGR term. To evaluate generalizability over time, we also constructed a temporal holdout for use in both the single gene and gene set benchmarks. GO annotations were dived into two parts: those released prior to 2024-03, and those released by 2025-03-16. We then identified terms that gained at least 5 genes between these versions, resulting in a final set of 19 GO terms for the temporal evaluation.

Methods › Gene set comparison benchmark:
> For the temporal holdout task, we used the 19 GO terms rom the gene-level temporal benchmark. Each term’s annotations were divided into those annotated before 2024 and those annotated from 2024 onward (that had not been previously annotated). ANDES was then run to compare each pre– and post-2024 cohort.

Supplementary figure caption:
> Figure S4. Gene set matching results for GO temporal holdout. ANDES ranks for each GO term with itself pre and post the 2024 holdout partition. Lower ranks indicate a better match between old and new terms.

Summary of what the text says (quotes only):
- **Cutoff:** "annotations that were released after the latest embedding method we evaluate was trained (March 2024)"; "classifiers were trained on pre-2024 annotations and tested on genes annotated after March 2024".
- **GO releases:** "we used the 2024-09-08 release of GO BP"; "GO annotations were dived into two parts: those released prior to 2024-03, and those released by 2025-03-16".
- **New annotations:** "terms that gained at least 5 genes between these versions, resulting in a final set of 19 GO terms"; for gene sets, "annotated from 2024 onward (that had not been previously annotated)".
- **GenePT numbers:** "mean AUROC GenePT-Model3: 0.774; Mashup: 0.777; Node2Vec 0.704; PPI-Raw 0.682"; "mean AUROC across all methods dropped from 0.79 to 0.62".
- **Masking of summary words:** not described in any of these paragraphs.

Background retry log (`handback/logs/biorxiv_retry.log`):
```
2026-10-04T06:04:31Z zhong 429 17
2026-10-04T06:04:36Z supp 429 17
2026-10-04T06:06:36Z zhong 429 17
2026-10-04T06:06:41Z supp 429 17
2026-10-04T06:08:31Z zhong_full HTTP 200 size=285825 retry-after=none
2026-10-04T06:09:11Z supp HTTP 429 size=17 retry-after=103
2026-10-04T06:12:05Z supp HTTP 429 size=17 retry-after=78
2026-10-04T06:14:33Z supp HTTP 429 size=17 retry-after=109
2026-10-04T06:17:32Z supp HTTP 429 size=17 retry-after=68
2026-10-04T06:19:51Z supp HTTP 429 size=17 retry-after=69
2026-10-04T06:22:10Z supp HTTP 429 size=17 retry-after=104
2026-10-04T06:25:04Z supp HTTP 429 size=17 retry-after=77
2026-10-04T06:27:32Z supp HTTP 429 size=17 retry-after=57
(stopped manually at about 06:29Z, before the 25-minute limit, after 8 consecutive 429s on the supplement page)
```

### 3.6 Step 5: downloads
- GO release listing near the cutoff (from `https://release.geneontology.org/`, fetched 2026-10-04T05:53:16Z): `2024-01-17, 2024-03-28, 2024-04-18, 2024-04-24, 2024-06-10, 2024-06-17, 2024-09-08, 2024-11-03`. **Chosen baseline: 2024-03-28**, the earliest dated folder on or after 2024-03-18. Both files exist, with HTTP/2 200:
  - `goa_human.gaf.gz`: content-length 11886126, last-modified Sun, 07 Apr 2024 23:54:49 GMT
  - `go-basic.obo`: content-length 31267219, last-modified Sun, 07 Apr 2024 23:59:07 GMT
- All curl calls returned rc=0. gene2pubmed `pipestatus=0 0 0`; `wc -l`: `3339674 data/raw/gene2pubmed_human.tsv`.
- Gene2vec: the given raw GitHub URL worked (no 404).
- GenePT zip: `GenePT_emebdding_v2.zip 574395233 md5:3f6ce4317e3a0091978ae5cb8fbf05a3` on the Zenodo API; local `md5 -q` = `3f6ce4317e3a0091978ae5cb8fbf05a3` (**match**).
- The zip had a nested folder `GenePT_emebdding_v2/`; I extracted flat with `unzip -j`. No renames: the file `GenePT_gene_protein_embedding_model_3_text.pickle.` keeps its trailing dot as in the zip.
- `df -h .` after extraction showed 23Gi free (later 20Gi), so the zip was **kept**.
- `data/README.md` written (table in §8). `.gitignore` contains `data/` and `.venv/`.

Zenodo record description (`https://zenodo.org/api/records/10833191`, fetched 2026-10-04T06:02:01Z), relevant to text/embedding pairing:
> These are the pulled NCBI (and UniProt, when applicable) summaries of genes, as well as the corresponding OpenAI text embeddings (text-embedding-ada-002 and text-embedding-3-large) computed on the summaries. See methods details in Chen and Zou (2024+).
> The unzipped folder contains four different files:
> NCBI_summary_of_genes.json (NCBI gene card summary of human genes)
> NCBI_UniProt_summary_of_genes.json (NCBI gene card and UniProt protein (when applicable) summary of human genes)
> GenePT_gene_embedding_ada_text.pickle (a dictionary of numpy array where gene names (upper case) are keys and text-embedding-ada-002 embeddings of the summary in 1. are the values)
> GenePT_gene_protein_embedding_model_3_text.pickle (a dictionary of numpy array where gene names (upper case) are keys and text-embedding-3-large embeddings of the summary in 1. are the values)

Record publication_date: 2024-03-18; license cc-by-4.0.

### 3.7 Step 6: embedding models
```
sentence-transformers/all-MiniLM-L6-v2 (1, 384)
BAAI/bge-small-en-v1.5 (1, 384)
```
Snapshots: all-MiniLM-L6-v2 `1110a243fdf4706b3f48f1d95db1a4f5529b4d41`; bge-small-en-v1.5 `5c38ec7c405ec4b44b94cc5a9bb96e735b38267a`. stderr held only an unauthenticated HF Hub warning and weight-loading progress bars.

### 3.8 Step 7: `handback/stage1_data_report.txt` (full, verbatim)
Decision table: the report ends with `END OF REPORT` and has 0 lines containing `FAILED` → go to Step 8.
```
== GENEPT FILES ==
Archive:  data/downloads/GenePT_emebdding_v2.zip
  Length      Date    Time    Name
---------  ---------- -----   ----
        0  03-18-2024 13:22   GenePT_emebdding_v2/
 19163456  02-19-2024 14:38   GenePT_emebdding_v2/NCBI_UniProt_summary_of_genes.json
460797248  02-15-2024 08:44   GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle
909175027  02-26-2024 14:42   GenePT_emebdding_v2/GenePT_gene_protein_embedding_model_3_text.pickle.
 10982188  02-19-2024 13:41   GenePT_emebdding_v2/NCBI_summary_of_genes.json
---------                     -------
1400117919                     5 files
GenePT_gene_embedding_ada_text.pickle	460797248 bytes
  pickle dict entries: 93800; vector dims (first 2000): {1536: 2000}
GenePT_gene_protein_embedding_model_3_text.pickle.	909175027 bytes
NCBI_UniProt_summary_of_genes.json	19163456 bytes
  json dict entries: 33703
NCBI_summary_of_genes.json	10982188 bytes
  json dict entries: 33703
zip_listing.txt	598 bytes

== TEXT / EMBEDDING KEY PAIRING ==
NCBI_UniProt_summary_of_genes.json vs GenePT_gene_embedding_ada_text.pickle: text=33703 emb=93800 both=33703 text_only=0 emb_only=60097
NCBI_summary_of_genes.json vs GenePT_gene_embedding_ada_text.pickle: text=33703 emb=93800 both=33703 text_only=0 emb_only=60097
NCBI_UniProt_summary_of_genes.json vs NCBI_summary_of_genes.json: common keys=33703 identical text=14361 NCBI_summary_of_genes.json longer=0

== SUMMARY TEXT DESCRIPTION ==
NCBI_UniProt_summary_of_genes.json: n=33703 words min=3 q25=3 median=57 q75=131 max=1410
  literal substring counts: {'Alliance of Genome Resources': 0, 'provided by RefSeq': 0, 'provided by': 54, 'Predicted to': 4544, 'Gene Ontology': 0, 'GO:': 0}
NCBI_summary_of_genes.json: n=33703 words min=3 q25=3 median=31 q75=74 max=518
  literal substring counts: {'Alliance of Genome Resources': 0, 'provided by RefSeq': 0, 'provided by': 0, 'Predicted to': 4541, 'Gene Ontology': 0, 'GO:': 0}

== ALLIANCE TEMPLATE CLASSIFICATION ==
NCBI_UniProt_summary_of_genes.json: alliance_template=375 (1.1%), curated=26882 (79.8%), mixed=6446 (19.1%)
  [alliance_template] PLGLB2: Gene Symbol PLGLB2 Predicted to enable serine-type endopeptidase activity. Predicted to be involved in proteolysis. Predicted to be located in extracellular region.
  [alliance_template] RGPD4: Gene Symbol RGPD4 Predicted to contribute to GTPase activator activity. Predicted to be involved in NLS-bearing protein import into nucleus. Predicted to be part of nuclear pore. Predicted to be active in cytoplasm. Protein summary: One of the 8 copies of RANBP2 clustered close to the chromosome 2 centromere.
  [alliance_template] BBX: Gene Symbol BBX Enables sequence-specific double-stranded DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to act upstream of or within bone development. Located in cytosol and nucleoplasm. Protein summary: Transcription factor that is necessary for cell cycle progression from G1 to S phase.
  [alliance_template] LMNTD2: Gene Symbol LMNTD2 Predicted to be a structural constituent of chromatin. Predicted to be involved in regulation of chromatin assembly. Predicted to act upstream of or within positive regulation of mRNA splicing, via spliceosome. Predicted to be active in lamin filament. Protein summary: Computationally mapped potential isoform sequences There are 4 potential isoforms mapped to this entry Align Ad...
  [alliance_template] SNRNP48: Gene Symbol SNRNP48 Predicted to enable metal ion binding activity. Predicted to be involved in RNA splicing. Located in cytosol and nucleoplasm. Part of U12-type spliceosomal complex. Protein summary: Likely involved in U12-type 5' splice site recognition.
  [mixed] RAB39A: Gene Symbol RAB39A Predicted to enable GTPase activity. Involved in phagosome acidification and phagosome-lysosome fusion. Located in phagocytic vesicle. Protein summary: Plays a role in the maturation and acidification of phagosomes that engulf pathogens, such as S.aureus and M.tuberculosis. Plays a role in vesicular trafficking. Plays a role in the fusion of phagosomes with lysosomes. Negatively...
  [mixed] NSF: Gene Symbol NSF Enables PDZ domain binding activity and ionotropic glutamate receptor binding activity. Involved in intracellular protein transport; positive regulation of protein catabolic process; and positive regulation of receptor recycling. Located in Golgi apparatus; cytosol; and plasma membrane. Implicated in developmental and epileptic encephalopathy. Protein summary: Required for vesicle-...
  [mixed] ZSCAN18: Gene Symbol ZSCAN18 Predicted to enable DNA-binding transcription factor activity, RNA polymerase II-specific and RNA polymerase II cis-regulatory region sequence-specific DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to be located in nucleus. Protein summary: May be involved in transcriptional regulation.
  [mixed] KLF18: Gene Symbol KLF18 Predicted to enable DNA-binding transcription factor activity, RNA polymerase II-specific and RNA polymerase II cis-regulatory region sequence-specific DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to be located in nucleus. Protein summary: Belongs to the krueppel C2H2-type zinc-finger protein family .
  [mixed] PSME4: Gene Symbol PSME4 Predicted to enable lysine-acetylated histone binding activity; peptidase activator activity; and proteasome binding activity. Predicted to be involved in DNA repair; proteasomal ubiquitin-independent protein catabolic process; and spermatogenesis, exchange of chromosomal proteins. Located in nucleoplasm. Protein summary: Associated component of the proteasome that specifically r...
  [curated] MIR3925: Gene Symbol MIR3925 microRNAs (miRNAs) are short (20-24 nt) non-coding RNAs that are involved in post-transcriptional regulation of gene expression in multicellular organisms by affecting both the stability and translation of mRNAs. miRNAs are transcribed by RNA polymerase II as part of capped and polyadenylated primary transcripts (pri-miRNAs) that can be either protein-coding or non-coding. The ...
  [curated] RNA5SP87: Gene Symbol RNA5SP87
  [curated] HPYR1: Gene Symbol HPYR1
  [curated] POU5F1P3: Gene Symbol POU5F1P3
  [curated] EFNB3: Gene Symbol EFNB3 EFNB3, a member of the ephrin gene family, is important in brain development as well as in its maintenance. Moreover, since levels of EFNB3 expression were particularly high in several forebrain subregions compared to other brain subregions, it may play a pivotal role in forebrain function. The EPH and EPH-related receptors comprise the largest subfamily of receptor protein-tyros...
NCBI_summary_of_genes.json: alliance_template=2578 (7.6%), curated=28436 (84.4%), mixed=2689 (8.0%)
  [alliance_template] ONECUT3: Gene Symbol ONECUT3 Enables sequence-specific double-stranded DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to be part of chromatin. Predicted to be active in nucleus.
  [alliance_template] FOLH1B: Gene Symbol FOLH1B Enables metallocarboxypeptidase activity. Involved in C-terminal protein deglutamylation. Located in extracellular exosome.
  [alliance_template] DDX52: Gene Symbol DDX52 Enables RNA binding activity. Predicted to be involved in maturation of SSU-rRNA. Located in membrane.
  [alliance_template] ZNF850: Gene Symbol ZNF850 Predicted to enable DNA-binding transcription activator activity, RNA polymerase II-specific and RNA polymerase II cis-regulatory region sequence-specific DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to be active in nucleus.
  [alliance_template] MTRNR2L10: Gene Symbol MTRNR2L10 Predicted to enable receptor antagonist activity. Predicted to be involved in negative regulation of execution phase of apoptosis. Predicted to be located in cytoplasm and extracellular region.
  [mixed] SYCP1: Gene Symbol SYCP1 Enables double-stranded DNA binding activity. Involved in protein homotetramerization. Predicted to be located in synaptonemal complex. Predicted to be active in central element; male germ cell nucleus; and transverse filament.
  [mixed] TTPAL: Gene Symbol TTPAL Predicted to enable phosphatidylinositol bisphosphate binding activity. Located in membrane.
  [mixed] EPHA6: Gene Symbol EPHA6 Predicted to enable transmembrane-ephrin receptor activity. Predicted to be involved in axon guidance; positive regulation of kinase activity; and transmembrane receptor protein tyrosine kinase signaling pathway. Located in nucleoplasm.
  [mixed] MID1IP1: Gene Symbol MID1IP1 Predicted to enable identical protein binding activity and protein C-terminus binding activity. Predicted to be involved in several processes, including negative regulation of microtubule depolymerization; positive regulation of fatty acid biosynthetic process; and protein polymerization. Predicted to be located in cytoplasm and microtubule cytoskeleton. Predicted to be active ...
  [mixed] COX7B2: Gene Symbol COX7B2 Predicted to enable cytochrome-c oxidase activity. Predicted to be involved in electron transport chain; oxidative phosphorylation; and proton transmembrane transport. Predicted to be located in mitochondrial respirasome. Predicted to be integral component of membrane. Predicted to be part of respiratory chain complex IV.
  [curated] SH2D3C: Gene Symbol SH2D3C This gene encodes an adaptor protein and member of a cytoplasmic protein family involved in cell migration. The encoded protein contains a putative Src homology 2 (SH2) domain and guanine nucleotide exchange factor-like domain which allows this signaling protein to form a complex with scaffolding protein Crk-associated substrate. Multiple transcript variants encoding different i...
  [curated] CCDC200: Gene Symbol CCDC200
  [curated] WDR64: Gene Symbol WDR64
  [curated] RPL21P65: Gene Symbol RPL21P65
  [curated] LINC02474: Gene Symbol LINC02474

== HGNC AND ID OVERLAP ==
HGNC symbols=45187 protein-coding=19299 prev/alias entries=58844
gene2vec first line tokens=201 (dim=200), genes=24447
text:NCBI_UniProt_summary_of_genes.json: keys=33703 protein-coding direct=18754 (97.2% of HGNC pc) rescued via prev/alias=157 not_in_HGNC=353
text:NCBI_summary_of_genes.json: keys=33703 protein-coding direct=18754 (97.2% of HGNC pc) rescued via prev/alias=157 not_in_HGNC=353
emb:GenePT_gene_embedding_ada_text.pickle: keys=93800 protein-coding direct=18821 (97.5% of HGNC pc) rescued via prev/alias=14516 not_in_HGNC=60324
gene2vec: keys=24447 protein-coding direct=16617 (86.1% of HGNC pc) rescued via prev/alias=1250 not_in_HGNC=5558
protein-coding genes present in ALL sources (direct match): 16532

== GO RELEASES ==
baseline: GAF header ['!gaf-version: 2.2', '!generated-by: GOC', '!date-generated: 2024-04-03T21:10', '!generated-by: GOC', '!date-generated: 2024-03-29T22:12', '!generated-by: GOC', '!date-generated: 2024-03-29T10:27', '!date-generated: 2024-02-09 09:07', '!generated-by: UniProt', '!go-version: http://purl.obolibrary.org/obo/go/releases/2024-01-28/extensions/go-plus.owl', '!generated-by: PANTHER', '!date-generated: 2024-02-13']
  rows=640678 BP rows=162073 BP NOT-qualified=590 experimental BP (gene,term) pairs, unpropagated=45016 genes=10032
  OBO ['format-version: 1.2', 'data-version: releases/2024-03-28']; terms=47817 BP terms=30709 obsolete=5505
current: GAF header ['!gaf-version: 2.2', '!generated-by: GOC', '!date-generated: 2026-05-21T08:10', '!generated-by: GOC', '!date-generated: 2026-05-20T11:45', '!generated-by: UniProt', '!date-generated: 2026-04-30 09:29', '!go-version: http://purl.obolibrary.org/obo/go/releases/2026-04-27/extensions/go-plus.ofn', '!generated-by: PANTHER', '!date-generated: 2026-01-05']
  rows=906445 BP rows=204919 BP NOT-qualified=621 experimental BP (gene,term) pairs, unpropagated=50751 genes=11629
  OBO ['format-version: 1.2', 'data-version: releases/2026-07-26']; terms=48340 BP terms=30885 obsolete=10248
experimental BP pairs in current but not baseline (unpropagated)=9887; in baseline but not current=4152
terms gaining >=10 new genes (unpropagated): 186

== GENE2PUBMED ==
human GeneIDs with >=1 PubMed link=163238; median papers=1; max=20402

END OF REPORT
```

Extra descriptive check (not requested; §5 D6). The report skips the model_3 pickle because its name ends in `.` rather than `.pickle`. I described it with the report's own `load_pickle_embeddings`, with no change to the script (`handback/logs/extra_model3_pickle_describe.txt`):
```
data/raw/genept/GenePT_gene_protein_embedding_model_3_text.pickle.: pickle dict entries: 133736; vector dims (first 2000): {3072: 2000}
NCBI_summary_of_genes.json vs model_3 pickle: text=33703 emb=133736 both=33703 text_only=0 emb_only=100033
NCBI_UniProt_summary_of_genes.json vs model_3 pickle: text=33703 emb=133736 both=33703 text_only=0 emb_only=100033
ada pickle vs model_3 pickle: ada=93800 model3=133736 both=93800 ada_only=0 model3_only=39936
```

### 3.9 Step 8: NCBI attribution spot-check
- **Gene choice:** the report prints samples for two JSON files. The cap is 9 genes, so I took "the first 3 genes of each class" in report order, which come from the first file printed, `NCBI_UniProt_summary_of_genes.json` (§7 Q6).
- **Source:** the NCBI web URL (`https://www.ncbi.nlm.nih.gov/gene/?term=PLGLB2[sym]+AND+human[orgn]`, 2026-10-04T05:59:31Z) returned "Checking your browser - reCAPTCHA". So I read the same Gene records through NCBI E-utilities (`esearch` with `<SYM>[sym] AND human[orgn]`, then `efetch db=gene retmode=xml`, field `<Entrezgene_summary>`), fetched 2026-10-04T06:00:38Z–06:00:58Z.

Verbatim NCBI summaries:
- **PLGLB2** (GeneID 5342): "Predicted to be located in extracellular region. [provided by Alliance of Genome Resources, Jul 2025]"
- **RGPD4** (GeneID 285190): "Predicted to enable SUMO transferase activity. Predicted to be involved in NLS-bearing protein import into nucleus and nuclear export. Predicted to be part of nuclear pore. [provided by Alliance of Genome Resources, Jun 2026]"
- **BBX** (GeneID 56987): "Enables sequence-specific double-stranded DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to act upstream of or within bone development. Located in nucleoplasm. [provided by Alliance of Genome Resources, Jun 2026]"
- **RAB39A** (GeneID 54734): "Enables G protein activity. Involved in autophagosome-lysosome fusion; phagosome acidification; and phagosome-lysosome fusion. Located in autolysosome membrane and phagocytic vesicle. [provided by Alliance of Genome Resources, Jun 2026]"
- **NSF** (GeneID 4905): "Enables PDZ domain binding activity and ionotropic glutamate receptor binding activity. Involved in intracellular protein transport; positive regulation of protein catabolic process; and positive regulation of receptor recycling. Located in cytosol and plasma membrane. Implicated in cocaine dependence and developmental and epileptic encephalopathy 96. [provided by Alliance of Genome Resources, Jun 2026]"
- **ZSCAN18** (GeneID 65982): "Predicted to enable DNA-binding transcription factor activity, RNA polymerase II-specific and RNA polymerase II cis-regulatory region sequence-specific DNA binding activity. Predicted to be involved in regulation of transcription by RNA polymerase II. Predicted to be located in nucleus. [provided by Alliance of Genome Resources, Jul 2025]"
- **MIR3925** (GeneID 100500885): "microRNAs (miRNAs) are short (20-24 nt) non-coding RNAs that are involved in post-transcriptional regulation of gene expression in multicellular organisms by affecting both the stability and translation of mRNAs. miRNAs are transcribed by RNA polymerase II as part of capped and polyadenylated primary transcripts (pri-miRNAs) that can be either protein-coding or non-coding. The primary transcript is cleaved by the Drosha ribonuclease III enzyme to produce an approximately 70-nt stem-loop precursor miRNA (pre-miRNA), which is further cleaved by the cytoplasmic Dicer ribonuclease to generate the mature miRNA and antisense miRNA star (miRNA*) products. The mature miRNA is incorporated into a RNA-induced silencing complex (RISC), which recognizes target mRNAs through imperfect base pairing with the miRNA and most commonly results in translational inhibition or destabilization of the target mRNA. The RefSeq represents the predicted microRNA stem-loop. [provided by RefSeq, Sep 2009]"
- **RNA5SP87** (GeneID 100873320): no `<Entrezgene_summary>` element (no summary).
- **HPYR1** (GeneID 93668): no `<Entrezgene_summary>` element (no summary).

| Gene | Our class (NCBI_UniProt file) | NCBI attribution (current record) |
|---|---|---|
| PLGLB2 | alliance_template | [provided by Alliance of Genome Resources, Jul 2025] |
| RGPD4 | alliance_template | [provided by Alliance of Genome Resources, Jun 2026] |
| BBX | alliance_template | [provided by Alliance of Genome Resources, Jun 2026] |
| RAB39A | mixed | [provided by Alliance of Genome Resources, Jun 2026] |
| NSF | mixed | [provided by Alliance of Genome Resources, Jun 2026] |
| ZSCAN18 | mixed | [provided by Alliance of Genome Resources, Jul 2025] |
| MIR3925 | curated | [provided by RefSeq, Sep 2009] |
| RNA5SP87 | curated | no summary on NCBI |
| HPYR1 | curated | no summary on NCBI |

### 3.10 Final disk check
```
Filesystem      Size    Used   Avail Capacity iused ifree %iused  Mounted on
/dev/disk3s5   460Gi   402Gi    20Gi    96%    3.2M  206M    2%   /System/Volumes/Data
```

## 4. Code changes
**none.** `diff -ru -x __pycache__ handback/stage1_code_before code` printed nothing (rc=0). Without `-x`, the only output is `Only in code: __pycache__` and `Only in code/tests: __pycache__`, created by pytest. The SHA-256 of both code files is unchanged (`b6bc90ad…`, `72e5612b…`).

## 5. Deviations
- **D1. Step 1 stop and resume.** The first disk check showed 2.9Gi (< 5 GB), so I stopped and reported, as the prompt says. The user freed space and said "go". I re-ran `df -h .` (28Gi) and continued from Step 2. The code snapshot from the first attempt was kept.
- **D2. Step 4b source.** PMC served a reCAPTCHA and academic.oup.com returned 403, so I quoted Kishore et al. from the Europe PMC full-text XML of the same PMC record (PMC7304461).
- **D3. Step 8 source.** NCBI Gene web pages served a reCAPTCHA, so I read the same Gene records' `Entrezgene_summary` via NCBI E-utilities. The text is NCBI's own record, including the bracketed attribution.
- **D4. PDF text extraction** used `pypdf` 6.16.2 from the system Python (`/opt/anaconda3`), not the project venv, so the venv matches `requirements.txt`. One page was rendered with macOS `sips` to read the Fig. B3 legend.
- **D5. Extra pages fetched beyond the prompt's list** (read-only): the IDASB flyer PDF (linked from the IDASB CFP page); BIBM CyberChair `ws_submit.php` and the IDASB S44 page (linked from the BIBM menu, but on wi-lab.com); the Zenodo record API (for the MD5 check and the file description).
- **D6. Extra descriptive check** of the model_3 pickle (§3.8): key count, vector dimension and key overlap only, using the unchanged script's loader. No vectors were compared and nothing was computed from them.
- **D7. Steps 3, 4 and 5 ran in parallel.** Downloads ran in the background while I fetched pages.

## 6. Problems and uncertainties
- **P1. Venue deadline conflict.** The IDASB website and the BIBM workshop list say Sept 27, 2026; the IDASB proposal flyer says October 15, 2026; CyberChair says "to be announced soon" with an active "Submit a New Paper" link. Which one is in force is not determinable from these pages. Not stated anywhere for IDASB: review type, presentation format, extra-page fee. The BIBM main CFP's "double blind" sentence is on the main-conference page. The BIBM workshop proposal call lets organizers choose "Full onsite, full online, or hybrid".
- **P2. bioRxiv HTTP 429.**
  - The GenePT v2 supplementary-material page was never fetched: 13 attempts, all HTTP 429, with Retry-After 57–109 s where reported; log in §3.5. Its 17-byte error body (`error code: 1015`) is kept as `data/papers/genept_biorxiv_v2_supplementary-material.HTTP429-body.txt`. So I cannot say whether v2 has supplementary files, but B.3 is in the main PDF.
  - Zhong v2 also returned 429 at first; it succeeded at 06:08:31Z.
  - The GenePT v2 PDF fetched fine before the rate limit started.
- **P3. The model_3 pickle has a trailing dot** (`...model_3_text.pickle.`) inside the zip itself. The report's `.pickle` check therefore skips it silently. This is not a `FAILED` line, so the Step 7 table did not trigger.
- **P4. Pairing is not resolvable from keys.**
  - Both JSONs have the same 33,703 keys.
  - Both pickles contain all of them: ada 93,800 keys (1536-d); model_3 133,736 keys (3072-d).
  - Zenodo says both pickles embed "the summary in 1." (`NCBI_summary_of_genes.json`), but the model_3 file name says `gene_protein`.
  - The two JSONs share 33,703 keys with 14,361 identical texts, and the NCBI-only text is never longer.
  - The pickles have 60,097 (ada) and 100,033 (model_3) keys with no summary in either JSON; their origin is unexplained.
- **P5. Every summary begins with "Gene Symbol <SYM> "**, glued to the first sentence, e.g. `Gene Symbol ONECUT3 Enables sequence-specific ...`.
  - `classify_summary` matches openings only at the sentence start, so the first template sentence of every summary is never counted. A summary of k template sentences gets at most k−1 hits.
  - In the UniProt file, a trailing "Protein summary: ..." block adds non-template sentences.
  - Both effects lower the `alliance_template` share and raise `mixed` / `curated`. I did not change the code; the code matches its written spec ("Matched at the start of a sentence").
- **P6. Many near-empty summaries.** Both JSONs have `q25=3` words. At least a quarter of entries look like `Gene Symbol <SYM>` with no summary (e.g. RNA5SP87, HPYR1, CCDC200, WDR64). I did not count them exactly, since that would be a new computation not in the prompt.
- **P7. Baseline release vs summary snapshot.** The rule gave 2024-03-28.
  - The JSONs inside the zip are dated 02-19-2024, and GenePT v2 was posted March 5, 2024.
  - The 2024-03-28 GAF header shows `date-generated: 2024-03-29` / `2024-04-03` and `go-version ... releases/2024-01-28`.
  - So the baseline is about 5–6 weeks after the summary file timestamps. 2024-01-17 is the last listed release before them.
  - STAGES.md describes the baseline as "one before the summary snapshot". I followed the prompt's fixed rule.
- **P8. "Current" GO files are not from the same date.** The current GAF header has `go-version .../releases/2026-04-27/...` and `date-generated: 2026-05-21T08:10`, while the current OBO has `data-version: releases/2026-07-26`.
- **P9. The spot-check reflects today's NCBI, not the 2024 snapshot.** Attribution dates are Jul 2025 / Jun 2026, and some texts changed. For example, the GenePT text for RGPD4 says "Predicted to contribute to GTPase activator activity", while NCBI now says "Predicted to enable SUMO transferase activity". The attribution shows these genes' current summaries come from the Alliance; it does not directly prove the 2024 text's source, though the 2024 text has the same template form.
- **P10. Template wording drift.** Kishore et al. 2020 describe MF sentences as "Exhibits ...", while current NCBI text uses "Enables ...". The script's opening list contains both.
- **P11. ID overlap.** Protein-coding genes present in all sources (direct match): 16532, against HGNC protein-coding = 19299. Gene2vec is the limiting source: 16617 direct (86.1% of HGNC pc).
- **P12. Shell mistakes during the run, all caught and redone:**
  - zsh treated `$g[sym]` as an array subscript in the first Step 8 loop, so all 9 genes returned the same EGFR record. I deleted those files and re-ran with `${g}`; only the corrected run is reported.
  - `echo ======` failed in zsh.
  - One curl URL with `[...]` needed `-g`.
- **P13. Free disk varied** (28Gi → 23Gi → 20Gi) by more than this stage's ~3.5 GB of files. Other processes on the Mac may be using space.

## 7. Questions for the reviewer
- **Q1. Venue.** Which IDASB deadline applies: Sept 27 (website and BIBM list) or Oct 15 (proposal flyer), with CyberChair "to be announced soon" and the link open? *Recommended:* the author emails the IDASB contact listed by BIBM (h.zheng@ulster.ac.uk), subject including "IDASB 2026", asking for the paper deadline, presentation format and review type. Name a fallback venue in the same Stage 2 prompt in case the answer is "closed".
- **Q2. Text/embedding pairing.** *Recommended:*
  - Pair `GenePT_gene_embedding_ada_text.pickle` with `NCBI_summary_of_genes.json` for C0-ada. This follows Zenodo ("embeddings of the summary in 1.") and GenePT methods (ada on NCBI summaries).
  - Leave the model_3 pickle out of the study unless needed, because its file name contradicts Zenodo's description.
  - Treat the extra pickle keys (no summary text) as outside the gene universe.
- **Q3. Template classifier and the "Gene Symbol <SYM> " prefix.** Should the Stage 2 spec strip the `Gene Symbol <SYM> ` prefix (and, for the UniProt file, everything from "Protein summary:" on) before sentence splitting? Then the template shares would be re-reported. *Recommended:* yes. Strip the prefix in a spec change you write, and report the NCBI-only file as primary for the curated vs Alliance split.
- **Q4. GO baseline.** Keep 2024-03-28 (prompt rule; after the JSON timestamps of 2024-02-19), or switch to 2024-01-17 (the last release before them)? *Recommended:* switch to 2024-01-17 if the baseline must predate the summary snapshot, as STAGES.md says; otherwise keep 2024-03-28 and state the overlap in the paper.
- **Q5. The model_3 file name with a trailing dot.** *Recommended:* never rename it (it matches the zip and its SHA-256 is recorded); have any Stage 2 loader take explicit file paths rather than matching extensions.
- **Q6. Step 8 sample choice.** I used the first file printed (NCBI_UniProt) to stay within 9 genes. *Recommended:* if you want the NCBI-only samples too (ONECUT3, FOLH1B, DDX52 / SYCP1, TTPAL, EPHA6 / SH2D3C, CCDC200, WDR64), add them to the Stage 2 prompt.
- **Q7. GenePT v2 supplementary-material page (blocked by bioRxiv 429).** *Recommended:* the author opens `https://www.biorxiv.org/content/10.1101/2023.10.16.562533v2.supplementary-material` in a browser and says in the thread whether any files are listed. Since B.3 was found in the main PDF, this can wait and need not block Stage 2.
- **Q8. "Current" GO pair.** *Recommended:* pin the current GAF and OBO to one dated folder on release.geneontology.org (e.g. the folder matching the GAF's go-version) for reproducibility, rather than `current.geneontology.org`.

## 8. Files

Data files (also in `data/README.md`). Sizes are in bytes. All downloaded 2026-10-04.

| Path | Size | SHA-256 |
|---|---|---|
| `data/downloads/GenePT_emebdding_v2.zip` | 574395233 | `6193575dcbd7bbf214c8ca3eb518cb3d13272443f98ecd6c26402411a7ca745e` |
| `data/raw/genept/GenePT_gene_embedding_ada_text.pickle` | 460797248 | `fd297510ddd3040744033fde0b0f2cf15a40ac8b2fd2fb02f10667295e55c862` |
| `data/raw/genept/GenePT_gene_protein_embedding_model_3_text.pickle.` | 909175027 | `cd64611edab310b517cbb461295605a37821f6cd4f8ddae90ba18a2f861a6618` |
| `data/raw/genept/NCBI_summary_of_genes.json` | 10982188 | `3db721b7cfbc35795428c6a7ab4e8a9c59fb5f5ccb8c64af50dd3c97c2d4de2e` |
| `data/raw/genept/NCBI_UniProt_summary_of_genes.json` | 19163456 | `319897d7e2f8a0306bccdf2d022982b3c8321dbec8273739094021d57e414c31` |
| `data/raw/genept/zip_listing.txt` | 598 | `973ff705ca344ddb45e886923f971b49299b0b7b3244fed301d067c76d8ecfd2` |
| `data/raw/go/current/go-basic.obo` | 32227785 | `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc` |
| `data/raw/go/current/goa_human.gaf.gz` | 15041996 | `db472faff1785878521693af62646546cea6af4a386609dd01c72c0554a46a30` |
| `data/raw/go/baseline/go-basic.obo` | 31267219 | `036ac170c4cc38e35af1749e208b1066b8a9e42edec4cc3aa5ed8befddb1e672` |
| `data/raw/go/baseline/goa_human.gaf.gz` | 11886126 | `a561c363d3f15caf173f3cf6101988b6dbed088de645004746b2d9f3f1f3c264` |
| `data/raw/hgnc_complete_set.txt` | 16963116 | `2b4224ea847df2fc6982f5b2a52804c5d92fbb8810b134afb636f6029452dc03` |
| `data/raw/Homo_sapiens.gene_info.gz` | 5191996 | `cdeb1bfc6c43e6d375ccc40ee56d0a08e7befff60dbc63b364e2b1bb825e3727` |
| `data/raw/gene2pubmed_human.tsv` | 65068222 | `e5458118101a0bc103e49ef9953690c03282945f9f1ca6b9f3c2fca205773001` |
| `data/raw/gene2vec_dim_200_iter_9.txt` | 56682049 | `3441d642d4b153e4366dac02d930b3d941b3002d5484a8b3b1b85974d4feeb61` |

Other new files:

| Path | Size (bytes) | SHA-256 |
|---|---|---|
| `CLAUDE.md` | 3218 | `d5bc164b35659afcf060f3a1150c0721745c3eb2095944c3b34e693d4b3242fa` |
| `STAGES.md` | 7501 | `1deb4621b9dfb4daaafc11ad174ce80328c8a4f4f6069e5c95353294bc2aa3ed` |
| `review.md` | 7157 | `cfa6bac7509212c6e7a6dd1dc7762eaf125cdb54d70fe416a8c2f1624514bc14` |
| `plan.md` | 21461 | `6f541e6904a2a00903da1846a5f6bfc7e8ba5347f836f439b900f868708605ad` |
| `requirements.txt` | 67 | `8028ae8eba9c8f6134ed73c78f7769849c9ea653b0d76bc38e34d53e35f5859e` |
| `.gitignore` | 13 | `0e8948407a9283311a46fde471fc663495d18e6238875b246301c00c3288e21d` |
| `code/stage1_data_report.py` | 12402 | `b6bc90ad80723450b53e9fcbf99f2ddffb64a6a80d653883b5de1a26f5ba45d3` |
| `code/tests/test_stage1_report.py` | 3678 | `72e5612b2f991ce9d12343004395f1e340998a5f08a91e756223c08f62c41fe9` |
| `data/README.md` | 3877 | `9a19206d18caa1dfddbb3ba62c0bead4c3c4eb53f4a59490c015ec11af366898` |
| `handback/stage1_data_report.txt` | 12347 | `2380e839de1901fdd0159d74a454ac95a51c1b95cffce37e92cad16e6c88321c` |
| `handback/stage1_pip_freeze.txt` | 789 | `5c59fdf5ea5cd8ed80aa7aee9bebe8bdda2a40c8ddde780fb7a1d9a8eb235a9d` |
| `handback/stage1_pytest.txt` | 98 | `d161a71b963882cba02d9bdbbed80589c4c38da80a1833e2bfba476f94a0aebb` |
| `handback/stage1_step6_models.txt` | 80 | `8744326aed9442ec59ddebe220b5333f791be177bb63881a0d79f5a90561579e` |
| `handback/logs/biorxiv_retry.log` | 671 | `e707aa48fcdd5574de1dfae35a9f128167683812209a5b98a8a77be15b9e32e3` |
| `handback/logs/biorxiv_retry2.sh` | 1706 | `db56361116b757a9b843737cda06207ae88a04df900d61984c9b77b5f3fa320d` |
| `handback/logs/data_checksums.tsv` | 1595 | `5340c32ac75fa662e320a3f82e7611cdb4ce398ea55b86708d2ea9af6bdb3984` |
| `handback/logs/download_baseline.log` | 836 | `3ed3923f14c8e8095b8648dbc3160b96ac39f6ad4b3ab79a26abff64c9267c96` |
| `handback/logs/download_main.log` | 20902 | `73b5021f148b91dccb11fcdaf1f9f00313cccc3ffc4f10ee6a0e61d14a789109` |
| `handback/logs/download_main.sh` | 1285 | `4446e3e72ba33224d9a2e4ab2a4e93c9baa043b8bba5c0958b30cb3daec4039c` |
| `handback/logs/extra_model3_pickle_describe.txt` | 452 | `cd280248e48eea0261d5618fadbc87fb4d1f5038724e0f6a4f317b94a854c063` |
| `handback/logs/go_release_listing.html` | 26290 | `13bd2354742beec0d5ef4e7dc365b31587ff24e0bcf1d541455d27238e87479e` |
| `handback/logs/html2text.py` | 1106 | `77b840f6cc2993e7bd915617f32da3180539ea8a8bd4be340de74de2cb130ba0` |
| `handback/logs/list_files.sh` | 932 | `cfaf6d677430cb68ba1d12e0c1a193e815a2f4cf36a260d876b21dec4b86439f` |
| `handback/logs/step6_stderr.txt` | 434 | `a3c38c2d608115eb79718c23006b291b6e8d748f0d4772436d547027ce1e6683` |
| `handback/logs/step8_ncbi_eutils.txt` | 3399 | `03944b9fc4b7144e5ae699955a6f9ef41f321e1f9eec5b91ecbb77c817a4c63e` |
| `data/papers/genept_biorxiv_v2_supplementary-material.HTTP429-body.txt` | 17 | `e395ad55dd18b1bf23944d00b156714accfae4ec834df208053575542dfae602` |
| `data/papers/genept_biorxiv_v2.full.pdf` | 15940010 | `cd8521a54746f50d33e726578b8ef3ffb48118edb973c304aeee819b0d9c5529` |
| `data/papers/genept_biorxiv_v2.full.txt` | 84309 | `23221f6698baeae8cc2bf92fd3a0ed9bee47b225bd8c27f5bc0120506d536d79` |
| `data/papers/kishore2020_europepmc_meta.json` | 10071 | `5226e2cb6ccbc05f375a35aac9bda7db94be660ed493013d491c6723029c861e` |
| `data/papers/kishore2020_oup.html` | 5843 | `e279a11d6fd394189133ae538a8a8d5c5fdb69ab08fae36ded1f016ee32d2712` |
| `data/papers/kishore2020_PMC7304461_fullText.txt` | 38890 | `a78084ffb9001c0eb45795a0f61c65f67eb9da7c32bd2622dba3b0a7ec366bc7` |
| `data/papers/kishore2020_PMC7304461_fullText.xml` | 121333 | `fbf4665ef20a410f3bd82be433077ce6679383dd4513c534f62a0ffd7df864ee` |
| `data/papers/kishore2020_PMC7304461.html` | 20411 | `5ebc854d37a1dfffa84d911486582185bc5d739ef4a3e7080c854ae5dd2099bd` |
| `data/papers/kishore2020_PMC7304461.txt` | 177 | `cb711a3c092be5873ed300f96cd13f276dc0904f1ef5954e5b362350c34655c8` |
| `data/papers/zhong2025_biorxiv_v2.full.html` | 285825 | `2a7a0f6cb3ddb67165f53e8083c46e263dcab6ee35da5e1452f38ea1edcb0b2a` |
| `data/papers/zhong2025_biorxiv_v2.full.txt` | 103336 | `3bebad0b279e13003ae8ff37b43f07d38b170f155ed6d7573cd012a411a27b4b` |
| `data/web/bibm_BIBM2026-registration0.html` | 7141 | `bfbc8fcc690f79ad0ec81514559714bd0fcc98d55b074ace05fa2a19a66c79a6` |
| `data/web/bibm_callforpaper2026.html` | 10469 | `2f39aae787dbd80c9f2023cf81dac9b233640e0823f76789ff0b03a7718b5ca0` |
| `data/web/bibm_CallforWorkshopProposals.html` | 13000 | `6cca0e31c6c5bea3507e1f681653f3d47a65f8c5482b1a79ce1027f4f98fd3f3` |
| `data/web/bibm_home.html` | 7770 | `62b8629f26337ac188dc7a4ea4ece4357795ef7e911d3069db94bd06a41d98fb` |
| `data/web/bibm_index.html` | 1545 | `134ea86a93403ec7bbd4dbdf10f883202e987ce15cc9ef57cfc6700936110436` |
| `data/web/bibm_menu.html` | 14297 | `19364523e9a56a67ef03715ed5814a25d69ce95ff2fc86e5b9200ce64662b819` |
| `data/web/bibm_workshopwebpges.html` | 65888 | `2530d92b545185d4c93073619b7def46d560f9aa3ac3dca7bc7d5248e654e131` |
| `data/web/bibm_workshopwebpges.txt` | 27918 | `1c7a76351ac3d987399450babbe7dfa323ab0bf998a13c84e1a9513837e4b0cf` |
| `data/web/bibm_ws_S44.html` | 5551 | `6e4a326cc2d3f0d54e6f5e2d3f767418a127c996fdf7644530011b63aef36f54` |
| `data/web/bibm_ws_submit.html` | 17541 | `973b8bfe0132042f8fb739a38979a666c1080d692b14394fe92c11b08ca86929` |
| `data/web/bibm_ws_submit.txt` | 12150 | `28a0b602f689ff9aedf0fdc07369078b00be28a69a6790165e7ea9f37c6f0467` |
| `data/web/idasb_app.5643eed8.js` | 66399 | `5fb5e14546779cc1fba40c2c0807b9b8d7ff1feeb2b2e9d0dba80ca4b5c6ed21` |
| `data/web/idasb_app.5643eed8.js.strings.txt` | 8520 | `3a30ac501331ba00adb507cee56f2f2e982d171d6fb872924a83de65840aab0e` |
| `data/web/idasb_chunk-vendors.e26e98df.js` | 852563 | `75d6ea42811cfb88d24e814b145f71ef2f9d37fa538c49ccc999f10e11ff1178` |
| `data/web/idasb_demo.9ef88cf9.js` | 617919 | `2e2bdfe7ea799e2c37fed34c47ed1907b9d2ec28289c97ba4e9f9a82361f8978` |
| `data/web/idasb_demo.9ef88cf9.js.strings.txt` | 112949 | `6f5277419f198af67231793c719c7b4bae966e9b773bfc16ef77501820d75a85` |
| `data/web/idasb_IDASB_2026.pdf` | 105034 | `55d8b83b656ac5dbe7407d1dd0c603d0cfe7130874dabda21729bc1cea014ba8` |
| `data/web/idasb_index.html` | 1553 | `e156ffa9f153ea0d33c839ee16b956dfbd9cae0d784918d628ad5239d8cf4ddb` |
| `data/web/idasb_textnodes.txt` | 13486 | `f2c07b89e686b695a7b97a34eb145a6789c2a6b585fc187e4dad6153b446087a` |
| `data/web/ncbi/efetch_BBX.xml` | 2522475 | `2f7108ba5ede7d24fbc76a37065071e1ecf8848dbf6aac16b4d8b30ebb13c598` |
| `data/web/ncbi/efetch_HPYR1.xml` | 70863 | `382aa6669fc49cea3bcd074b29430733379bae6d3160fe21f2f54d121053e2a1` |
| `data/web/ncbi/efetch_MIR3925.xml` | 74929 | `a58dfec927abc03e5164b37570c7f65038df69ffccd00f7efd5993866b7892cc` |
| `data/web/ncbi/efetch_NSF.xml` | 1661427 | `9cf383cc743818924254ee311c0a7f2d23c1e6ab4943735df883154836507871` |
| `data/web/ncbi/efetch_PLGLB2.xml` | 130940 | `5048e735f40560ba120cb4177b395328d43218002b9847e0f259278f89b6e3f2` |
| `data/web/ncbi/efetch_RAB39A.xml` | 362370 | `710c5f13056ab316d0fc42a562017b3e487bef1a651a1bea11a59360c2b17be6` |
| `data/web/ncbi/efetch_RGPD4.xml` | 767287 | `3dbf14a031cbe93a3cd8b4a70acfc56453722d2f92d13d68a28dde8627b7f9c8` |
| `data/web/ncbi/efetch_RNA5SP87.xml` | 59804 | `4b90451017b0cc2b49522a1547912fb7386db4ef01655c494166d645fa68d5dd` |
| `data/web/ncbi/efetch_ZSCAN18.xml` | 754427 | `3d719c548d5d7ac01b588accfccf3739c6bdade0951dd18aa7576a346e3f3458` |
| `data/web/ncbi/esearch_BBX.json` | 438 | `4831d408783760c799a872a2851d5335fe53bc406ce59df0d68d816bae78fbf8` |
| `data/web/ncbi/esearch_HPYR1.json` | 439 | `476e5764365849cf53cd797b8cb6fc6db37205fb815a5896b741a7b60583949e` |
| `data/web/ncbi/esearch_MIR3925.json` | 447 | `57e4e8842e64e4df6eb9b96ee9ab474b5fb6a1c90295c024770595be02510e3f` |
| `data/web/ncbi/esearch_NSF.json` | 436 | `f7b54f09f5a858161ccbdd183219062ec4a012b76d60374ee0d19f77a81ee7a5` |
| `data/web/ncbi/esearch_PLGLB2.json` | 440 | `c573cce9278d1a69cc7fd22ce551459a8d886531c7414102adaa148a0466bf3a` |
| `data/web/ncbi/esearch_RAB39A.json` | 443 | `b7b807847b0badde337dd368f0eadaae675d15996805ab490b656f2ce9ec77c0` |
| `data/web/ncbi/esearch_RGPD4.json` | 442 | `648e9c70ed209f74afc162a80073f869b6a7b71611b0d0edd2f14edfbf2f603d` |
| `data/web/ncbi/esearch_RNA5SP87.json` | 449 | `9d8218e2c673b16a65d86f7e2716f037844e49efd653478b3ff5b02861ed47dc` |
| `data/web/ncbi/esearch_ZSCAN18.json` | 445 | `1349a97670ec8695945071888273d7cd133b03bc620035344f51f0d6c1803ceb` |
| `data/web/ncbi/PLGLB2.html` | 20620 | `4200a81b54efe32790205b93ae7c8713bb021ee18bc007e95c2cdb0dadfc48c1` |
| `data/web/zenodo_10833191.json` | 5295 | `ba1359b40a6f86e8bc34e75b2d7cfb6fb9f841acb5bcfc469f42bba48e16df1b` |

Directories: `handback/stage1_code_before/` (copy of `code/` before the stage), `data/models/` (220972 KiB, Hugging Face cache), `.venv/` (1166124 KiB).
