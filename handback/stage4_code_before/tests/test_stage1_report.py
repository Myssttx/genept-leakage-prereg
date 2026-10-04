import gzip
import json
import os
import pickle
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from stage1_data_report import classify_summary  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(HERE), "stage1_data_report.py")


def test_classify_alliance_template():
    t = ("Predicted to enable DNA binding activity. Involved in DNA repair. "
         "Located in nucleus.")
    assert classify_summary(t)[0] == "alliance_template"


def test_classify_curated():
    t = ("This gene encodes a tumor suppressor protein containing transcriptional activation, "
         "DNA binding, and oligomerization domains. Mutations are associated with cancers.")
    assert classify_summary(t)[0] == "curated"


def test_classify_mixed():
    t = ("This gene encodes a kinase. It is expressed in brain. Involved in signal transduction.")
    assert classify_summary(t)[0] == "mixed"


def _gaf_line(sym, go, code, aspect="P", qual="involved_in"):
    r = ["UniProtKB", "P1", sym, qual, go, "PMID:1", code, "", aspect, "", "", "protein",
         "taxon:9606", "20240101", "GO_Central", "", ""]
    return "\t".join(r) + "\n"


def _make_fixture(root):
    raw = os.path.join(root, "data", "raw")
    os.makedirs(os.path.join(raw, "genept"))
    json.dump({"TP53": "This gene encodes a tumor suppressor.",
               "ABC1": "Predicted to enable ATP binding activity. Involved in transport. Located in membrane."},
              open(os.path.join(raw, "genept", "NCBI_summary_of_genes.json"), "w"))
    pickle.dump({"TP53": [0.1] * 8, "ABC1": [0.2] * 8},
                open(os.path.join(raw, "genept", "emb.pickle"), "wb"))
    open(os.path.join(raw, "genept", "zip_listing.txt"), "w").write("fake listing\n")
    with open(os.path.join(raw, "hgnc_complete_set.txt"), "w") as f:
        f.write("hgnc_id\tsymbol\tlocus_group\tprev_symbol\talias_symbol\n")
        f.write("HGNC:1\tTP53\tprotein-coding gene\t\t\n")
        f.write("HGNC:2\tABCA1\tprotein-coding gene\tABC1\t\n")
    with open(os.path.join(raw, "gene2vec_dim_200_iter_9.txt"), "w") as f:
        f.write("TP53 0.1 0.2 0.3 0.4\nABCA1 0.1 0.2 0.3 0.4\n")
    for rel, extra in [("baseline", []), ("current", [("ABCA1", "GO:0000002", "IDA")])]:
        d = os.path.join(raw, "go", rel)
        os.makedirs(d)
        with gzip.open(os.path.join(d, "goa_human.gaf.gz"), "wt") as f:
            f.write("!gaf-version: 2.2\n!date-generated: 2024-03-20\n")
            f.write(_gaf_line("TP53", "GO:0000001", "IMP"))
            f.write(_gaf_line("TP53", "GO:0000003", "IEA"))
            f.write(_gaf_line("TP53", "GO:0000004", "IDA", aspect="F"))
            for s, g, c in extra:
                f.write(_gaf_line(s, g, c))
        with open(os.path.join(d, "go-basic.obo"), "w") as f:
            f.write("format-version: 1.2\ndata-version: releases/2024-03-20\n\n[Term]\nid: GO:0000001\n"
                    "namespace: biological_process\n")
    with open(os.path.join(raw, "gene2pubmed_human.tsv"), "w") as f:
        f.write("#tax_id\tGeneID\tPubMed_ID\n9606\t7157\t1\n9606\t7157\t2\n9606\t19\t3\n")


def test_report_runs_on_fixture(tmp_path):
    _make_fixture(str(tmp_path))
    out = subprocess.run([sys.executable, SCRIPT, "--root", str(tmp_path)],
                         capture_output=True, text=True, check=True).stdout
    assert "FAILED" not in out, out
    assert "END OF REPORT" in out
    assert "alliance_template=1 (50.0%)" in out
    assert "rescued via prev/alias=1" in out
    assert "experimental BP pairs in current but not baseline (unpropagated)=1" in out
    assert "human GeneIDs with >=1 PubMed link=2" in out
