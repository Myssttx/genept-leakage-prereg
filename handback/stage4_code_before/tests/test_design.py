import gzip
import json
import os
import pickle
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import design_lib as D  # noqa: E402


# ---------------------------------------------------------------- text rules

def test_mask_keeps_prefix():
    text = "Gene Symbol REPAIR Involved in DNA repair. Repairs damage."
    out, n = D.apply_mask(text, frozenset({"repair", "dna"}))
    assert out.startswith("Gene Symbol REPAIR ")
    assert n == 3
    assert out == "Gene Symbol REPAIR Involved in . damage."


def test_mask_without_prefix_and_empty_body():
    assert D.split_prefix("Gene Symbol ABC1") == ("Gene Symbol ABC1", "")
    assert D.apply_mask("Gene Symbol ABC1", frozenset({"abc1"})) == ("Gene Symbol ABC1", 0)
    assert D.split_prefix("no prefix here") == ("", "no prefix here")


def test_plural_handling():
    assert D.norm("kinases") == "kinase"
    assert D.norm("Cells") == "cell"
    assert D.norm("process") == "process"      # ends in ss
    assert D.norm("status") == "status"        # ends in us
    assert D.norm("analysis") == "analysis"    # ends in is
    assert D.norm("has") == "has"              # not longer than 3 letters
    assert D.norm("gas") == "gas"


def test_content_tokens_and_mask_set():
    assert D.content_tokens("Regulation of DNA-templated transcription in cells") == ["dna", "templated", "transcription"]
    assert D.content_tokens("a b 1 processes") == []         # short tokens and (normalised) stopwords
    m = D.mask_set("alpha signaling", ["alpha pathways"])
    assert m == frozenset({"alpha", "signaling", "pathway"})


def test_deletion_count_matches_apply_mask():
    texts = ["Gene Symbol A1 Kinases bind kinase substrates.", "Gene Symbol B2 kinase", "Gene Symbol C3"]
    mask = frozenset({"kinase", "substrate"})
    total = sum(D.apply_mask(t, mask)[1] for t in texts)
    counts = sum((D.body_token_counts(t) for t in texts), start=D.Counter())
    assert total == sum(counts[x] for x in mask) == 4


# ---------------------------------------------------------------- selection rules

def test_jaccard_filter_and_off_target():
    pos = {"T:A": set("abcdef"), "T:B": set("abcdeg"), "T:C": set("mnopqr"), "T:D": set("stuvwx")}
    mask = {"T:A": frozenset({"aa"}), "T:B": frozenset({"bb"}), "T:C": frozenset({"cc"}), "T:D": frozenset({"dd"})}
    dels = {t: 10 for t in pos}
    anc = {t: set() for t in pos}
    assert D.jaccard(pos["T:A"], pos["T:B"]) > 0.5
    sel = D.select_terms(list(pos), pos, mask, dels, anc, n_total=4)
    terms = [t for t, _ in sel]
    assert len(terms) == 3 and not {"T:A", "T:B"} <= set(terms)
    for t, o in sel:
        assert o != t and D.jaccard(pos[t], pos[o]) <= 0.05 and not (mask[t] & mask[o])


def test_off_target_rules():
    pos = {"P": set("abcdefgh"), "T": set("abcd"), "X": set("wxyz"), "Y": set("qrst"), "Z": set("ijkl")}
    mask = {"P": frozenset({"pp"}), "T": frozenset({"tt"}), "X": frozenset({"tt", "xx"}),
            "Y": frozenset({"yy"}), "Z": frozenset({"zz"})}
    dels = {"P": 40, "T": 40, "X": 40, "Y": 60, "Z": 49}
    anc = {"P": set(), "T": {"P"}, "X": set(), "Y": set(), "Z": set()}
    # P is an ancestor, X shares a mask token, Y is outside max(5, 25% of 40)=10; Z (diff 9) qualifies
    assert D.find_off_target("T", list(pos), pos, mask, dels, anc) == "Z"
    dels["Z"] = 51
    assert D.find_off_target("T", list(pos), pos, mask, dels, anc) is None


def test_pick_genept_files(tmp_path):
    for n in ["NCBI_summary_of_genes.json", "NCBI_UniProt_summary_of_genes.json",
              "GenePT_gene_embedding_ada_text.pickle", "GenePT_gene_protein_embedding_model_3_text.pickle."]:
        (tmp_path / n).write_text("x")
    j, p = D.pick_genept_files(str(tmp_path))
    assert j.endswith("NCBI_summary_of_genes.json") and p.endswith("GenePT_gene_embedding_ada_text.pickle")
    (tmp_path / "ncbi_extra.json").write_text("x")
    with pytest.raises(D.PairingError):
        D.pick_genept_files(str(tmp_path))


# ---------------------------------------------------------------- end-to-end on a tiny fake dataset

LEAVES = {  # term: (name, parent relation, parent, genes)
    "GO:0000101": ("alpha signaling", "is_a", "GO:0000010", range(1, 7)),
    "GO:0000102": ("beta transport", "is_a", "GO:0000010", range(7, 13)),
    "GO:0000103": ("gamma folding", "is_a", "GO:0000010", range(13, 19)),
    "GO:0000104": ("delta repair", "part_of", "GO:0000020", range(19, 25)),
    "GO:0000105": ("epsilon export", "part_of", "GO:0000020", range(25, 31)),
    "GO:0000106": ("zeta import", "part_of", "GO:0000020", range(31, 37)),
}


def g(i):
    return f"G{i:02d}"


def _obo(path, obsolete=()):
    s = ["format-version: 1.2", "data-version: releases/test", "",
         "[Term]", "id: GO:0008150", "name: biological_process", "namespace: biological_process", "",
         "[Term]", "id: GO:0000010", "name: parent one", "namespace: biological_process", "is_a: GO:0008150", "",
         "[Term]", "id: GO:0000020", "name: parent two", "namespace: biological_process", "is_a: GO:0008150", "",
         "[Term]", "id: GO:0000200", "name: unrelated thing", "namespace: biological_process", "is_a: GO:0008150", ""]
    for t, (name, rel, par, _) in LEAVES.items():
        s += ["[Term]", f"id: {t}", f"name: {name}", "namespace: biological_process"]
        if t == "GO:0000101":
            s += ['synonym: "alpha pathways" EXACT []', 'synonym: "omega thing" RELATED []']
        if t in obsolete:
            s += ["is_obsolete: true"]
        else:
            s += [f"is_a: {par}" if rel == "is_a" else f"relationship: part_of {par}"]
        s.append("")
    s += ["[Typedef]", "id: part_of", "name: part of", ""]
    path.write_text("\n".join(s))


def _gaf_line(sym, go, code="IDA", qual="involved_in"):
    return "\t".join(["UniProtKB", "P" + sym, sym, qual, go, "PMID:1", code, "", "P", "", "", "protein",
                      "taxon:9606", "20240101", "X", "", ""]) + "\n"


def _gaf(path, extra=()):
    with gzip.open(path, "wt") as f:
        f.write("!gaf-version: 2.2\n")
        for t, (_, _, _, genes) in LEAVES.items():
            for i in genes:
                f.write(_gaf_line(g(i), t))
        for i in range(37, 51):
            f.write(_gaf_line(g(i), "GO:0000200"))
        f.write(_gaf_line("G40", "GO:0000102", qual="NOT|involved_in"))   # NOT row: dropped
        f.write(_gaf_line("G41", "GO:0000103", code="IEA"))              # non-experimental: dropped
        for sym, go in extra:
            f.write(_gaf_line(sym, go))


def make_fake(root):
    raw = root / "raw"
    (raw / "genept").mkdir(parents=True)
    genes = [g(i) for i in range(1, 61)]
    with open(raw / "hgnc.txt", "w") as f:
        f.write("hgnc_id\tsymbol\tname\tlocus_group\tprev_symbol\talias_symbol\n")
        for x in genes:
            prev = "OLD60" if x == "G60" else ""
            alias = "AMB" if x in ("G58", "G59") else ""
            f.write(f"HGNC:{x}\t{x}\tn\tprotein-coding gene\t{prev}\t{alias}\n")
        f.write("HGNC:R1\tR1\tn\tnon-coding RNA\t\t\n")
    words = {}
    for t, (name, _, _, gs) in LEAVES.items():
        for i in gs:
            words[g(i)] = f"Involved in {name}."
    summ = {}
    for x in genes:
        key = "OLD60" if x == "G60" else x
        body = words.get(x, "Unrelated text here.")
        summ[key] = f"Gene Symbol {key} {body}"
    summ["AMB"] = "Gene Symbol AMB ambiguous alias."
    summ["R1"] = "Gene Symbol R1 non-coding."
    json.dump(summ, open(raw / "genept" / "NCBI_summary_of_genes.json", "w"))
    json.dump(summ, open(raw / "genept" / "NCBI_UniProt_summary_of_genes.json", "w"))
    pickle.dump({k: [0.0] * 4 for k in list(summ) + ["EXTRA1"]},
                open(raw / "genept" / "GenePT_gene_embedding_ada_text.pickle", "wb"))
    with open(raw / "g2v.txt", "w") as f:
        for x in genes:
            if x != "G56":                     # G56 missing from Gene2vec -> not in universe
                f.write(f"{x}\t0.1 0.2 0.3\n")
    for rel in ("baseline", "current"):
        (raw / rel).mkdir()
    _obo(raw / "baseline" / "go-basic.obo")
    _obo(raw / "current" / "go-basic.obo", obsolete=("GO:0000106",))
    _gaf(raw / "baseline" / "goa.gaf.gz")
    _gaf(raw / "current" / "goa.gaf.gz",
         extra=[("G37", "GO:0000101"), ("G38", "GO:0000101"), ("G39", "GO:0000101")]
         + [(g(i), "GO:0000200") for i in range(51, 56)])
    return {"hgnc": str(raw / "hgnc.txt"), "genept_dir": str(raw / "genept"), "gene2vec": str(raw / "g2v.txt"),
            "base_gaf": str(raw / "baseline" / "goa.gaf.gz"), "base_obo": str(raw / "baseline" / "go-basic.obo"),
            "cur_gaf": str(raw / "current" / "goa.gaf.gz"), "cur_obo": str(raw / "current" / "go-basic.obo")}


@pytest.fixture
def fake(tmp_path):
    paths = make_fake(tmp_path)
    R = D.build_design(paths, n_main=4, n_pilot=2, min_pos=4, max_pos=10, max_neg=10, min_new=2)
    return tmp_path, R


def test_universe_and_labels(fake):
    _, R = fake
    U = set(R["universe"])
    assert "G60" in U and R["keys"]["G60"][0] == "OLD60"       # unique previous symbol is mapped
    assert not {"G56", "R1", "AMB"} & U                        # missing from Gene2vec / non-coding / ambiguous alias
    assert len(U) == 59
    assert R["stats"]["hgnc_ambiguous_aliases_dropped"] == 1    # AMB points to G58 and G59
    assert R["stats"]["summary_map"] == {"direct": 60, "via_alias": 1, "unresolved": 1}  # 59 G-keys + R1; OLD60; AMB
    assert R["stats"]["base_annotated_genes"] == 50            # G01-G50 (G41's IEA row adds nothing)
    assert R["stats"]["cur_annotated_genes"] == 49            # +G51-G55, -G31-G36 (their only term is obsolete)


def test_selection_and_no_overlap(fake):
    _, R = fake
    rows, design = R["rows"], R["design"]
    assert [r["set"] for r in rows].count("main") == 4 and [r["set"] for r in rows].count("pilot") == 2
    assert {r["term"] for r in rows} == set(LEAVES)
    by = {r["term"]: r for r in rows}
    assert by["GO:0000101"]["mask"] == ["alpha", "pathway", "signaling"]
    assert "G40" not in design["GO:0000102"]["pos"]                     # NOT row dropped
    for r in rows:
        d = design[r["term"]]
        P, N, NP, TN = set(d["pos"]), set(d["neg"]), set(d["new_pos"]), set(d["tneg"])
        assert not P & N and not P & NP and not P & TN and not N & TN and not NP & TN
        assert len(N) <= 10 and len(TN) <= 10
        assert not set(r["mask"]) & set(r["off_mask"])
        assert r["deletions"] == 12 and r["off_deletions"] == 12
    # negatives exclude genes annotated to a direct parent (part_of counts as a parent)
    sib = {g(i) for i in range(1, 19)}
    assert not set(design["GO:0000101"]["neg"]) & sib
    assert design["GO:0000101"]["new_pos"] == ["G37", "G38", "G39"] and design["GO:0000101"]["eligible"]
    assert not design["GO:0000106"]["eligible"]                        # obsolete in current
    assert design["GO:0000106"]["temporal_status"] == "obsolete_or_missing_in_current"


def test_deterministic_and_refuses_overwrite(fake, tmp_path):
    root, R = fake
    R2 = D.build_design(make_fake(tmp_path / "again"), n_main=4, n_pilot=2, min_pos=4, max_pos=10,
                        max_neg=10, min_new=2)
    assert R2["rows"] == R["rows"] and R2["design"] == R["design"]
    out = root / "frozen"
    sums = D.write_frozen(str(out), R)
    assert len(sums) == 5 and (out / "terms.json").exists() and (out / "SHA256SUMS").exists()
    assert (out / "stoplist.txt").read_text().split() == list(D.STOPWORDS)
    with pytest.raises(FileExistsError):
        D.write_frozen(str(out), R)
