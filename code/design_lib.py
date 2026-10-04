"""Stage 2 design library: gene universe, GO labels, summary tokens and masks, term selection,
off-target matches, negatives and temporal sets. Text and labels only; nothing here trains a
model or computes a score on real data."""
import gzip
import hashlib
import json
import os
import pickle
import random
import re
from collections import Counter, defaultdict

ROOT_BP = "GO:0008150"
BP = "biological_process"
EXPERIMENTAL = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}

# Stoplist exactly as given in the Stage 2 instructions.
STOPWORDS = tuple("""a an and are as at be been but by can for from has have in into is it its may of on or
that the their these this those to via was were which with within without other than such also both either
each not process processes regulation regulating positive negative cellular cell cells response responses
biological organismal multicellular organism involved establishment maintenance protein proteins gene genes
level levels activity""".split())

TOKEN = re.compile(r"[A-Za-z0-9]+")
TOKEN_WITH_SPACE = re.compile(r"([A-Za-z0-9]+)( ?)")
PREFIX = re.compile(r"^Gene Symbol \S+(?: |$)")


class PairingError(RuntimeError):
    """Raised when the GenePT text or embedding file cannot be picked unambiguously."""


# ---------------------------------------------------------------- text

def norm(tok):
    """Lowercase; strip a plural 's' if the word is longer than 3 letters and does not end in ss/us/is."""
    t = tok.lower()
    if len(t) > 3 and t.endswith("s") and not t.endswith(("ss", "us", "is")):
        t = t[:-1]
    return t


# Stopwords are matched after normalisation, so both the listed form and its normalised form are stopped
# (e.g. "processes" normalises to "processe").
STOP = frozenset(STOPWORDS) | frozenset(norm(w) for w in STOPWORDS)


def is_content(t):
    """t is an already-normalised token."""
    return len(t) >= 2 and t not in STOP


def split_prefix(text):
    """Split 'Gene Symbol X ' off the front. Returns (prefix, body); prefix is '' if absent."""
    m = PREFIX.match(text)
    if not m:
        return "", text
    return text[:m.end()], text[m.end():]


def content_tokens(text):
    return [t for t in (norm(x) for x in TOKEN.findall(text)) if is_content(t)]


def mask_set(name, exact_synonyms):
    toks = set(content_tokens(name))
    for s in exact_synonyms:
        toks.update(content_tokens(s))
    return frozenset(toks)


def apply_mask(text, mask):
    """Delete every body token whose normalised form is in `mask` (the prefix is never touched).
    A deleted token takes one following space with it. Returns (masked_text, n_deleted)."""
    prefix, body = split_prefix(text)
    n = 0

    def rep(m):
        nonlocal n
        if norm(m.group(1)) in mask:
            n += 1
            return ""
        return m.group(0)

    new_body = TOKEN_WITH_SPACE.sub(rep, body)
    return prefix + new_body, n


def body_token_counts(text):
    """Counter of normalised body tokens (prefix excluded)."""
    return Counter(norm(x) for x in TOKEN.findall(split_prefix(text)[1]))


# ---------------------------------------------------------------- genes

def read_hgnc(path):
    """Returns (approved symbols, protein-coding symbols, unique prev/alias -> approved map, n ambiguous aliases)."""
    approved, pc, amap = set(), set(), defaultdict(set)
    with open(path, encoding="utf-8") as f:
        idx = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        for line in f:
            r = line.rstrip("\n").split("\t")
            sym = r[idx["symbol"]]
            approved.add(sym)
            if r[idx["locus_group"]] == "protein-coding gene":
                pc.add(sym)
            for col in ("prev_symbol", "alias_symbol"):
                i = idx.get(col)
                if i is None or i >= len(r):
                    continue
                for a in r[i].strip('"').split("|"):
                    a = a.strip()
                    if a:
                        amap[a].add(sym)
    alias = {a: next(iter(s)) for a, s in amap.items() if len(s) == 1 and a not in approved}
    n_ambiguous = sum(1 for a, s in amap.items() if len(s) > 1 and a not in approved)
    return approved, pc, alias, n_ambiguous


def resolve(key, approved, alias):
    if key in approved:
        return key
    return alias.get(key)


def map_source(keys, approved, alias):
    """Map a source's keys to approved symbols. A key equal to the approved symbol wins; if only
    alias keys map to a symbol and there is more than one, the symbol is dropped as ambiguous.
    Returns ({symbol: source key}, stats Counter)."""
    by_sym, st = defaultdict(list), Counter()
    for k in keys:
        s = resolve(k, approved, alias)
        if s is None:
            st["unresolved"] += 1
        else:
            by_sym[s].append(k)
    out = {}
    for s, ks in by_sym.items():
        if s in ks:
            out[s] = s
            st["direct"] += 1
        elif len(ks) == 1:
            out[s] = ks[0]
            st["via_alias"] += 1
        else:
            st["ambiguous_multi_alias"] += 1
    return out, st


def pick_genept_files(genept_dir):
    """The one .json with 'ncbi' and not 'uniprot' in its name, and the one pickle with 'ada' and none of
    'uniprot', 'protein', 'model_3' in its name. Anything else raises PairingError."""
    names = sorted(os.listdir(genept_dir))
    low = {n: n.lower() for n in names}
    js = [n for n in names if low[n].endswith(".json") and "ncbi" in low[n] and "uniprot" not in low[n]]
    pk = [n for n in names if (".pickle" in low[n] or low[n].endswith(".pkl")) and "ada" in low[n]
          and not any(x in low[n] for x in ("uniprot", "protein", "model_3"))]
    if len(js) != 1 or len(pk) != 1:
        raise PairingError(f"ambiguous GenePT files: json candidates={js}, ada pickle candidates={pk}")
    return os.path.join(genept_dir, js[0]), os.path.join(genept_dir, pk[0])


def read_gene2vec_symbols(path):
    syms = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if parts:
                syms.append(parts[0])
    return syms


# ---------------------------------------------------------------- GO

def read_obo(path):
    """{id: {name, ns, obsolete, parents (is_a + part_of), exact (synonyms), alt}} for [Term] stanzas."""
    terms, cur, in_term = {}, None, False
    syn = re.compile(r'^"((?:[^"\\]|\\.)*)"\s+(\w+)')
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("["):
                in_term = line.strip() == "[Term]"
                cur = {"name": "", "ns": "", "obsolete": False, "parents": set(), "exact": [], "alt": []} if in_term else None
                continue
            if not in_term or not line or cur is None:
                continue
            key, _, val = line.partition(": ")
            if key == "id":
                terms[val.strip()] = cur
            elif key == "name":
                cur["name"] = val
            elif key == "namespace":
                cur["ns"] = val.strip()
            elif key == "is_obsolete":
                cur["obsolete"] = val.strip() == "true"
            elif key == "is_a":
                cur["parents"].add(val.split()[0])
            elif key == "relationship":
                p = val.split()
                if len(p) >= 2 and p[0] == "part_of":
                    cur["parents"].add(p[1])
            elif key == "alt_id":
                cur["alt"].append(val.strip())
            elif key == "synonym":
                m = syn.match(val)
                if m and m.group(2) == "EXACT":
                    cur["exact"].append(m.group(1).replace('\\"', '"'))
    return terms


def obo_data_version(path):
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("data-version:"):
                return line.split(":", 1)[1].strip()
            if line.startswith("["):
                break
    return None


def bp_terms(terms):
    return {t for t, d in terms.items() if d["ns"] == BP and not d["obsolete"]}


def bp_parents(terms, t):
    """Direct is_a/part_of parents of t that are live BP terms."""
    live = terms
    return {p for p in terms[t]["parents"] if p in live and live[p]["ns"] == BP and not live[p]["obsolete"]}


def ancestors(terms):
    """{BP term: set of all BP ancestors via is_a/part_of (excluding itself)}."""
    bp = bp_terms(terms)
    memo = {}

    def anc(t):
        if t in memo:
            return memo[t]
        memo[t] = set()  # guards against cycles
        out = set()
        for p in bp_parents(terms, t):
            out.add(p)
            out |= anc(p)
        memo[t] = out
        return out

    for t in sorted(bp):
        anc(t)
    return {t: memo[t] for t in bp}


def alt_map(terms):
    return {a: t for t, d in terms.items() if not d["obsolete"] for a in d["alt"]}


def read_gaf(path, terms, approved, alias, universe):
    """Experimental, non-NOT BP annotations of universe genes: {gene: set(direct GO ids)}, stats."""
    alts = alt_map(terms)
    direct, st = defaultdict(set), Counter()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!"):
                continue
            r = line.rstrip("\n").split("\t")
            if len(r) < 15 or r[8] != "P":
                continue
            if r[6] not in EXPERIMENTAL:
                continue
            st["exp_bp_rows"] += 1
            if "NOT" in r[3].split("|"):
                st["dropped_NOT"] += 1
                continue
            go = alts.get(r[4], r[4])
            if go != r[4]:
                st["alt_id_remapped"] += 1
            d = terms.get(go)
            if d is None or d["obsolete"] or d["ns"] != BP:
                st["dropped_obsolete_or_unknown_term"] += 1
                continue
            sym = resolve(r[2], approved, alias)
            if sym is None:
                st["dropped_unresolved_symbol"] += 1
                continue
            if sym not in universe:
                st["dropped_not_in_universe"] += 1
                continue
            direct[sym].add(go)
            st["kept"] += 1
    return dict(direct), st


def propagate(direct, anc):
    """{term: set(genes)} with true-path propagation over is_a/part_of."""
    pos = defaultdict(set)
    for g, ts in direct.items():
        full = set()
        for t in ts:
            full.add(t)
            full |= anc.get(t, set())
        for t in full:
            pos[t].add(g)
    return pos


# ---------------------------------------------------------------- design

def jaccard(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def find_off_target(t, candidates, pos, mask, dels, anc, off_jaccard=0.05):
    """Candidate term (not t, not an ancestor/descendant of t, no shared mask token, positive Jaccard
    <= off_jaccard) whose deletion count is closest to t's, within max(5, 25% of t's). Ties broken
    with random.Random('off-' + t). Returns None if there is no such term."""
    d = dels[t]
    tol = max(5, 0.25 * d)
    best, ties = None, []
    for o in sorted(candidates):
        if o == t or o in anc[t] or t in anc[o]:
            continue
        if mask[o] & mask[t]:
            continue
        if jaccard(pos[o], pos[t]) > off_jaccard:
            continue
        diff = abs(dels[o] - d)
        if diff > tol:
            continue
        if best is None or diff < best:
            best, ties = diff, [o]
        elif diff == best:
            ties.append(o)
    if not ties:
        return None
    return random.Random("off-" + t).choice(sorted(ties))


def select_terms(candidates, pos, mask, dels, anc, n_total, max_jaccard=0.5, off_jaccard=0.05, seed=0):
    """Shuffle candidates with random.Random(seed); accept a term if its positive Jaccard with every
    accepted term is <= max_jaccard and it has an off-target match. Returns [(term, off_term)]."""
    order = sorted(candidates)
    random.Random(seed).shuffle(order)
    accepted = []
    for t in order:
        if len(accepted) == n_total:
            break
        if any(jaccard(pos[t], pos[a]) > max_jaccard for a, _ in accepted):
            continue
        off = find_off_target(t, candidates, pos, mask, dels, anc, off_jaccard)
        if off is None:
            continue
        accepted.append((t, off))
    return accepted


def sample(pool, k, seed_str):
    pool = sorted(pool)
    if len(pool) <= k:
        return pool
    return sorted(random.Random(seed_str).sample(pool, k))


def build_design(paths, n_main=150, n_pilot=10, min_pos=20, max_pos=200, max_jaccard=0.5,
                 off_jaccard=0.05, max_neg=2000, min_new=10, seed=0):
    """Build the full design from files. `paths` keys: hgnc, genept_dir, gene2vec,
    base_gaf, base_obo, cur_gaf, cur_obo. Returns a dict with everything needed to freeze and report."""
    R = {"stats": {}}
    S = R["stats"]
    approved, pc, alias, n_amb = read_hgnc(paths["hgnc"])
    S["hgnc_approved"], S["hgnc_protein_coding"], S["hgnc_unique_aliases"], S["hgnc_ambiguous_aliases_dropped"] = \
        len(approved), len(pc), len(alias), n_amb

    json_path, ada_path = pick_genept_files(paths["genept_dir"])
    R["files"] = {"summaries": json_path, "ada": ada_path}
    with open(json_path, encoding="utf-8") as f:
        summaries = json.load(f)
    with open(ada_path, "rb") as f:
        ada_keys = list(pickle.load(f).keys())
    g2v_keys = read_gene2vec_symbols(paths["gene2vec"])

    maps = {}
    for name, keys in (("summary", summaries.keys()), ("ada", ada_keys), ("gene2vec", g2v_keys)):
        m, st = map_source(keys, approved, alias)
        maps[name] = {s: k for s, k in m.items() if s in pc}
        S[f"{name}_keys"] = len(set(keys))
        S[f"{name}_map"] = dict(st)
        S[f"{name}_protein_coding"] = len(maps[name])
    universe = sorted(set(maps["summary"]) & set(maps["ada"]) & set(maps["gene2vec"]))
    U = set(universe)
    R["universe"] = universe
    R["keys"] = {g: (maps["summary"][g], maps["ada"][g], maps["gene2vec"][g]) for g in universe}
    text = {g: summaries[maps["summary"][g]] for g in universe}
    R["text"] = text
    S["universe"] = len(universe)
    S["universe_empty_body"] = sum(1 for g in universe if not split_prefix(text[g])[1].strip())

    out = {}
    for rel in ("base", "cur"):
        terms = read_obo(paths[f"{rel}_obo"])
        anc = ancestors(terms)
        direct, st = read_gaf(paths[f"{rel}_gaf"], terms, approved, alias, U)
        pos = propagate(direct, anc)
        out[rel] = {"terms": terms, "anc": anc, "pos": pos, "annotated": set(direct)}
        S[f"{rel}_obo_version"] = obo_data_version(paths[f"{rel}_obo"])
        S[f"{rel}_gaf"] = dict(st)
        S[f"{rel}_annotated_genes"] = len(direct)
    B, C = out["base"], out["cur"]

    # candidates
    masks, cand = {}, []
    for t in sorted(bp_terms(B["terms"])):
        if t == ROOT_BP:
            continue
        n = len(B["pos"].get(t, ()))
        if not (min_pos <= n <= max_pos):
            continue
        d = B["terms"][t]
        m = mask_set(d["name"], d["exact"])
        if not m:
            continue
        masks[t] = m
        cand.append(t)
    S["candidates"] = len(cand)
    total = Counter()
    for g in sorted(B["annotated"]):
        total.update(body_token_counts(text[g]))
    dels = {t: sum(total[x] for x in masks[t]) for t in cand}
    pos_c = {t: B["pos"][t] for t in cand}
    sel = select_terms(cand, pos_c, masks, dels, B["anc"], n_main + n_pilot, max_jaccard, off_jaccard, seed)
    S["selected"] = len(sel)

    rows, design = [], {}
    for i, (t, o) in enumerate(sel):
        pos = B["pos"][t]
        excl = set(pos)
        for p in bp_parents(B["terms"], t):
            excl |= B["pos"].get(p, set())
        neg = sample(B["annotated"] - excl, max_neg, "neg-" + t)
        ct = C["terms"].get(t)
        if ct is None or ct["obsolete"] or ct["ns"] != BP:
            new_pos, tneg, eligible, why = [], [], False, "obsolete_or_missing_in_current"
        else:
            pc_t = C["pos"].get(t, set())
            new_pos = sorted(pc_t - pos)
            texcl = set(pc_t) | pos | set(neg)
            for p in bp_parents(C["terms"], t):
                texcl |= C["pos"].get(p, set())
            tneg = sample(C["annotated"] - texcl, max_neg, "tneg-" + t)
            eligible = len(new_pos) >= min_new
            why = "ok" if eligible else "too_few_new_pos"
        rows.append({
            "set": "main" if i < n_main else "pilot", "term": t, "name": B["terms"][t]["name"],
            "n_pos": len(pos), "mask": sorted(masks[t]), "deletions": dels[t],
            "off_term": o, "off_name": B["terms"][o]["name"], "off_mask": sorted(masks[o]),
            "off_deletions": dels[o]})
        design[t] = {"pos": sorted(pos), "neg": neg, "new_pos": new_pos, "tneg": tneg, "eligible": eligible,
                     "temporal_status": why}
    R["rows"], R["design"], R["masks"] = rows, design, masks
    return R


def write_frozen(out_dir, R):
    """Write frozen/ files. Refuses (FileExistsError) if terms.json already exists."""
    if os.path.exists(os.path.join(out_dir, "terms.json")):
        raise FileExistsError(f"{out_dir}/terms.json already exists; refusing to overwrite")
    os.makedirs(out_dir, exist_ok=True)
    files = {}
    files["stoplist.txt"] = "\n".join(STOPWORDS) + "\n"
    files["universe.txt"] = "\n".join(R["universe"]) + "\n"
    files["universe_keys.tsv"] = "symbol\tsummary_key\tada_key\tgene2vec_key\n" + "".join(
        f"{g}\t{a}\t{b}\t{c}\n" for g, (a, b, c) in sorted(R["keys"].items()))
    files["design.json"] = json.dumps(R["design"], indent=1, sort_keys=True) + "\n"
    files["terms.json"] = json.dumps(R["rows"], indent=1) + "\n"  # written last before SHA256SUMS
    sums = []
    for name in ("stoplist.txt", "universe.txt", "universe_keys.tsv", "design.json", "terms.json"):
        with open(os.path.join(out_dir, name), "w", encoding="utf-8") as f:
            f.write(files[name])
        sums.append(f"{hashlib.sha256(files[name].encode('utf-8')).hexdigest()}  {name}\n")
    with open(os.path.join(out_dir, "SHA256SUMS"), "w", encoding="utf-8") as f:
        f.writelines(sums)
    return sums
