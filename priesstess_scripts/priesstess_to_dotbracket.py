"""
priesstess_to_dotbracket.py

Converts PRIESSTESS/STREME sequence-structure motifs into a
(sequence, dot-bracket) pair that can be pasted directly into forna
(http://rna.tbi.univie.ac.at/forna) or fed to VARNA / RNAplot.

Background reading:
  Laverty et al. 2022, "PRIESSTESS: interpretable, high-performing models
  of the sequence and structure preferences of RNA-binding proteins",
  Nucleic Acids Research 50(19):e111. DOI: 10.1093/nar/gkac694
  (defines the 2-, 4- and 7-letter structural alphabets used below)

  Kerpedjiev, Hammer & Hofacker 2015, "forna (force-directed RNA):
  simple and effective online RNA secondary structure diagrams",
  Bioinformatics 31(20):3377-9. DOI: 10.1093/bioinformatics/btv372
"""

import os
import re

# ---------------------------------------------------------------------------
# Real PRIESSTESS_output layout (confirmed directly from the PRIESSTESS
# bash script, not assumed):
#
#   PRIESSTESS_output/
#     PRIESSTESS_model.sav
#     PRIESSTESS_model_weights.tab
#     seq-4/            streme.txt, PFM-1.txt, PFM-2.txt, ...
#     seq-struct-8/      streme.txt, PFM-1.txt, ...
#     seq-struct-16/     streme.txt, PFM-1.txt, ...
#     seq-struct-28/     streme.txt, PFM-1.txt, ...
#     struct-2/          streme.txt, PFM-1.txt, ...
#     struct-4/          streme.txt, PFM-1.txt, ...
#     struct-7/          streme.txt, PFM-1.txt, ...
#
# One subdirectory per alphabet (named exactly as below), each holding its
# own streme.txt (which defines that alphabet) and PFM-N.txt files. The
# model files sit at the top level, not inside any alphabet folder. A run
# only has the subdirectories for whichever alphabets were requested with
# -alph, so not all seven will necessarily be present.
# ---------------------------------------------------------------------------
KNOWN_ALPHABET_DIRS = [
    "seq-4",
    "seq-struct-8",
    "seq-struct-16",
    "seq-struct-28",
    "struct-2",
    "struct-4",
    "struct-7",
]

# ---------------------------------------------------------------------------
# Structural code tables (Laverty et al. 2022, Methods: "Folding and
# annotating structures"). Note that the SAME letter means different
# things in different tables -- e.g. "L" is a hairpin loop in the 4-letter
# alphabet but a 5' paired strand in the 7-letter alphabet. The alphabet
# must always be identified before a letter is interpreted.
# ---------------------------------------------------------------------------
NT_CHARS = set("ACGU")

STRUCT2 = {"U": "unpaired", "P": "paired_ambiguous"}
STRUCT4 = {
    "P": "paired_ambiguous",
    "L": "loop_hairpin",
    "U": "unpaired",
    "M": "loop_multi_internal_bulge",
}
STRUCT7 = {
    "E": "unpaired",
    "B": "loop_bulge",
    "H": "loop_hairpin",
    "L": "paired_5prime",
    "R": "paired_3prime",
    "T": "loop_internal",
    "M": "loop_multi",
}


def parse_alphabet_block(streme_path):
    """
    Read the ALPHABET block at the top of a STREME output file and return
    (alphabet_name, {symbol: label_or_None}).

    The label format varies by alphabet type, and isn't always quoted:
      combined:        A "A-P" 660000        -> label = "A-P"
      sequence-only:    A 660000               -> no label at all (None)
      structure-only:   L "Hairpin loop" 660000 -> label = "Hairpin loop"
                                                    (may contain spaces)
    The hex colour at the end of each line is discarded -- it's not needed
    to decode anything.
    """
    text = open(streme_path).read()
    match = re.search(r'ALPHABET\s+"([^"]+)".*?\n(.*?)\n\?\s*=', text, re.S)
    if not match:
        raise ValueError(f"Could not find an ALPHABET block in {streme_path}")
    name, body = match.group(1), match.group(2)

    line_pattern = re.compile(r'^(\S+)(?:\s+"([^"]*)")?\s+[0-9A-Fa-f]{6}\s*$')
    symbol_map = {}
    for line in body.strip().splitlines():
        line_match = line_pattern.match(line.strip())
        if not line_match:
            continue  # skips the "****" divider line etc.
        symbol, label = line_match.group(1), line_match.group(2)
        symbol_map[symbol] = label
    return name, symbol_map


def _pick_struct_table(struct_codes_seen):
    """Match a set of structure codes to the 2-, 4- or 7-letter semantic table it belongs to."""
    if struct_codes_seen == {"U", "P"}:
        return STRUCT2
    elif struct_codes_seen <= set(STRUCT4):
        return STRUCT4
    elif struct_codes_seen <= set(STRUCT7):
        return STRUCT7
    raise ValueError(f"Unrecognised structural code set: {struct_codes_seen}")


_COMBINED_LABEL = re.compile(r"^[ACGU]-[A-Z]$")


def classify_and_decode(symbol_map):
    """
    Turn {symbol: label_or_None} into {symbol: [nucleotide_or_None, semantic_context_or_None]}.

    Handles all three PRIESSTESS alphabet shapes, distinguished by looking
    at the alphabet as a whole rather than any single symbol in isolation
    (a single symbol's label can be genuinely ambiguous -- see the "U"
    nucleotide-vs-unpaired-code collision noted below):

      - sequence-only   : the full symbol set is exactly {A, C, G, U}.
                           There's no structural label at all -- the bare
                           symbol *is* the nucleotide.
      - combined        : every label is a two-part "nucleotide-code"
                           string, e.g. "A-P".
      - structure-only  : anything else. Here the symbol itself (not its
                           English label like "Hairpin loop") is the
                           structural code -- struct-2/4/7's symbols
                           (P, U, L, M, B, E, H, R, T) are exactly the
                           keys of STRUCT2/STRUCT4/STRUCT7, so the label
                           text doesn't need to be parsed at all.
    """
    symbols = set(symbol_map.keys())

    if symbols == NT_CHARS:
        return {symbol: [symbol, None] for symbol in symbol_map}

    labels = list(symbol_map.values())
    if all(label and _COMBINED_LABEL.match(label) for label in labels):
        decoded, struct_codes_seen = {}, set()
        for symbol, label in symbol_map.items():
            nt, struct_code = label.split("-")
            decoded[symbol] = [nt, struct_code]
            struct_codes_seen.add(struct_code)
        table = _pick_struct_table(struct_codes_seen)
    else:
        table = _pick_struct_table(symbols)
        decoded = {symbol: [None, symbol] for symbol in symbol_map}

    for symbol, (nt, struct_code) in decoded.items():
        decoded[symbol][1] = table[struct_code] if (struct_code and table) else None

    return decoded


def parse_pfm_file(pfm_path):
    """
    Parse a PRIESSTESS PFM-N.txt file: one row per alphabet symbol,
    tab-separated probabilities per motif position. Returns
    ({symbol: [prob_pos1, prob_pos2, ...]}, motif_width).
    """
    letter_probs, width = {}, None
    for line in open(pfm_path):
        parts = line.strip().split("\t")
        letter_probs[parts[0]] = [float(x) for x in parts[1:]]
        width = len(letter_probs[parts[0]])
    return letter_probs, width


def consensus_symbols(letter_probs, width):
    """Most probable alphabet symbol at each motif position."""
    return [max(letter_probs, key=lambda s: letter_probs[s][pos]) for pos in range(width)]


def build_seq_and_structure(symbols, decoded):
    """
    Turn a list of consensus symbols into (sequence, dot_bracket, context_notes, ambiguous_flag).

    Dot-bracket rules:
      paired_5prime      -> "("
      paired_3prime      -> ")"
      paired_ambiguous   -> "x"   (2- or 4-letter alphabets don't record which
                                    strand a paired base is on -- this can't be
                                    resolved from a single motif; see note below)
      anything else       -> "."   (all unpaired/loop flavours collapse to "."
                                    since dot-bracket itself can't distinguish
                                    hairpin loop vs bulge vs multiloop vs external)
      no structural info  -> "."   (sequence-only alphabet)
    """
    sequence, context_notes, dot_bracket = [], [], []
    ambiguous = False
    for symbol in symbols:
        nt, context = decoded[symbol]
        sequence.append(nt or "N")
        context_notes.append(context)
        if context == "paired_5prime":
            dot_bracket.append("(")
        elif context == "paired_3prime":
            dot_bracket.append(")")
        elif context == "paired_ambiguous":
            dot_bracket.append("x")
            ambiguous = True
        else:
            dot_bracket.append(".")
    return "".join(sequence), "".join(dot_bracket), context_notes, ambiguous


def motif_to_dotbracket(pfm_path, streme_path):
    """
    Convenience wrapper: given one PFM-N.txt file and the streme.txt that
    defines its alphabet, return a dict with sequence, dot-bracket, per-
    position context, and any warnings.
    """
    _, symbol_map = parse_alphabet_block(streme_path)
    decoded = classify_and_decode(symbol_map)
    letter_probs, width = parse_pfm_file(pfm_path)
    symbols = consensus_symbols(letter_probs, width)
    sequence, dot_bracket, context_notes, ambiguous = build_seq_and_structure(symbols, decoded)

    warnings = []
    if all(c is None for c in context_notes):
        warnings.append(
            "Sequence-only alphabet: no structural information exists for this "
            "motif. The dot-bracket string is all dots by default, not a "
            "genuine 'unstructured' prediction."
        )
    if ambiguous:
        warnings.append(
            "This alphabet marks paired bases without recording which strand "
            "(5' vs 3') they're on ('x' in the dot-bracket string). To resolve "
            "this, look for a second retained motif whose consensus is the "
            "reverse complement -- if one exists, they are likely the two "
            "sides of the same stem (see the L1L2 rung-matching example in "
            "chat). Otherwise this can't be auto-resolved from one motif alone."
        )

    return {
        "sequence": sequence,
        "dot_bracket": dot_bracket,
        "context_per_position": context_notes,
        "warnings": warnings,
    }


def parse_model_weights(weights_path):
    """
    Parse a PRIESSTESS_model_weights.tab file. Returns a list of
    (alphabet_PFM_name, weight) tuples for only the NONZERO entries
    (i.e. the motifs LASSO actually retained), sorted by |weight| descending.
    """
    retained = []
    for line in open(weights_path):
        name, weight = line.strip().split("\t")
        weight = float(weight)
        if weight != 0:
            retained.append((name, weight))
    retained.sort(key=lambda pair: abs(pair[1]), reverse=True)
    return retained


def discover_alphabet_dirs(priesstess_output_dir):
    """
    Walk a real PRIESSTESS_output directory and return
    {alphabet_name: (streme_path, pfm_dir)} for every alphabet subfolder
    that's actually present (a run doesn't necessarily use all seven).
    """
    pfm_lookup = {}
    for alphabet_name in KNOWN_ALPHABET_DIRS:
        alphabet_dir = os.path.join(priesstess_output_dir, alphabet_name)
        streme_path = os.path.join(alphabet_dir, "streme.txt")
        if os.path.isdir(alphabet_dir) and os.path.isfile(streme_path):
            pfm_lookup[alphabet_name] = (streme_path, alphabet_dir)
    return pfm_lookup


def summarize_retained_motifs(weights_path, pfm_lookup=None):
    """
    Report which motifs the final logistic regression model kept.

    pfm_lookup is optional: a dict mapping the alphabet name portion of each
    weights.tab entry (e.g. "seq-struct-16") to a (streme_path, pfm_dir)
    tuple, so that any retained motif you have files for gets fully decoded.
    Entries with no matching lookup are still reported by name/weight, just
    without a sequence/structure decode. In practice, build this dict with
    discover_alphabet_dirs() rather than by hand -- see
    process_priesstess_output() below.
    """
    retained = parse_model_weights(weights_path)
    if not retained:
        print("No motifs were retained with nonzero weight in this model.")
        return []

    results = []
    for name, weight in retained:
        alphabet_name, pfm_id = name.rsplit("_", 1)  # e.g. "seq-struct-16", "PFM-2"
        print(f"{name}  (weight={weight:.4f})")

        if pfm_lookup and alphabet_name in pfm_lookup:
            streme_path, pfm_dir = pfm_lookup[alphabet_name]
            pfm_path = f"{pfm_dir}/{pfm_id}.txt"
            try:
                decoded = motif_to_dotbracket(pfm_path, streme_path)
                print(f"   sequence:    {decoded['sequence']}")
                print(f"   dot-bracket: {decoded['dot_bracket']}")
                for w in decoded["warnings"]:
                    print(f"   note: {w}")
                results.append({"name": name, "weight": weight, **decoded})
            except FileNotFoundError:
                print(f"   (PFM file not found at {pfm_path} -- skipped)")
        else:
            print("   (no streme.txt/PFM directory supplied for this alphabet -- skipped)")
        print()

    return results


def write_dotbracket_file(results, output_path):
    """
    Write only fully-resolved motifs to one file in Vienna dot-bracket
    format:

        >name weight=...
        SEQUENCE
        STRUCTURE

    This is the format forna, VARNA and RNAplot all expect. "Resolved"
    means the dot_bracket string contains no 'x' -- i.e. it's either
    genuinely unpaired throughout, or it came from an alphabet (7-letter
    struct-7 / seq-struct-28) that records which strand a paired base is
    on. Motifs with any unresolved ambiguous pairing are deliberately left
    out of this file -- see write_unresolved_file -- because writing them
    here with 'x' silently turned into '.' would tell visualization
    software "this is unpaired" when the true answer is "paired, we just
    don't know which side," which is a materially different and
    misleading claim, not a safe simplification.
    """
    resolved = [r for r in results if "x" not in r["dot_bracket"]]
    with open(output_path, "w") as f:
        for result in resolved:
            header = f">{result['name']} weight={result['weight']:.4f}"
            if result["warnings"]:
                header += "  NOTE: " + " | ".join(result["warnings"])
            f.write(header + "\n")
            f.write(result["sequence"] + "\n")
            f.write(result["dot_bracket"] + "\n")
    return output_path, resolved


def write_unresolved_file(results, output_path):
    """
    Write motifs with unresolved ambiguous pairing to a separate,
    deliberately NON-dot-bracket file, so they can never be mistaken for
    a finished structure and accidentally fed into visualization software
    as-is. Each entry shows the sequence, the raw notation (with 'x' left
    in place to mark exactly which positions are paired-but-unknown-side),
    and a reminder of what would need to be done to resolve it.
    """
    unresolved = [r for r in results if "x" in r["dot_bracket"]]
    if not unresolved:
        return output_path, unresolved
    with open(output_path, "w") as f:
        f.write(
            "# Motifs below have at least one position PRIESSTESS annotated as\n"
            "# paired, but the alphabet used doesn't record which strand (5' vs\n"
            "# 3') that base is on -- 'x' marks those positions. These are NOT\n"
            "# valid dot-bracket and must not be fed to forna/VARNA/RNAplot as-is.\n"
            "# To resolve: check whether a same-run motif from a 7-letter alphabet\n"
            "# (struct-7 / seq-struct-28) covers the same region, or look for\n"
            "# another retained motif whose sequence is this one's reverse\n"
            "# complement (a likely stem partner).\n\n"
        )
        for result in unresolved:
            f.write(f">{result['name']} weight={result['weight']:.4f}\n")
            f.write(result["sequence"] + "\n")
            f.write(result["dot_bracket"] + "  (x = paired, side unknown)\n\n")
    return output_path, unresolved


def process_priesstess_output(priesstess_output_dir):
    """
    The main entry point for a real run: point this at a PRIESSTESS_output
    directory and it finds the weights file, discovers whichever alphabet
    subfolders are present, decodes every retained motif it can, prints a
    summary to the terminal as before, and saves the results into TWO
    files in that same output directory:

      PRIESSTESS_dotbracket_structures.dbn    -- fully resolved, ready
                                                   to visualize as-is
      PRIESSTESS_dotbracket_UNRESOLVED.txt    -- ambiguous-paired motifs
                                                   that need manual review
                                                   before visualizing
    """
    weights_path = os.path.join(priesstess_output_dir, "PRIESSTESS_model_weights.tab")
    pfm_lookup = discover_alphabet_dirs(priesstess_output_dir)
    print(f"Found alphabet folders: {', '.join(pfm_lookup) or '(none)'}")
    print()

    results = summarize_retained_motifs(weights_path, pfm_lookup=pfm_lookup)

    if not results:
        print("Nothing decoded, so no output files were written.")
        return results

    dotbracket_path = os.path.join(priesstess_output_dir, "PRIESSTESS_dotbracket_structures.dbn")
    unresolved_path = os.path.join(priesstess_output_dir, "PRIESSTESS_dotbracket_UNRESOLVED.txt")
    _, resolved = write_dotbracket_file(results, dotbracket_path)
    _, unresolved = write_unresolved_file(results, unresolved_path)

    print(f"Saved {len(resolved)} resolved motif(s) to: {dotbracket_path}")
    if unresolved:
        print(f"Saved {len(unresolved)} unresolved (ambiguous-paired) motif(s) to: {unresolved_path}")
        print("These need manual review before visualizing -- see the file's header comment.")

    return results


if __name__ == "__main__":
    import sys

    # Usage: python3 priesstess_to_dotbracket.py /path/to/PRIESSTESS_output
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "PRIESSTESS_output"
    process_priesstess_output(output_dir)
