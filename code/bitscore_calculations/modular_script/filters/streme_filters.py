"""
filters/streme_filters.py
=========================
Parse a streme.xml file and provide three motif-allowlist filters:

    filter_by_evalue(path, threshold)      → set of motif names
    filter_by_pvalue(path, threshold)      → set of motif names
    filter_by_num_sites(path, min_sites)   → set of motif names

Each function returns a set of motif name strings (e.g. {"10-GCUGGAGCUGGGRU",
"8-GGCUUCCUGCYCAMC"}) — the same names that appear in column 3 of the
narrowPeak file. Pass the returned set to read_best_site() as `allowed_motifs`
to restrict scoring to motifs that pass the chosen threshold.

Why streme.xml rather than streme.html?
----------------------------------------
The XML file is produced alongside streme.html in every STREME output
directory and is explicitly designed for machine parsing. The HTML layout
can change between MEME Suite versions; the XML schema is stable. Python's
built-in xml.etree.ElementTree is used — no extra dependencies needed.

streme.xml motif element attributes used here
----------------------------------------------
  alt            → the consensus name (e.g. "10-GCUGGAGCUGGGRU")
                   This matches column 3 of the narrowPeak / FIMO output.
  test_evalue    → E-value from the hold-out test set
  test_pvalue    → p-value from the hold-out test set
  train_pos_count → number of positive training sequences with a site
                    (used as the "number of sites" filter)

Note: if STREME ran without a hold-out set (very small input), test_evalue
and test_pvalue may be absent. The functions warn and skip those motifs when
the requested attribute is missing.
"""

import sys
import xml.etree.ElementTree as ET


# ---------------------------------------------------------------------------
# INTERNAL HELPER — parse streme.xml once, return a list of motif dicts
# ---------------------------------------------------------------------------
# All three public filter functions call this; it is not exported.
#
# Each returned dict has the shape:
#   {
#     "name":       str,    # the consensus name used in narrowPeak col 3
#     "evalue":     float or None,
#     "pvalue":     float or None,
#     "num_sites":  int   or None,
#   }
#
# How the XML is navigated:
#   ET.parse()  → loads the whole file into memory as a tree
#   root.iter("motif")  → yields every <motif> element anywhere in the tree
#   element.get("attr") → reads one XML attribute (returns None if absent)

def _parse_streme_xml(path):
    """
    Parse streme.xml and return a list of per-motif stat dicts.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    motifs : list[dict]
    """
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        sys.exit(f"ERROR: Could not parse XML file '{path}': {exc}")

    root = tree.getroot()
    motifs = []

    for elem in root.iter("motif"):
        name = elem.get("id")
        if name is None:
            name = elem.get("alt", "UNKNOWN")

        def _float(attr):
            val = elem.get(attr)
            return float(val) if val is not None else None

        def _int(attr):
            val = elem.get(attr)
            return int(val) if val is not None else None

        motifs.append({
            "name":      name,
            "evalue":    _float("test_evalue"),
            "pvalue":    _float("test_pvalue"),
            "num_sites": _int("train_pos_count"),
        })

    if not motifs:
        print(
            f"  [WARNING] No <motif> elements found in '{path}'. "
            "Check that this is a streme.xml file.",
            file=sys.stderr,
        )

    return motifs


# ---------------------------------------------------------------------------
# PUBLIC FILTER 1 — E-value threshold
# ---------------------------------------------------------------------------
# Returns the set of motif names whose test E-value is <= threshold.
# Lower E-value = more significant.

def filter_by_evalue(path, threshold):
    """
    Return motif names whose E-value is at or below ``threshold``.

    Parameters
    ----------
    path : str or Path
        Path to streme.xml.
    threshold : float
        Maximum allowed E-value (e.g. 0.05).

    Returns
    -------
    allowed : set[str]
        Motif names passing the filter. Pass to read_best_site(allowed_motifs=).
    """
    motifs  = _parse_streme_xml(path)
    allowed = set()
    skipped = 0

    for m in motifs:
        if m["evalue"] is None:
            print(
                f"  [WARNING] Motif '{m['name']}' has no E-value — skipping "
                "(STREME may not have had enough sequences for a hold-out set).",
                file=sys.stderr,
            )
            skipped += 1
            continue
        if m["evalue"] <= threshold:
            allowed.add(m["name"])

    print(
        f"  E-value filter (<= {threshold}): "
        f"{len(allowed)} / {len(motifs) - skipped} motifs pass"
        + (f" ({skipped} had no E-value)" if skipped else ""),
    )
    return allowed


# ---------------------------------------------------------------------------
# PUBLIC FILTER 2 — P-value threshold
# ---------------------------------------------------------------------------

def filter_by_pvalue(path, threshold):
    """
    Return motif names whose P-value is at or below ``threshold``.

    Parameters
    ----------
    path : str or Path
        Path to streme.xml.
    threshold : float
        Maximum allowed P-value (e.g. 0.05).

    Returns
    -------
    allowed : set[str]
    """
    motifs  = _parse_streme_xml(path)
    allowed = set()
    skipped = 0

    for m in motifs:
        if m["pvalue"] is None:
            print(
                f"  [WARNING] Motif '{m['name']}' has no P-value — skipping.",
                file=sys.stderr,
            )
            skipped += 1
            continue
        if m["pvalue"] <= threshold:
            allowed.add(m["name"])

    print(
        f"  P-value filter (<= {threshold}): "
        f"{len(allowed)} / {len(motifs) - skipped} motifs pass"
        + (f" ({skipped} had no P-value)" if skipped else ""),
    )
    return allowed


# ---------------------------------------------------------------------------
# PUBLIC FILTER 3 — Minimum number of sites
# ---------------------------------------------------------------------------

def filter_by_num_sites(path, min_sites):
    """
    Return motif names found in at least ``min_sites`` positive sequences.

    Parameters
    ----------
    path : str or Path
        Path to streme.xml.
    min_sites : int
        Minimum number of positive training sequences that must contain a
        site for the motif to be retained.

    Returns
    -------
    allowed : set[str]
    """
    motifs  = _parse_streme_xml(path)
    allowed = set()
    skipped = 0

    for m in motifs:
        if m["num_sites"] is None:
            print(
                f"  [WARNING] Motif '{m['name']}' has no site count — skipping.",
                file=sys.stderr,
            )
            skipped += 1
            continue
        if m["num_sites"] >= min_sites:
            allowed.add(m["name"])

    print(
        f"  Site-count filter (>= {min_sites}): "
        f"{len(allowed)} / {len(motifs) - skipped} motifs pass"
        + (f" ({skipped} had no site count)" if skipped else ""),
    )
    return allowed
