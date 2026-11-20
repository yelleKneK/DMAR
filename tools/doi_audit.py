#!/usr/bin/env python3
"""Compare every DOI in the roxygen references with its registered metadata.

For each \\doi{} in the roxygen comments of R/*.R, the script takes the
reference paragraph that contains it, asks Crossref for the metadata
registered for that DOI (and, when Crossref does not hold the DOI, the
doi.org resolver for citation-style JSON), and compares the first author's
surname, the year, the title, the volume, and the first and last page with
the text of the reference. A reference that disagrees in any of these is
listed with what was found. Only metadata is requested; no article is
downloaded.

This is a maintainer tool. It is not shipped (tools/ is excluded from the
tarball) and it never fails a build: the exit status is 0 whatever it
finds, and the last line of output is

    DOI audit: <n> of <m> references disagree with the registered metadata

which tools/release_gate.R prints. A disagreement is a prompt to open the
article, not a verdict: registries hold errors too (a DOI that resolves to
a different item of the same issue, a title in another capitalization, an
online-first year).

Usage, from the package root:

    python3 tools/doi_audit.py                 # print the audit
    python3 tools/doi_audit.py --csv out.csv   # also write every row
    python3 tools/doi_audit.py --cache c.json  # keep responses between runs
    python3 tools/doi_audit.py --root /path/to/package

Crossref asks automated clients to identify themselves with a contact
address. The script sends one only when the environment variable
DOI_AUDIT_CONTACT holds it (for example an address of the maintainer's
choosing); nothing that identifies a person is written in this file, and
without the variable the requests carry a generic client name alone.
"""
import argparse
import csv
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import unicodedata

CLIENT = "doi-audit/1.0"


def user_agent():
    contact = os.environ.get("DOI_AUDIT_CONTACT", "").strip()
    return "%s (mailto:%s)" % (CLIENT, contact) if contact else CLIENT


def norm(s):
    s = unicodedata.normalize("NFKD", s or "")
    s = "".join(c for c in s if not unicodedata.combining(c))
    return s.lower()


def strip_rd(s):
    s = re.sub(r"\\(emph|code|eqn|strong|dQuote|sQuote|bold|url|href)\{", "{", s)
    s = s.replace("\\&", "&").replace("\\%", "%")
    s = re.sub(r"[{}]", "", s)
    return re.sub(r"\s+", " ", s).strip()


def collect(root):
    """Every (file, line, doi, reference paragraph) in the roxygen of R/."""
    refs = []
    for f in sorted(glob.glob(os.path.join(root, "R", "*.R"))):
        with open(f, encoding="utf-8") as fh:
            lines = fh.read().split("\n")
        for i, ln in enumerate(lines):
            if not ln.startswith("#'"):
                continue
            for m in re.finditer(r"\\doi\{([^}]+)\}", ln):
                # Walk back to the start of the paragraph that holds the DOI.
                j = i
                while (j > 0 and lines[j - 1].startswith("#'")
                       and lines[j - 1].strip() != "#'"
                       and "@references" not in lines[j - 1]
                       and "\\doi{" not in lines[j - 1]):
                    j -= 1
                para = " ".join(l[2:].strip() for l in lines[j:i + 1])
                para = re.sub(r"^@references\s*", "", para)
                refs.append((os.path.basename(f), i + 1, m.group(1).strip(),
                             strip_rd(para)))
    return refs


def get(url, accept=None):
    # -g: a DOI can hold brackets, which curl would otherwise read as a range.
    cmd = ["curl", "-s", "-g", "-L", "-m", "30", "-A", user_agent()]
    if accept:
        cmd += ["-H", "Accept: " + accept]
    return subprocess.run(cmd + [url], capture_output=True, text=True,
                          encoding="utf-8", errors="replace").stdout


def years_of(m, keys):
    out = set()
    for k in keys:
        for dp in ((m.get(k) or {}).get("date-parts") or []):
            if dp and dp[0]:
                out.add(str(dp[0]))
    return sorted(out)


def fetch(doi, cache):
    if doi in cache and cache[doi].get("source") not in (None, "FAILED"):
        return cache[doi]
    raw = ""
    try:
        raw = get("https://api.crossref.org/works/" + doi)
        m = json.loads(raw)["message"]
        out = {"source": "crossref",
               "authors": [a.get("family") or a.get("name") or ""
                           for a in m.get("author", [])],
               "title": (m.get("title") or [""])[0],
               "container": (m.get("container-title") or [""])[0],
               "volume": m.get("volume"), "issue": m.get("issue"),
               "page": m.get("page"),
               "years": years_of(m, ("issued", "published-print",
                                     "published-online", "published")),
               "type": m.get("type")}
    except Exception:
        try:
            raw = get("https://doi.org/" + doi,
                      "application/vnd.citationstyles.csl+json")
            m = json.loads(raw)
            out = {"source": "doi.org",
                   "authors": [a.get("family") or a.get("literal") or ""
                               for a in m.get("author", [])],
                   "title": m.get("title") or "",
                   "container": (m.get("container-title")
                                 or m.get("publisher") or ""),
                   "volume": m.get("volume"), "issue": m.get("issue"),
                   "page": m.get("page"),
                   "years": years_of(m, ("issued", "published")),
                   "type": m.get("type")}
        except Exception as e:
            # Not cached: the next run asks again.
            return {"source": "FAILED", "error": str(e)[:80], "raw": raw[:120]}
    cache[doi] = out
    time.sleep(0.15)
    return out


def compare(cite, md):
    """The list of fields in which the reference text and the metadata differ."""
    if md.get("source") in (None, "FAILED"):
        return ["unresolved"]
    flags = []
    c = norm(cite)
    first = norm((md.get("authors") or [""])[0])
    first_last = first.split()[-1] if first else ""
    if first and first not in c and first_last not in c:
        flags.append("first_author(%s)" % first)
    cy = re.search(r"\((\d{4})[a-z]?(?:/\d{4})?[^)]*\)", cite)
    if cy and md.get("years") and cy.group(1) not in md["years"]:
        # A year one off is usual for an article published online first.
        near = any(abs(int(cy.group(1)) - int(y)) <= 1 for y in md["years"])
        flags.append("%s(cited %s; registered %s)" % (
            "year_one_off" if near else "year", cy.group(1),
            "/".join(md["years"])))
    title = re.sub(r"<[^>]+>", "", md.get("title") or "")
    words = [w for w in re.findall(r"[a-z0-9]+", norm(title)) if len(w) > 3]
    if words:
        share = sum(1 for w in words if w in c) / len(words)
        if share < 0.6:
            flags.append("title(%.2f of its words found: %s)" % (share, title[:70]))
    v = md.get("volume")
    if v and not re.search(r"(?<!\d)%s(?!\d)" % re.escape(str(v)), cite):
        flags.append("volume(registered %s)" % v)
    p = md.get("page")
    if p:
        for q in re.split(r"[-\u2013]+", p)[:2]:
            q = q.strip()
            if q and not re.search(r"(?<![\dA-Za-z])%s(?!\d)" % re.escape(q), cite):
                flags.append("page(registered %s)" % p)
                break
    return flags


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=".", help="package root (default: .)")
    ap.add_argument("--cache", default=None,
                    help="JSON file of responses to reuse and extend "
                         "(default: a file in the session's temporary folder)")
    ap.add_argument("--csv", default=None, help="write every row to this file")
    ap.add_argument("--quiet", action="store_true",
                    help="print the summary line only")
    args = ap.parse_args()

    if not os.path.exists(os.path.join(args.root, "DESCRIPTION")):
        print("DOI audit: not run (no DESCRIPTION in %s)" % args.root)
        return 0
    cache_path = args.cache or os.path.join(tempfile.gettempdir(),
                                            "doi_audit_cache.json")
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            cache = json.load(fh)

    refs = collect(args.root)
    rows = []
    for (f, ln, doi, para) in refs:
        md = fetch(doi, cache)
        rows.append({
            "file": f, "line": ln, "doi": doi,
            "flags": "; ".join(compare(para, md)), "cited": para,
            "registered": "%s | %s | %s | %s(%s), %s | %s | %s" % (
                ", ".join(md.get("authors") or [])[:80],
                "/".join(md.get("years") or []),
                re.sub(r"<[^>]+>", "", md.get("title") or "")[:120],
                md.get("volume"), md.get("issue"), md.get("page"),
                md.get("container"), md.get("source"))})
    with open(cache_path, "w") as fh:
        json.dump(cache, fh, indent=0)
    if args.csv and rows:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    # A reference repeated on several help pages counts once per distinct
    # (DOI, text) pair.
    distinct = {}
    for r in rows:
        distinct.setdefault((r["doi"], r["cited"]), r)
    flagged = [r for r in distinct.values() if r["flags"]]
    unresolved = [r for r in flagged if r["flags"] == "unresolved"]
    if not args.quiet:
        print("roxygen \\doi occurrences: %d; distinct references: %d; "
              "distinct DOIs: %d" % (len(rows), len(distinct),
                                     len({r["doi"] for r in rows})))
        for r in sorted(flagged, key=lambda r: (r["file"], r["line"])):
            print("\n[%s:%s] %s\n  differs in: %s\n  cited:      %s\n"
                  "  registered: %s" % (r["file"], r["line"], r["doi"],
                                        r["flags"], r["cited"][:330],
                                        r["registered"]))
        print()
    if unresolved:
        print("DOI audit: %d references could not be looked up (no network, "
              "or the DOI is not registered)" % len(unresolved))
    print("DOI audit: %d of %d references disagree with the registered "
          "metadata" % (len(flagged) - len(unresolved), len(distinct)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
