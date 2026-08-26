#!/usr/bin/env python3
"""Validate DOI-backed titles in the exported scientific-reference catalogue."""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


DOI_PATTERN = re.compile(r"^10\.\d{4,9}/\S+$", re.IGNORECASE)

# A record whose metadata contradicts Crossref is a different outcome from a
# catalogue that could not be read or a Crossref that could not be reached:
# the first is a defect in the repository, the second leaves the bibliography
# simply unverified.  Callers need to tell them apart.
EXIT_METADATA_FAILURE = 1
EXIT_CANNOT_CHECK = 2

# Crossref returns titles as the publisher stored them: with typographic
# markup, with Greek letters where the catalogue spells the letter name, and --
# for the older volume-numbered journals -- prefixed by the article number
# ("LXXIX." for Rutherford 1911).  None of that is a metadata discrepancy, so
# it is normalized away before the titles are compared.
MARKUP_PATTERN = re.compile(r"<[^>]*>")
ARTICLE_NUMBER_PATTERN = re.compile(r"^[ivxlcdm]+\.\s*")
GREEK_LETTER_NAMES = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
    "ε": "epsilon", "ζ": "zeta", "η": "eta", "θ": "theta",
    "ι": "iota", "κ": "kappa", "λ": "lambda", "μ": "mu",
    "ν": "nu", "ξ": "xi", "ο": "omicron", "π": "pi",
    "ρ": "rho", "σ": "sigma", "τ": "tau", "υ": "upsilon",
    "φ": "phi", "χ": "chi", "ψ": "psi", "ω": "omega",
}


def parse_sources(path: Path) -> list[dict[str, str]]:
    sources: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if line.startswith("[source:") and line.endswith("]"):
            current = {"id": line[8:-1], "line": str(line_number)}
            sources.append(current)
        elif line.startswith("["):
            current = None
        elif current is not None and "=" in line:
            key, value = line.split("=", 1)
            current[key] = value
    return sources


def display_title(value: str) -> str:
    """Publisher markup and line breaks removed, wording left untouched."""
    return " ".join(html.unescape(MARKUP_PATTERN.sub(" ", value)).split())


def normalized_title(value: str) -> str:
    value = html.unescape(MARKUP_PATTERN.sub(" ", value))
    value = unicodedata.normalize("NFKC", value).casefold()
    value = "".join(GREEK_LETTER_NAMES.get(character, character) for character in value)
    value = ARTICLE_NUMBER_PATTERN.sub("", value.strip())
    return " ".join(re.findall(r"[\w]+", value, flags=re.UNICODE))


class LookupUnavailable(Exception):
    """Crossref could not be consulted; says nothing about the metadata."""


def crossref_title(doi: str, timeout: float) -> str:
    encoded_doi = urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(
        f"https://api.crossref.org/works/{encoded_doi}",
        headers={
            "Accept": "application/json",
            "User-Agent": "CREM-reference-validator/1.0 (https://github.com/pytlotad/Positronium)",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise ValueError("DOI is not registered in Crossref") from error
        raise LookupUnavailable(f"Crossref returned HTTP {error.code}") from error
    except json.JSONDecodeError as error:
        raise LookupUnavailable(f"Crossref returned malformed JSON: {error}") from error
    except OSError as error:
        raise LookupUnavailable(f"Crossref unreachable: {error}") from error
    titles = payload.get("message", {}).get("title", [])
    if not titles or not isinstance(titles[0], str):
        raise ValueError("Crossref record has no title")
    return titles[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("catalogue", nargs="?", default="ScientificalReferences.txt")
    parser.add_argument("--timeout", type=float, default=15.0)
    args = parser.parse_args()

    try:
        sources = parse_sources(Path(args.catalogue))
    except OSError as error:
        print(f"cannot read catalogue: {error}", file=sys.stderr)
        return EXIT_CANNOT_CHECK

    failures = 0
    checked = 0
    skipped = 0
    unreachable = 0
    for source in sources:
        source_id = source.get("id", "<unknown>")
        doi = source.get("doi", "").strip()
        if not doi:
            skipped += 1
            print(f"SKIP {source_id}: no DOI")
            continue
        title = source.get("title", "").strip()
        expected_url = f"https://doi.org/{doi}"
        problems: list[str] = []
        if not DOI_PATTERN.fullmatch(doi):
            problems.append(f"invalid DOI syntax: {doi}")
        if source.get("url", "").strip() != expected_url:
            problems.append(f"URL must be {expected_url}")
        remote_title = ""
        unavailable = False
        if not problems:
            try:
                remote_title = crossref_title(doi, args.timeout)
                checked += 1
                if normalized_title(title) != normalized_title(remote_title):
                    problems.append(
                        f"title mismatch: catalogue={title!r}, "
                        f"Crossref={display_title(remote_title)!r}"
                    )
            except LookupUnavailable as error:
                unavailable = True
                unreachable += 1
                print(f"ERROR {source_id}: {error}", file=sys.stderr)
            except (ValueError, KeyError) as error:
                problems.append(str(error))
        if problems:
            failures += 1
            for problem in problems:
                print(f"FAIL {source_id}: {problem}", file=sys.stderr)
        elif not unavailable:
            print(f"PASS {source_id}: {doi} — {display_title(remote_title)}")

    print(f"Reference metadata: {checked} DOI checked, {skipped} without DOI, "
          f"{unreachable} unreachable, {failures} failed")
    if failures:
        return EXIT_METADATA_FAILURE
    if unreachable:
        print("Crossref could not be reached for every record above; their "
              "metadata is unverified, not wrong.", file=sys.stderr)
        return EXIT_CANNOT_CHECK
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
