#!/usr/bin/env python3
"""Query the cBioPortal REST API to list samples for a given cancer type.

Outputs a delimited table (TSV by default) with columns:
    Cancer Type, Patient ID, Sample ID, Study ID
"""

import argparse
import csv
import logging
import re
import sys

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DEFAULT_BASE_URL = "https://www.cbioportal.org/api"
DEFAULT_PAGE_SIZE = 5000
DEFAULT_TIMEOUT = 30.0

HEADERS = {
    "Accept": "application/json",
    "User-Agent": "cbioportal-samples-query/1.0",
}

GENERIC_QUERY_WORDS = {
    "cancer",
    "cancers",
    "tumor",
    "tumors",
    "neoplasm",
    "neoplasms",
    "malignancy",
    "of",
    "the",
    "and",
    "in",
}

logger = logging.getLogger(__name__)


class CBioPortalError(Exception):
    """Error raised when an API request, payload check, or lookup fails."""


def create_session():
    """Create a requests Session with retry backoff for transient HTTP errors."""
    session = requests.Session()
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def fetch_json(session, base_url, endpoint, params=None, timeout=DEFAULT_TIMEOUT):
    """Send a GET request to a cBioPortal REST endpoint and return parsed JSON."""
    url = f"{base_url.rstrip('/')}/{endpoint.lstrip('/')}"
    try:
        response = session.get(url, params=params, headers=HEADERS, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise CBioPortalError(f"HTTP request failed for {url}: {exc}")

    try:
        return response.json()
    except ValueError:
        raise CBioPortalError(f"Invalid JSON response from {url}")


def _match_exact_cancer_types(query_lower, all_cancer_types):
    """Match a lowercase query against cancerTypeId first, then name or shortName."""
    matched = set()
    for record in all_cancer_types:
        ct_id = str(record.get("cancerTypeId", "")).strip()
        if ct_id and ct_id.lower() == query_lower:
            matched.add(ct_id)
    if matched:
        return matched

    for record in all_cancer_types:
        ct_id = str(record.get("cancerTypeId", "")).strip()
        name = str(record.get("name", "")).strip().lower()
        short_name = str(record.get("shortName", "")).strip().lower()
        if ct_id and (name == query_lower or short_name == query_lower):
            matched.add(ct_id)
    return matched


def resolve_cancer_types(cancer_type_query, all_cancer_types, include_subtypes=True):
    """Resolve a user cancer-type string to a dict of {cancerTypeId: name}.

    Matching order:
    1. Exact match on cancerTypeId (case-insensitive).
    2. Exact match on name or shortName (case-insensitive).
    3. If input ends with 'cancer' or 'cancers' (e.g., 'Breast Cancer'), strip
       that trailing word and check for an exact OncoTree match.

    By default, includes child subtypes in the OncoTree hierarchy so queries
    for organ or mid-level categories (e.g., 'breast', 'nsclc', 'mel') include
    their specific histological subtypes.
    """
    raw_query = cancer_type_query.strip()
    query = raw_query.lower()
    if not query:
        raise CBioPortalError("Cancer type input cannot be empty.")

    id_to_name = {}
    children_by_parent = {}
    for record in all_cancer_types:
        ct_id = str(record.get("cancerTypeId", "")).strip()
        if not ct_id:
            continue
        name = str(record.get("name") or ct_id).strip()
        parent = str(record.get("parent", "")).strip()
        id_to_name[ct_id] = name
        if parent:
            children_by_parent.setdefault(parent, []).append(ct_id)

    matched_ids = _match_exact_cancer_types(query, all_cancer_types)

    # OncoTree names organ categories 'Breast', 'Lung', 'Prostate', etc. without 'Cancer'
    if not matched_ids:
        for suffix in (" cancers", " cancer"):
            if query.endswith(suffix):
                stripped = query[: -len(suffix)].strip()
                if stripped:
                    matched_ids = _match_exact_cancer_types(stripped, all_cancer_types)
                    if matched_ids:
                        logger.info(
                            "Matched '%s' to OncoTree entry: %s",
                            raw_query,
                            ", ".join(
                                f"{cid} ({id_to_name[cid]})"
                                for cid in sorted(matched_ids)
                            ),
                        )
                break

    if not matched_ids:
        tokens = [
            tok
            for tok in re.split(r"[^a-z0-9]+", query)
            if tok and tok not in GENERIC_QUERY_WORDS
        ]
        suggestions = []
        for cid, cname in sorted(id_to_name.items()):
            searchable = f"{cid} {cname}".lower()
            if query in searchable or (
                tokens and all(tok in searchable for tok in tokens)
            ):
                suggestions.append(f"{cid} ({cname})")
        hint = f" Close matches: {', '.join(suggestions[:8])}" if suggestions else ""
        raise CBioPortalError(
            f"Cancer type '{raw_query}' was not found in cBioPortal.{hint}"
        )

    explicitly_mixed = "mixed" in matched_ids
    if include_subtypes:
        stack = list(matched_ids)
        while stack:
            current_id = stack.pop()
            for child_id in children_by_parent.get(current_id, []):
                if child_id == "mixed" and not explicitly_mixed:
                    continue
                if child_id not in matched_ids:
                    matched_ids.add(child_id)
                    stack.append(child_id)

    return {cid: id_to_name[cid] for cid in sorted(matched_ids)}


def get_studies_for_cancer_types(all_studies, target_cancer_types):
    """Filter studies whose cancerTypeId belongs to target_cancer_types.

    Filters on study['cancerTypeId'] rather than /studies?keyword=... because
    keyword searches match substrings inside unrelated study names and IDs.
    """
    matched_studies = [
        study
        for study in all_studies
        if study.get("cancerTypeId") in target_cancer_types and study.get("studyId")
    ]
    return sorted(matched_studies, key=lambda s: str(s["studyId"]))


def find_ancestors_with_studies(target_ids, all_cancer_types, all_studies):
    """Find parent cancer types that have studies when a subtype has none."""
    parent_by_id = {}
    name_by_id = {}
    for rec in all_cancer_types:
        cid = str(rec.get("cancerTypeId", "")).strip()
        if cid:
            parent_by_id[cid] = str(rec.get("parent", "")).strip()
            name_by_id[cid] = str(rec.get("name") or cid).strip()

    study_counts = {}
    for study in all_studies:
        ct_id = str(study.get("cancerTypeId", "")).strip()
        if ct_id:
            study_counts[ct_id] = study_counts.get(ct_id, 0) + 1

    ancestors = []
    seen = set(target_ids)
    for start_id in sorted(target_ids):
        curr = parent_by_id.get(start_id, "")
        while curr and curr not in seen:
            seen.add(curr)
            count = study_counts.get(curr, 0)
            if count > 0 and curr != "other":
                label = name_by_id.get(curr, curr)
                noun = "study" if count == 1 else "studies"
                ancestors.append(f"{curr} ({label}, {count} {noun})")
            curr = parent_by_id.get(curr, "")
    return ancestors


def get_samples_for_study(
    session, base_url, study, page_size=DEFAULT_PAGE_SIZE, timeout=DEFAULT_TIMEOUT
):
    """Fetch all samples for one study via pageNumber/pageSize pagination."""
    if page_size <= 0:
        raise ValueError("page_size must be a positive integer.")

    study_id = study["studyId"]
    endpoint = f"/studies/{study_id}/samples"
    samples = []
    page_number = 0

    while True:
        params = {
            "projection": "SUMMARY",
            "pageSize": page_size,
            "pageNumber": page_number,
            "sortBy": "sampleId",
            "direction": "ASC",
        }
        page_items = fetch_json(
            session, base_url, endpoint, params=params, timeout=timeout
        )
        if not isinstance(page_items, list):
            raise CBioPortalError(
                f"Expected a JSON list of samples for study '{study_id}', "
                f"got {type(page_items).__name__}."
            )

        samples.extend(page_items)
        if len(page_items) < page_size:
            break
        page_number += 1

    expected_count = study.get("allSampleCount")
    if isinstance(expected_count, int) and len(samples) != expected_count:
        raise CBioPortalError(
            f"Sample count mismatch for study '{study_id}': "
            f"expected {expected_count} (allSampleCount), got {len(samples)}."
        )

    return samples


def collect_samples_for_cancer_type(
    session,
    cancer_type_query,
    base_url=DEFAULT_BASE_URL,
    include_subtypes=True,
    page_size=DEFAULT_PAGE_SIZE,
    timeout=DEFAULT_TIMEOUT,
):
    """Fetch cancer types, filter matching studies, and collect sample records."""
    all_cancer_types = fetch_json(
        session,
        base_url,
        "/cancer-types",
        params={"projection": "SUMMARY"},
        timeout=timeout,
    )
    if not isinstance(all_cancer_types, list):
        raise CBioPortalError("Expected a JSON list from /cancer-types.")

    target_cancer_types = resolve_cancer_types(
        cancer_type_query, all_cancer_types, include_subtypes=include_subtypes
    )
    preview_ids = list(target_cancer_types.items())[:5]
    preview_str = ", ".join(f"{cid} ({cname})" for cid, cname in preview_ids)
    if len(target_cancer_types) > 5:
        preview_str += f", ... ({len(target_cancer_types)} total)"
    logger.info("Resolved '%s' -> %s", cancer_type_query, preview_str)

    all_studies = fetch_json(
        session,
        base_url,
        "/studies",
        params={"projection": "SUMMARY"},
        timeout=timeout,
    )
    if not isinstance(all_studies, list):
        raise CBioPortalError("Expected a JSON list from /studies.")

    matching_studies = get_studies_for_cancer_types(all_studies, target_cancer_types)
    if not matching_studies:
        ancestors = find_ancestors_with_studies(
            target_cancer_types.keys(), all_cancer_types, all_studies
        )
        extra = ""
        if not include_subtypes:
            extra = " Re-run without --exact-type to include child subtypes."
        elif ancestors:
            extra = f" Broader parent types with studies: {', '.join(ancestors)}."
        logger.warning(
            "No studies found directly under cancer type ID(s): %s.%s",
            ", ".join(target_cancer_types.keys()),
            extra,
        )
        return []

    logger.info("Found %d matching studies.", len(matching_studies))

    rows = []
    for study in matching_studies:
        study_id = str(study["studyId"])
        ct_id = str(study.get("cancerTypeId", ""))
        cancer_type_label = target_cancer_types.get(ct_id, ct_id)

        study_samples = get_samples_for_study(
            session, base_url, study, page_size=page_size, timeout=timeout
        )
        logger.info("Study %s: retrieved %d samples.", study_id, len(study_samples))

        for sample in study_samples:
            rows.append(
                {
                    "Cancer Type": cancer_type_label,
                    "Patient ID": str(sample.get("patientId", "")),
                    "Sample ID": str(sample.get("sampleId", "")),
                    "Study ID": str(sample.get("studyId", study_id)),
                }
            )

    return rows


def write_samples_table(rows, output_stream, delimiter="\t"):
    """Write sample rows with the required columns to output_stream."""
    fieldnames = ["Cancer Type", "Patient ID", "Sample ID", "Study ID"]
    writer = csv.DictWriter(
        output_stream,
        fieldnames=fieldnames,
        delimiter=delimiter,
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Query the cBioPortal REST API for samples of a given cancer type."
    )
    parser.add_argument(
        "cancer_type",
        help=(
            "Cancer type code, short name, or name "
            "(e.g., 'acc', 'brca', 'Breast Cancer', 'Non-Small Cell Lung Cancer')."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file path (default: write to stdout).",
    )
    parser.add_argument(
        "--format",
        choices=["tsv", "csv"],
        default="tsv",
        help="Output delimiter format: 'tsv' (default) or 'csv'.",
    )
    parser.add_argument(
        "--exact-type",
        action="store_true",
        help="Match only the exact cancer type without including OncoTree child subtypes.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=DEFAULT_PAGE_SIZE,
        help=f"Samples fetched per API request page (default: {DEFAULT_PAGE_SIZE}).",
    )
    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help=f"cBioPortal API base URL (default: {DEFAULT_BASE_URL}).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_TIMEOUT,
        help=f"HTTP request timeout in seconds (default: {DEFAULT_TIMEOUT}).",
    )
    return parser.parse_args()


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )
    args = parse_args()
    delimiter = "," if args.format == "csv" else "\t"

    try:
        with create_session() as session:
            rows = collect_samples_for_cancer_type(
                session=session,
                cancer_type_query=args.cancer_type,
                base_url=args.base_url,
                include_subtypes=not args.exact_type,
                page_size=args.page_size,
                timeout=args.timeout,
            )

        if args.output:
            with open(args.output, "w", encoding="utf-8", newline="") as out_fh:
                write_samples_table(rows, out_fh, delimiter=delimiter)
            logger.info("Wrote %d sample rows to %s.", len(rows), args.output)
        else:
            write_samples_table(rows, sys.stdout, delimiter=delimiter)
            logger.info("Wrote %d sample rows to stdout.", len(rows))
    except (CBioPortalError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
