#!/usr/bin/env python3
import argparse
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from anarci import anarci

# ANARCI allowed species list for germline assignment
ALLOWED_SPECIES = {"human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca"}

# IMGT region ranges
IMGT_RANGES = {
    "fr1": (1, 26),
    "cdr1": (27, 38),
    "fr2": (39, 55),
    "cdr2": (56, 65),
    "fr3": (66, 104),
    "cdr3": (105, 117),
    "fr4": (118, 129),
}


def chain_allow(chain_val: str) -> Set[str]:
    c = str(chain_val).strip().upper()
    if c == "H":
        return {"H"}
    if c == "L":
        return {"K", "L"}
    if c == "K":
        return {"K"}
    # For other labels (e.g., PDB chain letters), allow all Ig/TCR chains
    return {"A", "B", "D", "G", "H", "K", "L"}


def normalize_species(val: object) -> Optional[str]:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return None
    s = str(val).strip().lower()
    if not s:
        return None
    if s in ALLOWED_SPECIES:
        return s
    # Common scientific names -> ANARCI species keys
    if "homo sapiens" in s:
        return "human"
    if "mus musculus" in s:
        return "mouse"
    if "rattus" in s:
        return "rat"
    if "oryctolagus" in s:
        return "rabbit"
    if "macaca" in s:
        return "rhesus"
    if "sus scrofa" in s or s == "pig":
        return "pig"
    if "vicugna" in s or "alpaca" in s:
        return "alpaca"
    return None


def aligned_seq(numbering: Sequence[Tuple[Tuple[int, str], str]]) -> str:
    return "".join(aa for ((pos, ins), aa) in numbering)


def region_seq(numbering: Sequence[Tuple[Tuple[int, str], str]], start: int, end: int) -> str:
    return "".join(
        aa
        for ((pos, ins), aa) in numbering
        if start <= pos <= end and aa != "-"
    )


def numbering_string(numbering: Sequence[Tuple[Tuple[int, str], str]]) -> str:
    parts: List[str] = []
    for (pos, ins), aa in numbering:
        if aa == "-":
            continue
        suffix = ins.strip()
        label = f"{pos}{suffix}" if suffix else f"{pos}"
        parts.append(f"{label}:{aa}")
    return " ".join(parts)


def annotate(
    df: pd.DataFrame,
    sequence_col: str,
    chain_col: str,
    species_col: Optional[str],
    include_numbering: bool,
) -> pd.DataFrame:
    rows: List[Dict] = []

    group_cols: List[str] = [chain_col]
    if species_col and species_col in df.columns:
        group_cols.append(species_col)

    for key, grp in df.groupby(group_cols, dropna=False):
        if isinstance(key, tuple):
            chain_val, species_val = key
        else:
            chain_val, species_val = key, None
        allow = chain_allow(chain_val)
        allowed_species = normalize_species(species_val)
        seqs: List[Tuple[str, str]] = []
        idxs: List[int] = []
        for idx, seq in grp[sequence_col].items():
            if not isinstance(seq, str) or not seq.strip():
                seq = ""
            seqs.append((str(idx), seq))
            idxs.append(idx)

        if allowed_species:
            numberings, hit_tables, _ = anarci(
                seqs,
                assign_germline=True,
                allow=allow,
                allowed_species=[allowed_species],
            )
        else:
            numberings, hit_tables, _ = anarci(
                seqs, assign_germline=True, allow=allow
            )

        for idx, seq, nums, hits in zip(
            idxs, grp[sequence_col], numberings, hit_tables
        ):
            base = df.loc[idx].to_dict()
            base["chain_input"] = chain_val

            if not nums or not hits:
                base.update(
                    {
                        "anarci_status": "no_hit",
                        "chain_assigned": None,
                        "species": None,
                        "v_gene": None,
                        "v_gene_species": None,
                        "v_gene_score": None,
                        "j_gene": None,
                        "j_gene_species": None,
                        "j_gene_score": None,
                        "imgt_alignment": None,
                        "fr1": None,
                        "cdr1": None,
                        "fr2": None,
                        "cdr2": None,
                        "fr3": None,
                        "cdr3": None,
                        "fr4": None,
                    }
                )
                if include_numbering:
                    base["imgt_numbering"] = None
                rows.append(base)
                continue

            numbering, start, end = nums[0]
            hit = hits[0]

            v_gene = hit.get("germlines", {}).get("v_gene")
            j_gene = hit.get("germlines", {}).get("j_gene")

            v_gene_name = v_gene[0][1] if v_gene and v_gene[0] else None
            v_gene_species = v_gene[0][0] if v_gene and v_gene[0] else None
            v_gene_score = v_gene[1] if v_gene and len(v_gene) > 1 else None

            j_gene_name = j_gene[0][1] if j_gene and j_gene[0] else None
            j_gene_species = j_gene[0][0] if j_gene and j_gene[0] else None
            j_gene_score = j_gene[1] if j_gene and len(j_gene) > 1 else None

            base.update(
                {
                    "anarci_status": "ok",
                    "chain_assigned": hit.get("chain_type"),
                    "species": hit.get("species"),
                    "v_gene": v_gene_name,
                    "v_gene_species": v_gene_species,
                    "v_gene_score": v_gene_score,
                    "j_gene": j_gene_name,
                    "j_gene_species": j_gene_species,
                    "j_gene_score": j_gene_score,
                    "imgt_alignment": aligned_seq(numbering),
                    "fr1": region_seq(numbering, *IMGT_RANGES["fr1"]),
                    "cdr1": region_seq(numbering, *IMGT_RANGES["cdr1"]),
                    "fr2": region_seq(numbering, *IMGT_RANGES["fr2"]),
                    "cdr2": region_seq(numbering, *IMGT_RANGES["cdr2"]),
                    "fr3": region_seq(numbering, *IMGT_RANGES["fr3"]),
                    "cdr3": region_seq(numbering, *IMGT_RANGES["cdr3"]),
                    "fr4": region_seq(numbering, *IMGT_RANGES["fr4"]),
                }
            )
            if include_numbering:
                base["imgt_numbering"] = numbering_string(numbering)
            rows.append(base)

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run ANARCI on a CSV of antibody/TCR chains and extract genes/CDRs."
    )
    parser.add_argument(
        "--input",
        default="data/ab_chains.csv",
        help="Input CSV path (default: data/ab_chains.csv)",
    )
    parser.add_argument(
        "--output",
        default="output/anarci_annotations.csv",
        help="Output CSV path (default: output/anarci_annotations.csv)",
    )
    parser.add_argument(
        "--sequence-col",
        default="Sequence",
        help="Column name containing sequences (default: Sequence)",
    )
    parser.add_argument(
        "--chain-col",
        default="Chain",
        help="Column name containing chain label (default: Chain)",
    )
    parser.add_argument(
        "--species-col",
        default="species",
        help=(
            "Optional column with species for germline filtering "
            "(default: species)"
        ),
    )
    parser.add_argument(
        "--include-numbering",
        action="store_true",
        help="Add an imgt_numbering column like '1:E 2:V ...'",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if args.sequence_col not in df.columns:
        raise SystemExit(f"Missing sequence column: {args.sequence_col}")
    if args.chain_col not in df.columns:
        raise SystemExit(f"Missing chain column: {args.chain_col}")

    out_df = annotate(
        df,
        sequence_col=args.sequence_col,
        chain_col=args.chain_col,
        species_col=args.species_col,
        include_numbering=args.include_numbering,
    )
    out_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
