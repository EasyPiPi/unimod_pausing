# ---------------------------------------------------------------------------
# Batch sliding-window GC content and GC skew for all .fa files in a directory
# ---------------------------------------------------------------------------

import argparse
from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VALID_BASES = frozenset("ACGT")
LOW_QUALITY_WINDOW = slice(499, 1499)   # region checked for valid-base ratio
LOW_QUALITY_THRESHOLD = 0.99


# ---------------------------------------------------------------------------
# Sequence metrics
# ---------------------------------------------------------------------------

def sliding_gc_metrics(seq: str, window_size: int) -> Tuple[List[float], List[float]]:
    """Return (gc_content, gc_skew) lists computed in sliding windows."""
    seq = seq.upper()
    n = len(seq)

    gc_content: List[float] = []
    gc_skew: List[float] = []

    for start in range(0, n - window_size + 1):
        sub = seq[start : start + window_size]
        g = sub.count("G")
        c = sub.count("C")

        gc_content.append((g + c) / window_size)
        gc_skew.append((g - c) / (g + c) if (g + c) > 0 else 0.0)

    return gc_content, gc_skew


def is_low_quality(sequence: str) -> bool:
    """Return True if the quality-check window has too many non-ACGT bases."""
    region = sequence[LOW_QUALITY_WINDOW]
    if not region:
        return True
    valid = sum(base in VALID_BASES for base in region)
    return (valid / len(region)) < LOW_QUALITY_THRESHOLD


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------

def process_fasta(fasta_path: Path, window_size: int, output_dir: Path) -> None:
    """Filter low-quality genes, then compute and write sliding-window metrics."""
    stem   = fasta_path.stem
    prefix = f"{stem}_{window_size}"

    gc_content_path  = output_dir / f"{prefix}_gc_content.csv"
    gc_skew_path     = output_dir / f"{prefix}_gc_skew.csv"
    low_quality_path = output_dir / f"{prefix}_low_quality_genes.csv"

    with open(gc_content_path, "w") as f_gc, \
         open(gc_skew_path,    "w") as f_skew, \
         open(low_quality_path, "w") as f_low:

        for record in SeqIO.parse(fasta_path, "fasta"):
            sequence  = str(record.seq).upper()
            gene_name = record.description.split("::", 1)[0]

            if is_low_quality(sequence):
                f_low.write(gene_name + "\n")
                continue

            gc_content, gc_skew = sliding_gc_metrics(sequence, window_size)

            f_gc.write(gene_name + "," + ",".join(f"{x:.4f}" for x in gc_content) + "\n")
            f_skew.write(gene_name + "," + ",".join(f"{x:.4f}" for x in gc_skew) + "\n")

    print(f"  Done: {fasta_path.name} -> {gc_content_path.name}, {gc_skew_path.name}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch sliding-window GC content and GC skew for all .fa files in a directory"
    )
    parser.add_argument("--input_dir",   required=True,           help="Directory containing input .fa files")
    parser.add_argument("--resolution",  required=True, type=int, help="Sliding-window size (bp)")
    parser.add_argument("--output_path", required=True,           help="Directory to write output CSV files")
    args = parser.parse_args()

    input_dir  = Path(args.input_dir)
    output_dir = Path(args.output_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    fa_files = sorted(input_dir.glob("*.fa"))
    if not fa_files:
        print(f"No .fa files found in {input_dir}")
        return

    print(f"Found {len(fa_files)} .fa file(s) in {input_dir}")
    for fasta_path in fa_files:
        process_fasta(fasta_path, args.resolution, output_dir)

    print("All jobs finished!")


if __name__ == "__main__":
    main()


# Example:
# python 3_3_calGCInfo.py \
#   --input_dir  /path/to/fasta_dir \
#   --resolution 100 \
#   --output_path /path/to/output_dir
