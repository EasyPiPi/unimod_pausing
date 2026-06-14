import os
import argparse

from Bio import SeqIO
from pathlib import Path
import csv

WINDOW_LEN = 200

def main():
    parser = argparse.ArgumentParser(description="Batch-compute GC content and GC skew for all .fa files in a directory")
    parser.add_argument("--input_dir", required=True, help="Directory containing input .fa files")
    parser.add_argument("--output_path", required=True, help="Directory to write output CSV files")

    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    fa_files = sorted(input_dir.glob("*.fa"))
    if not fa_files:
        print(f"No .fa files found in {input_dir}")
        return

    for fasta_path in fa_files:
        stem = fasta_path.stem
        out_file = output_path / f"{stem}_{WINDOW_LEN}_gc_info.csv"

        with open(out_file, "w") as f_gc:
            for record in SeqIO.parse(fasta_path, "fasta"):
                sequence = str(record.seq).upper()
                gene_name = record.description.split("::", 1)[0]

                tss = len(sequence) // 2
                start = tss
                end = min(tss + WINDOW_LEN, len(sequence))

                sub_sequence = sequence[start:end]
                g = sub_sequence.count("G")
                c = sub_sequence.count("C")

                gc_content = (g + c) / WINDOW_LEN
                gc_skew = round((g - c) / (g + c), 4) if (g + c) > 0 else 0

                f_gc.write(f"{gene_name},{gc_content},{gc_skew}\n")

        print(f"Done: {fasta_path.name} -> {out_file.name}")

    print("All jobs finished!")


if __name__ == "__main__":
    main()

# Example:
# python 3_2_calGCInfo_region.py \
#   --input_dir /Users/zeng/Desktop/project/YiXin_Likelihood/outputs/publish/betweenSpecies/sequence \
#   --output_path /Users/zeng/Desktop/project/YiXin_Likelihood/outputs/publish/betweenSpecies/sequence
