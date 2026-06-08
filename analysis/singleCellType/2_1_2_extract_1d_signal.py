import pandas as pd
import numpy as np
import pyranges as pr
import argparse 

def load_pairs_filter_mapq(pairs_path):
    """
    Read a .pairs or .pairs.gz file and keep rows with mapq1>=30 and mapq2>=30.
    """
    cols = ["readID","chrom1","pos1","chrom2","pos2","strand1","strand2","pair_type","mapq1","mapq2"]
    df = pd.read_table(
        pairs_path,
        compression="gzip" if pairs_path.endswith(".gz") else None,
        comment="#",
        header=None,
        names=cols,
        dtype={"chrom1": str, "chrom2": str, "strand1": str, "strand2": str}
    )
    df = df[(df["mapq1"] >= 30) & (df["mapq2"] >= 30)]
    return df

def to_1d_shift_by_strand(df, shift_bp=75):
    """
    Collapse both ends into a 1D table and shift downstream by 'shift_bp' on each end
    according to its own strand. Clamp negatives to 1.
    """
    left = df[["chrom1","pos1","strand1"]].rename(
        columns={"chrom1":"chrom", "pos1":"pos", "strand1":"strand"}
    )
    right = df[["chrom2","pos2","strand2"]].rename(
        columns={"chrom2":"chrom", "pos2":"pos", "strand2":"strand"}
    )
    ends = pd.concat([left, right], ignore_index=True)

    # Shift downstream by strand
    ends["pos_shift"] = (ends["pos"] + np.where(ends["strand"].eq("+"), shift_bp, -shift_bp))
    ends["pos_shift"] = ends["pos_shift"].clip(lower=1).astype(int)

    # Convert to 0-based half-open intervals for PyRanges (single-base interval)
    ends["Start"] = ends["pos_shift"] - 1
    ends["End"] = ends["pos_shift"]
    ends = ends[["chrom","Start","End","strand"]].rename(columns={"chrom":"Chromosome"})
    counts = (
        ends.groupby(["Chromosome","Start","End"])
        .size()
        .reset_index(name="Score")
        )
    count_pr = pr.PyRanges(counts)
    return count_pr

def load_tss_and_make_windows(tss_bed, flank=1000):
    """
    Load a TSS bed and build +/- flank windows around each TSS.

    If assume_point=True, treat the TSS as the midpoint of each BED row
    (works for point-like TSS beds or narrow intervals). Otherwise, keep
    original intervals and just extend both sides by 'flank'.
    """
    tss = pd.read_csv(tss_bed, sep='\t', header=None)  # columns: Chromosome, Start, End, (Name), (Score), (Strand) ...

    tss.columns = ['Chromosome', 'tss_start', 'tss_end', 'Name', 'gCount', 'Strand', 'pCount']
    center = tss["tss_start"].astype(int)
    tss["Start"] = (center - flank).clip(lower=0)
    tss["End"] = center + flank
    keep_cols = ['Chromosome', 'Start', 'End', 'tss_start', 'Name']
    return pr.PyRanges(tss[keep_cols])




def main():
    parser = argparse.ArgumentParser(description="Extract 1D signal from a pair.gz file")
    parser.add_argument("--pair_path", required=True, help="Path to pair.gz file.")
    parser.add_argument("--tss_file", required=True, help="Path to tss file.")
    parser.add_argument("--output_file", required=True, help="Path to tss file.")
    args = parser.parse_args()

    pairs_df = load_pairs_filter_mapq(args.pair_path)
    shifted_pr = to_1d_shift_by_strand(pairs_df, shift_bp=75)
    tss_pr = load_tss_and_make_windows(args.tss_file)

    shifted_tss_pr = shifted_pr.join(tss_pr, how=None)
    cols_to_keep = ["Chromosome", "Start", "End", "Score", "Name", "tss_start"]
    shifted_tss_pr[cols_to_keep].to_bed(args.output_file)
    print(shifted_tss_pr)
if __name__ == '__main__':
    main()

#python extract_1d_signal.py \
#--pair_path /grid/siepel/home/zeng/projects/Zhao_2025/comparaReg_analyses/resources/contacts/pairs_files/human_cd4_merged.pairs.gz \
# --tss_file  /grid/siepel/home/zeng/projects/Zhao_2025/result/STADyUM/cd4_tss_sort.bed \
#--output_file /grid/siepel/home/zeng/projects/Zhao_2025/result/STADyUM/cd4_tss_1d_signal.bed




