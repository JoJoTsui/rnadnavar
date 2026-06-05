import os
import re
import glob
import pandas as pd
import numpy as np


def merge_csv_files(folder):
    """
    Merge all CSV files in a folder into a single DataFrame.
    Extract the last number in each filename as the sample ID and insert
    it as the first column.
    """
    all_files = glob.glob(os.path.join(folder, "*.csv"))
    if not all_files:
        raise FileNotFoundError(f"No CSV files found in folder: '{folder}'")

    df_list = []

    for file in all_files:
        df = pd.read_csv(file, low_memory=False)

        basename = os.path.basename(file)
        stem, _ = os.path.splitext(basename)

        # Extract all numbers from the filename and use the last one as sample ID
        numbers = re.findall(r'\d+', stem)
        if not numbers:
            raise ValueError(f"Cannot extract sample number from filename: '{basename}'")
        sample = numbers[-1]  # take the last group of digits as sample

        df.insert(0, "Sample", sample)
        df_list.append(df)

    merged_df = pd.concat(df_list, ignore_index=True)
    print(f"Merged {len(all_files)} files, total {len(merged_df)} variants")
    return merged_df


def process_vcf_data(input_df, output_csv_path):
    """
    Process the merged VCF data:
    - Split AD columns into REF and ALT allele depths.
    - Calculate VAF (Variant Allele Frequency) for each caller.
    - Compute mean DP, REF, ALT, VAF across DNA tools and RNA tools.
    - Save the processed DataFrame to a CSV file.
    """
    df = input_df.copy()

    # Define which columns belong to DNA/RNA callers
    ad_columns = [
        "DNA_deepsomatic_AD", "DNA_mutect2_AD", "DNA_strelka_AD",
        "RNA_deepsomatic_AD", "RNA_mutect2_AD", "RNA_strelka_AD"
    ]
    dp_columns = [
        "DNA_deepsomatic_DP", "DNA_mutect2_DP", "DNA_strelka_DP",
        "RNA_deepsomatic_DP", "RNA_mutect2_DP", "RNA_strelka_DP"
    ]

    # -- Step 1: split AD, calculate VAF, handling missing or malformed data --
    for ad_col, dp_col in zip(ad_columns, dp_columns):
        if ad_col not in df.columns:
            print(f"Warning: Column '{ad_col}' not found. Skipping its processing.")
            continue  # do not create derived columns; they will be omitted in means

        # Convert AD column to string, split by comma into two parts (REF, ALT)
        ad_split = df[ad_col].astype(str).str.split(',', expand=True)

        # Ensure we have exactly two columns; if more or fewer, pad or trim
        if ad_split.shape[1] < 2:
            # not enough values – fill missing with NaN
            for extra in range(2 - ad_split.shape[1]):
                ad_split[ad_split.shape[1]] = np.nan
        elif ad_split.shape[1] > 2:
            # too many values – keep only first two
            ad_split = ad_split.iloc[:, :2]

        ad_split.columns = ['REF', 'ALT']

        # Convert to numeric, coercing errors to NaN
        ref_vals = pd.to_numeric(ad_split['REF'], errors='coerce')
        alt_vals = pd.to_numeric(ad_split['ALT'], errors='coerce')

        df[f"{ad_col}_REF"] = ref_vals
        df[f"{ad_col}_ALT"] = alt_vals

        # Get DP column; if it does not exist, VAF cannot be calculated
        if dp_col not in df.columns:
            print(f"Warning: DP column '{dp_col}' missing, cannot compute VAF for {ad_col}")
            continue

        dp_vals = pd.to_numeric(df[dp_col], errors='coerce')

        # Calculate VAF: ALT / DP, set to NaN if DP is 0 or missing
        with np.errstate(divide='ignore', invalid='ignore'):
            vaf = alt_vals / dp_vals
        vaf[dp_vals == 0] = np.nan
        vaf[dp_vals.isna()] = np.nan
        df[f"{ad_col}_VAF"] = vaf

    # -- Step 2: compute mean values for DNA callers --
    dna_tools = ["deepsomatic", "mutect2", "strelka"]

    def safe_mean_columns(df, prefix_list, suffix):
        """Compute row-wise mean for columns that exist; skip missing ones."""
        cols = [f"{p}{suffix}" for p in prefix_list]
        existing = [c for c in cols if c in df.columns]
        if not existing:
            return pd.Series(np.nan, index=df.index)
        return df[existing].mean(axis=1)

    # DNA mean columns
    dna_prefixes = [f"DNA_{tool}_" for tool in dna_tools]
    df["DNA_DP_mean"] = safe_mean_columns(df, dna_prefixes, "DP")
    df["DNA_REF_mean"] = safe_mean_columns(df, dna_prefixes, "AD_REF")
    df["DNA_ALT_mean"] = safe_mean_columns(df, dna_prefixes, "AD_ALT")
    df["DNA_VAF_mean"] = safe_mean_columns(df, dna_prefixes, "AD_VAF")

    # RNA mean columns
    rna_prefixes = [f"RNA_{tool}_" for tool in dna_tools]  # same tool names
    df["RNA_DP_mean"] = safe_mean_columns(df, rna_prefixes, "DP")
    df["RNA_REF_mean"] = safe_mean_columns(df, rna_prefixes, "AD_REF")
    df["RNA_ALT_mean"] = safe_mean_columns(df, rna_prefixes, "AD_ALT")
    df["RNA_VAF_mean"] = safe_mean_columns(df, rna_prefixes, "AD_VAF")

    # Save results
    df.to_csv(output_csv_path, index=False)
    print(f"Processed data saved to {output_csv_path}")
    # The function intentionally returns None; avoid assigning it to a variable


if __name__ == "__main__":
    input_folder = r"/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/variant_info"
    output_file = r"/t9k/mnt/WorkSpace/data/ngs/zhanlingmin/variant_info/set13_variant_info.csv"
    merged_df = merge_csv_files(input_folder)
    process_vcf_data(merged_df, output_file)  # no need to capture return value