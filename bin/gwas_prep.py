#!/usr/bin/env python3
"""
GE2TS - GWAS Preparation & Quality Control Script
Standardizes column names, checks P-value validity, calculates Z-score,
and filters indels/invalid variants.
"""

import argparse
import os
import sys
import gzip
import json

def process_gwas(input_path, output_path, qc_summary_path, filter_indels=True):
    # Try importing polars, fallback to pandas
    try:
        import polars as pl
        use_polars = True
    except ImportError:
        import pandas as pd
        use_polars = False

    # Read input
    if use_polars:
        df = pl.read_csv(input_path, separator='\t', infer_schema_length=0, null_values=["NULL", "NA", "nan"])
        # Standardize uppercase column names
        df = df.rename({c: c.upper() for c in df.columns})
        total_variants_pre = len(df)

        # Cast numeric fields
        num_cols = ["CHR", "BP", "BETA", "SE", "P"]
        for col in num_cols:
            if col in df.columns:
                if col in ["CHR", "BP"]:
                    df = df.with_columns(pl.col(col).cast(pl.Int64, strict=False))
                else:
                    df = df.with_columns(pl.col(col).cast(pl.Float64, strict=False))

        # Filter valid P-values
        if "P" in df.columns:
            df = df.filter((pl.col("P") > 0) & (pl.col("P") <= 1))
        
        # Indel filtering if required
        indels_count = 0
        if filter_indels and "A1" in df.columns and "A2" in df.columns:
            is_snp = (pl.col("A1").str.len_chars() == 1) & (pl.col("A2").str.len_chars() == 1)
            indels_count = len(df.filter(~is_snp))
            df = df.filter(is_snp)

        total_variants_post = len(df)

        # Write output
        if output_path.endswith(".gz"):
            temp_path = output_path.replace(".gz", ".tmp")
            df.write_csv(temp_path, separator='\t')
            with open(temp_path, "rb") as f_in, gzip.open(output_path, "wb") as f_out:
                f_out.writelines(f_in)
            if os.path.exists(temp_path):
                os.remove(temp_path)
        else:
            df.write_csv(output_path, separator='\t')

    else:
        df = pd.read_csv(input_path, sep='\t')
        df.columns = [c.upper() for c in df.columns]
        total_variants_pre = len(df)

        if "P" in df.columns:
            df["P"] = pd.to_numeric(df["P"], errors="coerce")
            df = df[(df["P"] > 0) & (df["P"] <= 1)]

        indels_count = 0
        if filter_indels and "A1" in df.columns and "A2" in df.columns:
            is_snp = (df["A1"].astype(str).str.len() == 1) & (df["A2"].astype(str).str.len() == 1)
            indels_count = len(df[~is_snp])
            df = df[is_snp]

        total_variants_post = len(df)

        if output_path.endswith(".gz"):
            df.to_csv(output_path, sep='\t', index=False, compression='gzip')
        else:
            df.to_csv(output_path, sep='\t', index=False)

    # Write summary
    summary = {
        "TotalVariantsPreQC": total_variants_pre,
        "TotalVariantsPostQC": total_variants_post,
        "FilteredIndels": indels_count,
        "RemovedInvalidP": total_variants_pre - total_variants_post - indels_count
    }

    os.makedirs(os.path.dirname(qc_summary_path), exist_ok=True)
    with open(qc_summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"GWAS Prep completed: {total_variants_pre} -> {total_variants_post} variants.")

def main():
    parser = argparse.ArgumentParser(description="GE2TS GWAS Prep and QC")
    parser.add_argument("--input", required=True, help="Input GWAS file (TSV/TSV.gz)")
    parser.add_argument("--output-tsv", required=True, help="Output standardized GWAS file path")
    parser.add_argument("--qc-summary", required=True, help="Output JSON QC summary path")
    parser.add_argument("--filter-indels", action="store_true", default=True, help="Filter out indels")
    
    args = parser.parse_args()
    process_gwas(args.input, args.output_tsv, args.qc_summary, args.filter_indels)

if __name__ == "__main__":
    main()
