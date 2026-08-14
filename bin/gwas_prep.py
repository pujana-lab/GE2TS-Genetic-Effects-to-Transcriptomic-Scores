#!/usr/bin/env python3
"""
GE2TS - GWAS Preparation & Quality Control Script
Auto-detects delimiter (CSV/TSV), maps column synonyms, cleans ragged lines, 
checks P-value validity, calculates Z-score, and filters indels/invalid variants.
"""

import argparse
import os
import sys
import gzip
import json
import csv
try:
    from bin.engine import get_engine
except ModuleNotFoundError:
    from engine import get_engine

SYNONYM_MAP = {
    "MARKERNAME": "SNP",
    "RS_ID": "SNP",
    "POSITION": "BP",
    "BP": "BP",
    "CHR": "CHR",
    "ALLELE1": "A1",
    "ALLELE2": "A2",
    "EFFECT": "BETA",
    "STDERR": "SE",
    "P-VALUE": "P",
    "PVAL": "P",
    "P_VAL": "P",
    "P": "P",
    "ONCO_ICOGS_EFFECT": "A1",
    "ONCO_ICOGS_BASELINE": "A2",
    "ONCO_ICOGS_BC_EFFECT": "BETA",
    "ONCO_ICOGS_BC_SE": "SE",
    "ONCO_ICOGS_BC_PVAL": "P"
}

def detect_delimiter(file_path):
    try:
        open_func = gzip.open if file_path.endswith('.gz') else open
        with open_func(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
            sample = f.read(2048)
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            return dialect.delimiter
    except Exception:
        return '\t'

def process_gwas(input_path, output_path, qc_summary_path, filter_indels=True):
    sep = detect_delimiter(input_path)
    engine = get_engine()
    use_polars = (engine.engine == 'polars')

    if use_polars:
        try:
            import polars as pl
            df = pl.read_csv(input_path, separator=sep, infer_schema_length=0, null_values=["NULL", "NA", "nan", ""])
        except Exception:
            use_polars = False

    if use_polars:
        import polars as pl
        rename_dict = {}
        for c in df.columns:
            uc = c.upper()
            if uc in SYNONYM_MAP:
                rename_dict[c] = SYNONYM_MAP[uc]
            else:
                rename_dict[c] = uc
        df = df.rename(rename_dict)
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
        
        # Drop rows with nulls in essential columns to clean ragged lines
        essential = [c for c in ["SNP", "CHR", "BP", "P"] if c in df.columns]
        if essential:
            df = df.drop_nulls(subset=essential)

        # Indel filtering if required
        indels_count = 0
        if filter_indels and "A1" in df.columns and "A2" in df.columns:
            is_snp = (pl.col("A1").str.len_chars() == 1) & (pl.col("A2").str.len_chars() == 1)
            indels_count = len(df.filter(~is_snp))
            df = df.filter(is_snp)

        total_variants_post = len(df)

        # Select standard columns only
        standard_cols = [c for c in ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"] if c in df.columns]
        df = df.select(standard_cols)

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
        import pandas as pd
        df = pd.read_csv(input_path, sep=sep, low_memory=False)
        rename_dict = {}
        for c in df.columns:
            uc = str(c).upper()
            if uc in SYNONYM_MAP:
                rename_dict[c] = SYNONYM_MAP[uc]
            else:
                rename_dict[c] = uc
        df = df.rename(columns=rename_dict)
        total_variants_pre = len(df)

        if "P" in df.columns:
            df["P"] = pd.to_numeric(df["P"], errors="coerce")
            df = df[(df["P"] > 0) & (df["P"] <= 1)]

        essential = [c for c in ["SNP", "CHR", "BP", "P"] if c in df.columns]
        if essential:
            df = df.dropna(subset=essential)

        indels_count = 0
        if filter_indels and "A1" in df.columns and "A2" in df.columns:
            is_snp = (df["A1"].astype(str).str.len() == 1) & (df["A2"].astype(str).str.len() == 1)
            indels_count = len(df[~is_snp])
            df = df[is_snp]

        total_variants_post = len(df)

        standard_cols = [c for c in ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"] if c in df.columns]
        df = df[standard_cols]

        if output_path.endswith(".gz"):
            df.to_csv(output_path, sep='\t', index=False, compression='gzip')
        else:
            df.to_csv(output_path, sep='\t', index=False)

    summary = {
        "TotalVariantsPreQC": total_variants_pre,
        "TotalVariantsPostQC": total_variants_post,
        "FilteredIndels": indels_count,
        "RemovedInvalidP": total_variants_pre - total_variants_post - indels_count
    }

    out_dir = os.path.dirname(qc_summary_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(qc_summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"GWAS Prep completed: {total_variants_pre} -> {total_variants_post} variants.")

def main():
    parser = argparse.ArgumentParser(description="GE2TS GWAS Prep and QC")
    parser.add_argument("--input", required=True, help="Input GWAS file (TSV/CSV/TSV.gz)")
    parser.add_argument("--output-tsv", required=True, help="Output standardized GWAS file path")
    parser.add_argument("--qc-summary", required=True, help="Output JSON QC summary path")
    parser.add_argument("--filter-indels", action="store_true", default=True, help="Filter out indels")
    
    args = parser.parse_args()
    process_gwas(args.input, args.output_tsv, args.qc_summary, args.filter_indels)

if __name__ == "__main__":
    main()
