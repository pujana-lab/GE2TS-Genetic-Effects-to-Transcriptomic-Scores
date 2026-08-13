#!/usr/bin/env python3
"""
GE2TS - Variant Matching Script
Matches GWAS summary statistics variants against 1000G PLINK reference .bim file,
mapping variant IDs to reference rsIDs.
"""

import argparse
import os
import sys
import csv
import gzip
import polars as pl

def match_variants(gwas_file, output_file, bim_file):
    print(f"Loading reference BIM file: {bim_file}")
    ref = pl.read_csv(
        bim_file,
        separator='\t',
        has_header=False,
        new_columns=['CHR', 'ID', 'CM', 'POS', 'A1', 'A2'],
        schema_overrides={'CHR': pl.String, 'POS': pl.Int64, 'ID': pl.String, 'A1': pl.String, 'A2': pl.String}
    )
    
    comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    ref = ref.with_columns([
        pl.concat_str([pl.col("CHR"), pl.col("POS"), pl.col("A1"), pl.col("A2")], separator=":").alias("key_exact"),
        pl.concat_str([pl.col("CHR"), pl.col("POS"), pl.col("A2"), pl.col("A1")], separator=":").alias("key_swapped"),
        pl.col("A1").map_elements(lambda x: "".join([comp_map.get(b, b) for b in x]), return_dtype=pl.String).alias("A1_c"),
        pl.col("A2").map_elements(lambda x: "".join([comp_map.get(b, b) for b in x]), return_dtype=pl.String).alias("A2_c")
    ]).with_columns(
        pl.concat_str([pl.col("CHR"), pl.col("POS"), pl.col("A1_c"), pl.col("A2_c")], separator=":").alias("key_strand")
    )

    id_map = dict(zip(ref['ID'], ref['ID']))
    exact_map = dict(zip(ref['key_exact'], ref['ID']))
    swap_map = dict(zip(ref['key_swapped'], ref['ID']))
    strand_map = dict(zip(ref['key_strand'], ref['ID']))
    
    print(f"Processing GWAS file: {gwas_file}")
    open_func = gzip.open if gwas_file.endswith('.gz') else open
    
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open_func(gwas_file, 'rt') as f_in, open(output_file, 'wt') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        fieldnames = reader.fieldnames
        if not fieldnames:
            return
            
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        matched_count = 0
        total_count = 0
        
        for row in reader:
            total_count += 1
            chr_val = str(row.get('CHR', ''))
            bp_val = str(row.get('BP', ''))
            a1_val = str(row.get('A1', '')).upper()
            a2_val = str(row.get('A2', '')).upper()
            orig_id = row.get('SNP')
            
            k = f"{chr_val}:{bp_val}:{a1_val}:{a2_val}"
            
            if orig_id and id_map.get(orig_id):
                row['SNP'] = id_map.get(orig_id)
                matched_count += 1
            elif exact_map.get(k):
                row['SNP'] = exact_map.get(k)
                matched_count += 1
            elif swap_map.get(k):
                row['SNP'] = swap_map.get(k)
                matched_count += 1
            elif strand_map.get(k):
                row['SNP'] = strand_map.get(k)
                matched_count += 1
            else:
                pass
            
            writer.writerow(row)

    print(f"Variant matching complete. Matched {matched_count} / {total_count} variants to reference BIM IDs.")

def main():
    parser = argparse.ArgumentParser(description="GE2TS GWAS Variant Matching")
    parser.add_argument("--gwas", required=True, help="Input standardized GWAS file")
    parser.add_argument("--bim", required=True, help="Path to PLINK .bim reference file")
    parser.add_argument("--output-tsv", required=True, help="Output matched GWAS file")
    
    args = parser.parse_args()
    match_variants(args.gwas, args.output_tsv, args.bim)

if __name__ == "__main__":
    main()
