#!/usr/bin/env python3
"""
GE2TS - Variant Matching Script
Matches GWAS summary statistics variants against 1000G PLINK reference .bim file,
mapping variant IDs to reference rsIDs.
"""

import argparse
import os
import sys
import gzip
import csv
try:
    from bin.engine import get_engine
except ModuleNotFoundError:
    from engine import get_engine

def load_reference_maps(bim_file):
    print(f"Loading reference BIM file: {bim_file}")
    engine = get_engine()
    return engine.load_bim_and_create_maps(bim_file)

def match_variant_row(row, id_map, exact_map, swap_map, strand_map):
    chr_val = str(row.get('CHR', ''))
    bp_val = str(row.get('BP', ''))
    a1_val = str(row.get('A1', '')).upper()
    a2_val = str(row.get('A2', '')).upper()
    orig_id = row.get('SNP')
    
    k = f"{chr_val}:{bp_val}:{a1_val}:{a2_val}"
    
    if orig_id and id_map.get(orig_id):
        row['SNP'] = id_map.get(orig_id)
        return row, True
    elif exact_map.get(k):
        row['SNP'] = exact_map.get(k)
        return row, True
    elif swap_map.get(k):
        row['SNP'] = swap_map.get(k)
        return row, True
    elif strand_map.get(k):
        row['SNP'] = strand_map.get(k)
        return row, True
    
    return row, False

def stream_and_match_gwas(gwas_file, output_file, maps):
    id_map, exact_map, swap_map, strand_map = maps
    print(f"Processing GWAS file: {gwas_file}")
    open_func = gzip.open if gwas_file.endswith('.gz') else open
    
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    matched_count = 0
    total_count = 0

    with open_func(gwas_file, 'rt') as f_in, open(output_file, 'wt') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        fieldnames = reader.fieldnames
        if not fieldnames:
            return matched_count, total_count
            
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for row in reader:
            total_count += 1
            updated_row, is_matched = match_variant_row(row, id_map, exact_map, swap_map, strand_map)
            if is_matched:
                matched_count += 1
            writer.writerow(updated_row)

    return matched_count, total_count

def match_variants(gwas_file, output_file, bim_file):
    maps = load_reference_maps(bim_file)
    matched, total = stream_and_match_gwas(gwas_file, output_file, maps)
    print(f"Variant matching complete. Matched {matched} / {total} variants to reference BIM IDs.")

def main():
    parser = argparse.ArgumentParser(description="GE2TS GWAS Variant Matching")
    parser.add_argument("--gwas", required=True, help="Input standardized GWAS file")
    parser.add_argument("--bim", required=True, help="Path to PLINK .bim reference file")
    parser.add_argument("--output-tsv", required=True, help="Output matched GWAS file")
    
    args = parser.parse_args()
    match_variants(args.gwas, args.output_tsv, args.bim)

if __name__ == "__main__":
    main()
