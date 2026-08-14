#!/usr/bin/env python3
"""
GE2TS - QC Audit, Transcriptomic Scoring & GSEA Pathway Analysis Script
Processes MAGMA gene scores, evaluates pathway enrichment against GMT databases,
and outputs transcriptomic scores and QC summary metrics.
"""

import argparse
import os
import sys
import numpy as np
from scipy.stats import norm

try:
    import polars as pl
    USE_POLARS = True
except ImportError:
    USE_POLARS = False
import pandas as pd

def parse_gmt(gmt_path):
    pathways = {}
    if not os.path.exists(gmt_path):
        return pathways

    with open(gmt_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                name = parts[0]
                genes = [g.strip().upper() for g in parts[2:] if g.strip()]
                pathways[name] = set(genes)
    return pathways

def load_genes_data(genes_out_path, gene_loc=None):
    if not os.path.exists(genes_out_path):
        raise FileNotFoundError(f"MAGMA genes file not found: {genes_out_path}")

    if USE_POLARS:
        try:
            genes_df = pl.read_csv(genes_out_path, separator='\t').to_pandas()
        except Exception:
            genes_df = pd.read_csv(genes_out_path, sep=r'\s+', engine='python')
    else:
        try:
            genes_df = pd.read_csv(genes_out_path, sep=r'\s+', engine='python')
        except Exception:
            genes_df = pd.read_csv(genes_out_path, sep='\t')

    genes_df.columns = [c.upper() for c in genes_df.columns]
    
    if "SYMBOL" not in genes_df.columns and gene_loc and os.path.exists(gene_loc):
        try:
            loc_df = pd.read_csv(gene_loc, sep=r'\s+', header=None)
            mapping = dict(zip(loc_df[0], loc_df[loc_df.columns[-1]]))
            genes_df["SYMBOL"] = genes_df["GENE"].map(mapping)
            print(f"[QC GSEA] Added SYMBOL column using {gene_loc}")
        except Exception as e:
            print(f"[QC GSEA] Failed to add SYMBOL column using {gene_loc}: {e}")
            
    if "SYMBOL" in genes_df.columns:
        genes_df["SYMBOL"] = genes_df["SYMBOL"].astype(str).str.upper()

    return genes_df

def calculate_pathway_zscores(genes_df, pathways):
    pathway_results = []
    
    for p_name, p_genes in pathways.items():
        matched = genes_df[genes_df["SYMBOL"].isin(p_genes)] if "SYMBOL" in genes_df.columns else pd.DataFrame()
        
        n_matched = len(matched)
        if n_matched > 0 and "ZSTAT" in matched.columns:
            mean_z = float(matched["ZSTAT"].mean())
            median_z = float(matched["ZSTAT"].median())
            max_z = float(matched["ZSTAT"].max())
            top_gene = matched.loc[matched["ZSTAT"].idxmax(), "SYMBOL"] if "SYMBOL" in matched.columns else "N/A"
            combined_z = mean_z * np.sqrt(n_matched)
            p_val = float(2 * norm.sf(abs(combined_z)))
        else:
            mean_z, median_z, max_z, top_gene, p_val = 0.0, 0.0, 0.0, "N/A", 1.0

        pathway_results.append({
            "PATHWAY": p_name,
            "N_GENES": len(p_genes),
            "N_MATCHED": n_matched,
            "MEAN_Z": round(mean_z, 4),
            "MEDIAN_Z": round(median_z, 4),
            "MAX_Z": round(max_z, 4),
            "TOP_GENE": top_gene,
            "P_VALUE": round(p_val, 6)
        })

    pathway_df = pd.DataFrame(pathway_results)
    if not pathway_df.empty:
        pathway_df = pathway_df.sort_values(by=["MEAN_Z", "P_VALUE"], ascending=[False, True])
    return pathway_df

def write_qc_summary(genes_df, pathway_df, phenotype, qc_summary_path):
    total_genes = len(genes_df)
    significant_genes = len(genes_df[genes_df["P"] < 0.05]) if "P" in genes_df.columns else 0

    qc_data = [{
        "Phenotype": phenotype,
        "TotalTestedGenes": total_genes,
        "SignificantGenes_P005": significant_genes,
        "PathwaysEvaluated": len(pathway_df) if not pathway_df.empty else 0,
        "TopEnrichedPathway": pathway_df.iloc[0]["PATHWAY"] if not pathway_df.empty else "N/A",
        "TopPathway_PValue": pathway_df.iloc[0]["P_VALUE"] if not pathway_df.empty else 1.0
    }]

    qc_df = pd.DataFrame(qc_data)
    qc_out_dir = os.path.dirname(qc_summary_path)
    if qc_out_dir:
        os.makedirs(qc_out_dir, exist_ok=True)

    qc_df.to_csv(qc_summary_path, sep='\t', index=False)

def run_gsea_and_qc(genes_out_path, gmt_path, out_dir, qc_summary_path, gene_loc=None):
    os.makedirs(out_dir, exist_ok=True)
    phenotype = os.path.basename(genes_out_path).replace(".genes.out", "")

    genes_df = load_genes_data(genes_out_path, gene_loc)
    pathways = parse_gmt(gmt_path)
    pathway_df = calculate_pathway_zscores(genes_df, pathways)

    gsea_tsv = os.path.join(out_dir, f"{phenotype}_gsea_pathways.tsv")
    pathway_df.to_csv(gsea_tsv, sep='\t', index=False)

    write_qc_summary(genes_df, pathway_df, phenotype, qc_summary_path)

    print(f"GSEA complete: Saved pathway scores to {gsea_tsv} and QC matrix to {qc_summary_path}.")

def main():
    parser = argparse.ArgumentParser(description="GE2TS GSEA and QC Audit")
    parser.add_argument("--genes-out", required=True, help="Path to MAGMA .genes.out file")
    parser.add_argument("--pathway-gmt", required=True, help="Path to GMT pathway database")
    parser.add_argument("--out-dir", required=True, help="Output directory for GSEA results")
    parser.add_argument("--qc-summary-tsv", required=True, help="Output path for QC summary TSV")
    parser.add_argument("--gene-loc", help="Path to Gene location file for SYMBOL mapping")

    args = parser.parse_args()
    run_gsea_and_qc(
        genes_out_path=args.genes_out,
        gmt_path=args.pathway_gmt,
        out_dir=args.out_dir,
        qc_summary_path=args.qc_summary_tsv,
        gene_loc=args.gene_loc
    )

if __name__ == "__main__":
    main()
