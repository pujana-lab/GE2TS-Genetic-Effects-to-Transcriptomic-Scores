#!/usr/bin/env python3
"""
GE2TS - QC Audit, Transcriptomic Scoring & GSEA Pathway Analysis Script
Processes MAGMA gene scores, evaluates pathway enrichment against GMT databases,
and outputs transcriptomic scores, QC summary metrics, and an interactive HTML report.
"""

import argparse
import os
import sys
import numpy as np
from scipy.stats import norm
import pandas as pd

try:
    from bin.engine import get_engine
except ModuleNotFoundError:
    from engine import get_engine

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

    engine = get_engine()
    genes_df = engine.read_genes_out(genes_out_path)

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

def load_html_template(template_path=None):
    default_template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates", "gsea_report.html")
    candidate_paths = [template_path, default_template_path]
    
    for candidate in candidate_paths:
        if candidate and os.path.exists(candidate):
            with open(candidate, "r", encoding="utf-8") as f:
                return f.read()

    # Minimal fallback template if file is missing
    return """<!DOCTYPE html><html><head><title>{{PHENOTYPE}}</title></head><body><h1>{{PHENOTYPE}}</h1><div>{{TABLE_ROWS}}</div></body></html>"""

def generate_gsea_html_report(phenotype, genes_df, pathway_df, html_out_path, template_path=None):
    total_genes = len(genes_df)
    sig_genes = len(genes_df[genes_df["P"] < 0.05]) if "P" in genes_df.columns else 0
    num_pathways = len(pathway_df) if not pathway_df.empty else 0
    top_pathway = pathway_df.iloc[0]["PATHWAY"] if not pathway_df.empty else "N/A"
    top_p = pathway_df.iloc[0]["P_VALUE"] if not pathway_df.empty else 1.0

    if not pathway_df.empty:
        table_rows = "\n".join([
            f"""
            <tr>
                <td><strong>{row['PATHWAY']}</strong></td>
                <td>{row['N_GENES']}</td>
                <td>{row['N_MATCHED']}</td>
                <td>{row['MEAN_Z']}</td>
                <td>{row['MEDIAN_Z']}</td>
                <td>{row['MAX_Z']}</td>
                <td><span class="gene-badge">{row['TOP_GENE']}</span></td>
                <td><code>{row['P_VALUE']}</code></td>
            </tr>
            """
            for _, row in pathway_df.iterrows()
        ])
    else:
        table_rows = "<tr><td colspan='8'>No pathway enrichment results available.</td></tr>"

    raw_template = load_html_template(template_path)
    html_content = (
        raw_template
        .replace("{{PHENOTYPE}}", str(phenotype))
        .replace("{{TOTAL_GENES}}", str(total_genes))
        .replace("{{SIG_GENES}}", str(sig_genes))
        .replace("{{NUM_PATHWAYS}}", str(num_pathways))
        .replace("{{TOP_PATHWAY}}", str(top_pathway))
        .replace("{{TOP_P}}", str(top_p))
        .replace("{{TABLE_ROWS}}", table_rows)
    )

    out_dir = os.path.dirname(html_out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(html_out_path, "w", encoding="utf-8") as f:
        f.write(html_content)

def run_gsea_and_qc(genes_out_path, gmt_path, out_dir, qc_summary_path, gene_loc=None, template_path=None):
    os.makedirs(out_dir, exist_ok=True)
    phenotype = os.path.basename(genes_out_path).replace(".genes.out", "")

    genes_df = load_genes_data(genes_out_path, gene_loc)
    pathways = parse_gmt(gmt_path)
    pathway_df = calculate_pathway_zscores(genes_df, pathways)

    gsea_tsv = os.path.join(out_dir, f"{phenotype}_gsea_pathways.tsv")
    pathway_df.to_csv(gsea_tsv, sep='\t', index=False)

    html_out = os.path.join(out_dir, f"{phenotype}_gsea_report.html")
    generate_gsea_html_report(phenotype, genes_df, pathway_df, html_out, template_path=template_path)

    write_qc_summary(genes_df, pathway_df, phenotype, qc_summary_path)

    print(f"GSEA complete: Saved pathway scores to {gsea_tsv}, HTML report to {html_out}, and QC matrix to {qc_summary_path}.")

def main():
    parser = argparse.ArgumentParser(description="GE2TS GSEA and QC Audit")
    parser.add_argument("--genes-out", required=True, help="Path to MAGMA .genes.out file")
    parser.add_argument("--pathway-gmt", required=True, help="Path to GMT pathway database")
    parser.add_argument("--out-dir", required=True, help="Output directory for GSEA results")
    parser.add_argument("--qc-summary-tsv", required=True, help="Output path for QC summary TSV")
    parser.add_argument("--gene-loc", help="Path to Gene location file for SYMBOL mapping")
    parser.add_argument("--template", help="Path to custom HTML report template file")

    args = parser.parse_args()
    run_gsea_and_qc(
        genes_out_path=args.genes_out,
        gmt_path=args.pathway_gmt,
        out_dir=args.out_dir,
        qc_summary_path=args.qc_summary_tsv,
        gene_loc=args.gene_loc,
        template_path=args.template
    )

if __name__ == "__main__":
    main()
