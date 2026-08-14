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
    from bin.engine import get_engine, HAS_POLARS
except ModuleNotFoundError:
    from engine import get_engine, HAS_POLARS

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
    if engine.engine == 'polars':
        import polars as pl
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

def generate_gsea_html_report(phenotype, genes_df, pathway_df, html_out_path):
    total_genes = len(genes_df)
    sig_genes = len(genes_df[genes_df["P"] < 0.05]) if "P" in genes_df.columns else 0
    num_pathways = len(pathway_df) if not pathway_df.empty else 0
    top_pathway = pathway_df.iloc[0]["PATHWAY"] if not pathway_df.empty else "N/A"
    top_p = pathway_df.iloc[0]["P_VALUE"] if not pathway_df.empty else 1.0

    table_rows = ""
    if not pathway_df.empty:
        for _, row in pathway_df.iterrows():
            table_rows += f"""
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
    else:
        table_rows = "<tr><td colspan='8'>No pathway enrichment results available.</td></tr>"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GE2TS GSEA Enrichment Report - {phenotype}</title>
    <style>
        :root {{
            --primary: #2563eb;
            --primary-dark: #1d4ed8;
            --bg-dark: #0f172a;
            --card-bg: #1e293b;
            --text-light: #f8fafc;
            --text-muted: #94a3b8;
            --accent: #10b981;
            --border: #334155;
        }}
        body {{
            font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background-color: var(--bg-dark);
            color: var(--text-light);
            margin: 0;
            padding: 2rem;
            line-height: 1.5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        header {{
            margin-bottom: 2rem;
            border-bottom: 1px solid var(--border);
            padding-bottom: 1rem;
        }}
        h1 {{
            margin: 0 0 0.5rem 0;
            font-size: 2rem;
            color: var(--text-light);
        }}
        .subtitle {{
            color: var(--text-muted);
            margin: 0;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1rem;
            margin-bottom: 2rem;
        }}
        .stat-card {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 0.5rem;
            padding: 1.25rem;
        }}
        .stat-card .label {{
            font-size: 0.875rem;
            color: var(--text-muted);
            margin-bottom: 0.5rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}
        .stat-card .value {{
            font-size: 1.75rem;
            font-weight: 700;
            color: var(--accent);
        }}
        .table-container {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 0.5rem;
            overflow-x: auto;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            text-align: left;
            font-size: 0.95rem;
        }}
        th {{
            background-color: rgba(15, 23, 42, 0.6);
            color: var(--text-muted);
            font-weight: 600;
            padding: 1rem;
            border-bottom: 1px solid var(--border);
            text-transform: uppercase;
            font-size: 0.8rem;
            letter-spacing: 0.05em;
        }}
        td {{
            padding: 1rem;
            border-bottom: 1px solid var(--border);
        }}
        tr:hover {{
            background-color: rgba(51, 65, 85, 0.4);
        }}
        .gene-badge {{
            background-color: rgba(37, 99, 235, 0.2);
            color: #60a5fa;
            padding: 0.2rem 0.5rem;
            border-radius: 0.25rem;
            font-family: monospace;
            font-size: 0.875rem;
        }}
        code {{
            font-family: monospace;
            color: #f43f5e;
        }}
        footer {{
            margin-top: 3rem;
            text-align: center;
            color: var(--text-muted);
            font-size: 0.875rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 GE2TS GSEA Enrichment Report</h1>
            <p class="subtitle">Phenotype: <strong>{phenotype}</strong> | Transcriptomic Pathway Score Analysis</p>
        </header>

        <div class="stats-grid">
            <div class="stat-card">
                <div class="label">Total Tested Genes</div>
                <div class="value">{total_genes}</div>
            </div>
            <div class="stat-card">
                <div class="label">Sig. Genes (P < 0.05)</div>
                <div class="value">{sig_genes}</div>
            </div>
            <div class="stat-card">
                <div class="label">Pathways Evaluated</div>
                <div class="value">{num_pathways}</div>
            </div>
            <div class="stat-card">
                <div class="label">Top Enriched Pathway</div>
                <div class="value" style="font-size: 1.1rem; word-break: break-word;">{top_pathway} (P={top_p})</div>
            </div>
        </div>

        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th>Pathway</th>
                        <th>N Genes</th>
                        <th>N Matched</th>
                        <th>Mean Z</th>
                        <th>Median Z</th>
                        <th>Max Z</th>
                        <th>Top Gene</th>
                        <th>P-Value</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>

        <footer>
            Generated automatically by GE2TS DSL2 Nextflow Pipeline
        </footer>
    </div>
</body>
</html>
"""
    out_dir = os.path.dirname(html_out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(html_out_path, "w", encoding="utf-8") as f:
        f.write(html_content)

def run_gsea_and_qc(genes_out_path, gmt_path, out_dir, qc_summary_path, gene_loc=None):
    os.makedirs(out_dir, exist_ok=True)
    phenotype = os.path.basename(genes_out_path).replace(".genes.out", "")

    genes_df = load_genes_data(genes_out_path, gene_loc)
    pathways = parse_gmt(gmt_path)
    pathway_df = calculate_pathway_zscores(genes_df, pathways)

    gsea_tsv = os.path.join(out_dir, f"{phenotype}_gsea_pathways.tsv")
    pathway_df.to_csv(gsea_tsv, sep='\t', index=False)

    html_out = os.path.join(out_dir, f"{phenotype}_gsea_report.html")
    generate_gsea_html_report(phenotype, genes_df, pathway_df, html_out)

    write_qc_summary(genes_df, pathway_df, phenotype, qc_summary_path)

    print(f"GSEA complete: Saved pathway scores to {gsea_tsv}, HTML report to {html_out}, and QC matrix to {qc_summary_path}.")

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
