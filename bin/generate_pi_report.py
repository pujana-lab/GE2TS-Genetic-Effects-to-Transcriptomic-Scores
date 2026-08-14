#!/usr/bin/env python3
"""
GE2TS - Consolidated Executive GSEA Matrix Report Generator for PI Review
Parses pathway enrichment results across all 9 execution folders (3 Subtypes x 3 Modes)
and produces a 50 Hallmark Pathways x 9 Datasets matrix HTML dashboard.
Differentially expressed / enriched pathways (P < 0.05) are highlighted in RED.
"""

import glob
import os
import datetime
import pandas as pd

def generate_consolidated_report(results_dir="/home/luis/Escritorio/gwas_to_gsea", output_html="/home/luis/Escritorio/gwas_to_gsea/ge2ts_cimba_gsea_summary_report.html"):
    runs = [
        ('BRCA1 10kb', 'results_BRCA1_1010'),
        ('BRCA1 100kb', 'results_BRCA1_100100'),
        ('BRCA1 10kb+HiC', 'results_BRCA1_10kb_hic'),
        ('BRCA2 10kb', 'results_BRCA2_1010'),
        ('BRCA2 100kb', 'results_BRCA2_100100'),
        ('BRCA2 10kb+HiC', 'results_BRCA2_10kb_hic'),
        ('TNBC 10kb', 'results_TNBC_1010'),
        ('TNBC 100kb', 'results_TNBC_100100'),
        ('TNBC 10kb+HiC', 'results_TNBC_10kb_hic'),
    ]

    matrices = {}
    all_pathways = set()

    for col_name, dir_key in runs:
        match = glob.glob(os.path.join(results_dir, dir_key, 'qc_gsea/GSEA/*.tsv'))
        if match:
            df = pd.read_csv(match[0], sep='\t')
            matrices[col_name] = df.set_index('PATHWAY')
            all_pathways.update(df['PATHWAY'].tolist())

    if not matrices:
        print(f"No GSEA TSV files found in {results_dir}")
        return

    # Sort pathways by overall maximum Z-score or minimum P-value across all runs
    pathway_max_z = {}
    pathway_min_p = {}
    for p in all_pathways:
        max_z = 0.0
        min_p = 1.0
        for col_name, _ in runs:
            if col_name in matrices and p in matrices[col_name].index:
                row = matrices[col_name].loc[p]
                max_z = max(max_z, float(row['MEAN_Z']))
                min_p = min(min_p, float(row['P_VALUE']))
        pathway_max_z[p] = max_z
        pathway_min_p[p] = min_p

    sorted_pathways = sorted(all_pathways, key=lambda x: (pathway_min_p[x], -pathway_max_z[x]))

    # Build Matrix Rows
    table_rows_html = ""
    sig_summary = {col_name: 0 for col_name, _ in runs}

    for p in sorted_pathways:
        row_html = f"<tr><td class='pathway-name'><strong>{p}</strong></td>"
        
        for col_name, _ in runs:
            if col_name in matrices and p in matrices[col_name].index:
                data_row = matrices[col_name].loc[p]
                mean_z = float(data_row['MEAN_Z'])
                p_val = float(data_row['P_VALUE'])
                top_gene = str(data_row['TOP_GENE'])
                n_matched = int(data_row['N_MATCHED'])
                n_genes = int(data_row['N_GENES'])

                p_str = f"{p_val:.6f}" if p_val >= 1e-4 else f"{p_val:.2e}"
                title_attr = f"title='{p} &#10;Dataset: {col_name} &#10;Mean Z: {mean_z:.4f} &#10;P-Value: {p_str} &#10;Top Gene: {top_gene} &#10;Matched: {n_matched}/{n_genes}'"

                if p_val < 0.05:
                    sig_summary[col_name] += 1
                    cell_html = f"<td class='cell-de' {title_attr}><span class='de-badge'>Z={mean_z:.3f}<br><small>(p={p_str})</small></span></td>"
                else:
                    cell_html = f"<td class='cell-normal' {title_attr}><span class='normal-text'>Z={mean_z:.3f}<br><small>(p={p_str})</small></span></td>"
            else:
                cell_html = "<td class='cell-missing'>-</td>"

            row_html += cell_html

        row_html += "</tr>\n"
        table_rows_html += row_html

    today_str = datetime.date.today().strftime("%d %b %Y")

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GE2TS - GSEA 9-Dataset Summary Matrix for PI Review</title>
    <style>
        :root {{
            --bg-dark: #0f172a;
            --card-bg: #1e293b;
            --card-header: #334155;
            --text-light: #f8fafc;
            --text-muted: #94a3b8;
            --accent-red: #ef4444;
            --accent-red-bg: rgba(239, 68, 68, 0.25);
            --accent-green: #10b981;
            --accent-blue: #3b82f6;
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
            max-width: 1500px;
            margin: 0 auto;
        }}
        header {{
            border-bottom: 2px solid var(--border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            display: flex;
            justify-content: space-between;
            align-items: flex-start;
            flex-wrap: wrap;
            gap: 1rem;
        }}
        h1 {{
            margin: 0 0 0.5rem 0;
            font-size: 2rem;
            color: var(--text-light);
        }}
        .subtitle {{
            color: var(--text-muted);
            margin: 0;
            font-size: 1rem;
        }}
        .header-meta {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 0.5rem;
            padding: 0.75rem 1.25rem;
            font-size: 0.85rem;
            color: var(--text-muted);
        }}
        .legend-box {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 0.5rem;
            padding: 1rem 1.5rem;
            margin-bottom: 2rem;
            display: flex;
            align-items: center;
            gap: 2rem;
            flex-wrap: wrap;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
            font-size: 0.9rem;
        }}
        .de-badge {{
            background-color: var(--accent-red-bg);
            color: #fca5a5;
            border: 1px solid var(--accent-red);
            padding: 0.25rem 0.5rem;
            border-radius: 0.375rem;
            font-weight: 700;
            font-size: 0.825rem;
            display: inline-block;
            text-align: center;
            line-height: 1.2;
        }}
        .normal-text {{
            color: #64748b;
            font-size: 0.8rem;
            display: inline-block;
            text-align: center;
            line-height: 1.2;
        }}
        .table-container {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 0.75rem;
            overflow-x: auto;
            margin-bottom: 2.5rem;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            text-align: center;
            font-size: 0.875rem;
        }}
        th {{
            background-color: var(--card-header);
            color: var(--text-light);
            font-weight: 600;
            padding: 0.85rem 0.5rem;
            border-bottom: 2px solid var(--border);
            border-right: 1px solid var(--border);
            font-size: 0.8rem;
            letter-spacing: 0.02em;
        }}
        th.pathway-col {{
            text-align: left;
            padding-left: 1rem;
            min-width: 260px;
        }}
        td {{
            padding: 0.6rem 0.4rem;
            border-bottom: 1px solid var(--border);
            border-right: 1px solid var(--border);
        }}
        td.pathway-name {{
            text-align: left;
            padding-left: 1rem;
            font-size: 0.85rem;
        }}
        tr:hover {{
            background-color: rgba(51, 65, 85, 0.6);
        }}
        .filter-controls {{
            display: flex;
            gap: 1rem;
            margin-bottom: 1.25rem;
            align-items: center;
        }}
        .filter-controls input {{
            background-color: var(--card-bg);
            border: 1px solid var(--border);
            color: var(--text-light);
            padding: 0.5rem 1rem;
            border-radius: 0.375rem;
            font-size: 0.9rem;
            width: 320px;
        }}
        .btn {{
            background-color: var(--accent-blue);
            color: white;
            border: none;
            padding: 0.5rem 1rem;
            border-radius: 0.375rem;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.875rem;
        }}
        .btn:hover {{
            background-color: #2563eb;
        }}
        footer {{
            margin-top: 3rem;
            border-top: 1px solid var(--border);
            padding-top: 1.5rem;
            text-align: center;
            color: var(--text-muted);
            font-size: 0.875rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div>
                <h1>🧬 GE2TS GSEA 9-Dataset Differential Enrichment Matrix</h1>
                <p class="subtitle">50 Hallmark Pathways Evaluated Across CIMBA BRCA1, BRCA2 & TNBC Meta-Analysis (3 Window Modes)</p>
            </div>
            <div class="header-meta">
                <div><strong>Report Target:</strong> PI Functional Validation Summary</div>
                <div><strong>Total Pathways:</strong> 50 Hallmark Sets</div>
                <div><strong>Date Generated:</strong> {today_str}</div>
            </div>
        </header>

        <!-- Legend -->
        <div class="legend-box">
            <div><strong>Matrix Key:</strong></div>
            <div class="legend-item">
                <span class="de-badge">Z=0.573<br><small>(p=0.000001)</small></span>
                <span><strong>Differentially Enriched / Significant ($P &lt; 0.05$)</strong> (Highlighted in RED)</span>
            </div>
            <div class="legend-item">
                <span class="normal-text">Z=0.134<br><small>(p=0.0961)</small></span>
                <span>Not Significant (P &ge; 0.05)</span>
            </div>
        </div>

        <!-- Filter Bar -->
        <div class="filter-controls">
            <input type="text" id="searchInput" onkeyup="filterMatrix()" placeholder="Search pathway name...">
            <button class="btn" onclick="toggleOnlyDE()">Toggle Only DE Pathways (P &lt; 0.05)</button>
        </div>

        <!-- Matrix Table (50 Pathways x 9 Datasets) -->
        <div class="table-container">
            <table id="matrixTable">
                <thead>
                    <tr>
                        <th class="pathway-col">Hallmark Pathway Name (50 Sets)</th>
                        <th>BRCA1<br>10kb</th>
                        <th>BRCA1<br>100kb</th>
                        <th>BRCA1<br>10kb+HiC</th>
                        <th>BRCA2<br>10kb</th>
                        <th>BRCA2<br>100kb</th>
                        <th>BRCA2<br>10kb+HiC</th>
                        <th>TNBC<br>10kb</th>
                        <th>TNBC<br>100kb</th>
                        <th>TNBC<br>10kb+HiC</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows_html}
                </tbody>
            </table>
        </div>

        <footer>
            Generated automatically by GE2TS Pipeline for PI Validation
        </footer>
    </div>

    <script>
        let onlyDE = false;

        function filterMatrix() {{
            const input = document.getElementById('searchInput').value.toUpperCase();
            const table = document.getElementById('matrixTable');
            const tr = table.getElementsByTagName('tr');

            for (let i = 1; i < tr.length; i++) {{
                const row = tr[i];
                const textContent = row.getElementsByTagName('td')[0].textContent.toUpperCase();
                const hasDE = row.getElementsByClassName('de-badge').length > 0;

                const matchesSearch = textContent.indexOf(input) > -1;
                const matchesDE = !onlyDE || hasDE;

                if (matchesSearch && matchesDE) {{
                    row.style.display = '';
                }} else {{
                    row.style.display = 'none';
                }}
            }}
        }}

        function toggleOnlyDE() {{
            onlyDE = !onlyDE;
            filterMatrix();
        }}
    </script>
</body>
</html>
"""

    out_dir = os.path.dirname(output_html)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Successfully regenerated consolidated PI matrix report at: {output_html}")

if __name__ == '__main__':
    generate_consolidated_report()
