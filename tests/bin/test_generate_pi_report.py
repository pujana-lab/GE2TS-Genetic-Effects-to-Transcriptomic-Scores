import os
import tempfile
import pandas as pd
from bin.generate_pi_report import generate_consolidated_report

def test_generate_consolidated_report_empty():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_html = os.path.join(tmpdir, "report.html")
        generate_consolidated_report(results_dir=tmpdir, output_html=out_html)
        assert not os.path.exists(out_html)

def test_generate_consolidated_report_with_data():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a mock run folder structure
        run_dir = os.path.join(tmpdir, "results_BRCA1_1010", "qc_gsea", "GSEA")
        os.makedirs(run_dir, exist_ok=True)
        
        tsv_path = os.path.join(run_dir, "sample_gsea_pathways.tsv")
        df = pd.DataFrame([{
            "PATHWAY": "HALLMARK_DNA_REPAIR",
            "N_GENES": 100,
            "N_MATCHED": 90,
            "MEAN_Z": 0.5,
            "MEDIAN_Z": 0.4,
            "MAX_Z": 2.5,
            "TOP_GENE": "RAD52",
            "P_VALUE": 0.01
        }])
        df.to_csv(tsv_path, sep='\t', index=False)
        
        out_html = os.path.join(tmpdir, "report.html")
        generate_consolidated_report(results_dir=tmpdir, output_html=out_html)
        
        assert os.path.exists(out_html)
        with open(out_html, "r") as f:
            html = f.read()
            assert "HALLMARK_DNA_REPAIR" in html
            assert "BRCA1" in html
