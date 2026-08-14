import os
import sys
import tempfile
import pytest
import pandas as pd
from bin import qc_gsea
from bin.qc_gsea import (
    run_gsea_and_qc, main, parse_gmt, load_genes_data,
    calculate_pathway_zscores, write_qc_summary, generate_gsea_html_report
)

def test_run_gsea_and_qc():
    dummy_genes_data = "GENE\tSYMBOL\tCHR\tSTART\tSTOP\tNSNPS\tNPARAM\tZSTAT\tP\n1\tGENE_A\t1\t5000\t35000\t10\t10\t3.45\t0.00028\n2\tGENE_B\t2\t10000\t30000\t8\t8\t2.15\t0.01578\n"
    gmt_path = "assets/test_data/dummy_pathways.gmt"

    with tempfile.TemporaryDirectory() as tmpdir:
        genes_out = os.path.join(tmpdir, "sample.genes.out")
        with open(genes_out, "w") as f:
            f.write(dummy_genes_data)

        out_dir = os.path.join(tmpdir, "GSEA")
        qc_tsv = os.path.join(tmpdir, "qc_summary.tsv")

        run_gsea_and_qc(
            genes_out_path=genes_out,
            gmt_path=gmt_path,
            out_dir=out_dir,
            qc_summary_path=qc_tsv
        )

        assert os.path.exists(qc_tsv)
        assert os.path.exists(os.path.join(out_dir, "sample_gsea_pathways.tsv"))
        assert os.path.exists(os.path.join(out_dir, "sample_gsea_report.html"))

        with open(os.path.join(out_dir, "sample_gsea_report.html"), "r") as f:
            html_text = f.read()
            assert "GE2TS GSEA Enrichment Report" in html_text
            assert "sample" in html_text

        pathway_df = pd.read_csv(os.path.join(out_dir, "sample_gsea_pathways.tsv"), sep='\t')
        assert "PATHWAY" in pathway_df.columns
        assert len(pathway_df) == 2

def test_parse_gmt_nonexistent():
    pathways = parse_gmt("non_existent.gmt")
    assert pathways == {}

def test_qc_gsea_symbol_fallback():
    dummy_genes_data = "GENE\tCHR\tSTART\tSTOP\tNSNPS\tNPARAM\tZSTAT\tP\n1\t1\t5000\t35000\t10\t10\t3.45\t0.00028\n"
    gmt_path = "assets/test_data/dummy_pathways.gmt"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        genes_out = os.path.join(tmpdir, "no_symbol.genes.out")
        with open(genes_out, "w") as f:
            f.write(dummy_genes_data)
        
        gene_loc = os.path.join(tmpdir, "mapping.gene.loc")
        with open(gene_loc, "w") as f:
            f.write("1\tCHR\tSTART\tSTOP\tGENE_A\n")
            
        out_dir = os.path.join(tmpdir, "GSEA")
        qc_tsv = os.path.join(tmpdir, "qc.tsv")
        
        run_gsea_and_qc(
            genes_out_path=genes_out,
            gmt_path=gmt_path,
            out_dir=out_dir,
            qc_summary_path=qc_tsv,
            gene_loc=gene_loc
        )
        
        pathway_df = pd.read_csv(os.path.join(out_dir, "no_symbol_gsea_pathways.tsv"), sep='\t')
        assert pathway_df["N_MATCHED"].iloc[0] > 0
        assert os.path.exists(os.path.join(out_dir, "no_symbol_gsea_report.html"))

def test_generate_gsea_html_report_empty():
    genes_df = pd.DataFrame()
    pathway_df = pd.DataFrame()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        html_out = os.path.join(tmpdir, "empty_report.html")
        generate_gsea_html_report("empty_pheno", genes_df, pathway_df, html_out)
        
        assert os.path.exists(html_out)
        with open(html_out, "r") as f:
            html = f.read()
            assert "No pathway enrichment results available" in html

def test_load_html_template_custom_and_fallback(monkeypatch):
    from bin.qc_gsea import load_html_template
    # Custom template file
    with tempfile.TemporaryDirectory() as tmpdir:
        custom_tpl = os.path.join(tmpdir, "custom.html")
        with open(custom_tpl, "w") as f:
            f.write("<html><body>{phenotype}</body></html>")
            
        tpl = load_html_template(custom_tpl)
        assert "{phenotype}" in tpl

    # Fallback template
    monkeypatch.setattr(os.path, "exists", lambda path: False)
    tpl_fallback = load_html_template("non_existent.html")
    assert "<html>" in tpl_fallback

def test_load_genes_data_nonexistent():
    with pytest.raises(FileNotFoundError):
        load_genes_data("non_existent.genes.out")

def test_load_genes_data_pandas_fallback(monkeypatch):
    from bin.engine import DataEngine
    monkeypatch.setattr(qc_gsea, "get_engine", lambda: DataEngine(engine_type='pandas'))
    dummy_genes_data = "GENE\tCHR\tSTART\tSTOP\tNSNPS\tNPARAM\tZSTAT\tP\n1\t1\t5000\t35000\t10\t10\t3.45\t0.00028\n"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        genes_out = os.path.join(tmpdir, "sample.genes.out")
        with open(genes_out, "w") as f:
            f.write(dummy_genes_data)
            
        df = load_genes_data(genes_out)
        assert "GENE" in df.columns

def test_calculate_pathway_zscores_no_match():
    genes_df = pd.DataFrame({"SYMBOL": ["GENE_X"], "ZSTAT": [1.0]})
    pathways = {"TEST_PATHWAY": {"GENE_A", "GENE_B"}}
    
    res = calculate_pathway_zscores(genes_df, pathways)
    assert len(res) == 1
    assert res.iloc[0]["N_MATCHED"] == 0
    assert res.iloc[0]["P_VALUE"] == 1.0

def test_write_qc_summary():
    genes_df = pd.DataFrame({"P": [0.01, 0.04, 0.20]})
    pathway_df = pd.DataFrame([{"PATHWAY": "PATH_1", "P_VALUE": 0.001}])
    
    with tempfile.TemporaryDirectory() as tmpdir:
        qc_path = os.path.join(tmpdir, "summary.tsv")
        write_qc_summary(genes_df, pathway_df, "pheno1", qc_path)
        
        assert os.path.exists(qc_path)
        df = pd.read_csv(qc_path, sep='\t')
        assert df["SignificantGenes_P005"].iloc[0] == 2
        assert df["TopEnrichedPathway"].iloc[0] == "PATH_1"

def test_qc_gsea_cli(monkeypatch):
    dummy_genes_data = "GENE\tSYMBOL\tCHR\tSTART\tSTOP\tNSNPS\tNPARAM\tZSTAT\tP\n1\tGENE_A\t1\t5000\t35000\t10\t10\t3.45\t0.00028\n"
    gmt_path = "assets/test_data/dummy_pathways.gmt"

    with tempfile.TemporaryDirectory() as tmpdir:
        genes_out = os.path.join(tmpdir, "cli_sample.genes.out")
        with open(genes_out, "w") as f:
            f.write(dummy_genes_data)

        out_dir = os.path.join(tmpdir, "GSEA_CLI")
        qc_tsv = os.path.join(tmpdir, "cli_qc_summary.tsv")

        test_args = [
            "qc_gsea.py",
            "--genes-out", genes_out,
            "--pathway-gmt", gmt_path,
            "--out-dir", out_dir,
            "--qc-summary-tsv", qc_tsv
        ]
        monkeypatch.setattr(sys, "argv", test_args)
        main()

        assert os.path.exists(qc_tsv)
        assert os.path.exists(os.path.join(out_dir, "cli_sample_gsea_pathways.tsv"))
        assert os.path.exists(os.path.join(out_dir, "cli_sample_gsea_report.html"))
