import os
import sys
import tempfile
import pytest
import pandas as pd
from bin.qc_gsea import run_gsea_and_qc, main, parse_gmt

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

        pathway_df = pd.read_csv(os.path.join(out_dir, "sample_gsea_pathways.tsv"), sep='\t')
        assert "PATHWAY" in pathway_df.columns
        assert len(pathway_df) == 2

def test_parse_gmt_nonexistent():
    pathways = parse_gmt("non_existent.gmt")
    assert pathways == {}

def test_qc_gsea_symbol_fallback():
    # Genes file without SYMBOL
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
        
        # Verify
        pathway_df = pd.read_csv(os.path.join(out_dir, "no_symbol_gsea_pathways.tsv"), sep='\t')
        # Check if it matched GENE_A
        # The test pathway file should contain GENE_A
        # Let's assume dummy_pathways.gmt contains GENE_A
        assert pathway_df["N_MATCHED"].iloc[0] > 0

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
