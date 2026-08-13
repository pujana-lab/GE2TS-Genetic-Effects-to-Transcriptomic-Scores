import os
import sys
import tempfile
import subprocess
import pandas as pd
from bin.magma_wrapper import run_magma, main

def test_run_magma_mock():
    gwas_path = "assets/test_data/sample_gwas.tsv"

    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "test_run")
        
        run_magma(
            gwas_path=gwas_path,
            bfile_prefix="dummy_ref",
            gene_loc="dummy.gene.loc",
            window="10,10",
            out_prefix=out_prefix,
            sample_size="10000",
            mock=True
        )

        genes_out = f"{out_prefix}.genes.out"
        genes_raw = f"{out_prefix}.genes.raw"

        assert os.path.exists(genes_out)
        assert os.path.exists(genes_raw)

        df = pd.read_csv(genes_out, sep='\t')
        assert "SYMBOL" in df.columns
        assert "P" in df.columns
        assert len(df) == 3

def test_run_magma_mock_with_hic():
    gwas_path = "assets/test_data/sample_gwas.tsv"
    hic_bed = "assets/test_data/dummy_ref.bed"

    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "test_hic")
        
        run_magma(
            gwas_path=gwas_path,
            bfile_prefix="dummy_ref",
            gene_loc="dummy.gene.loc",
            window="10,10",
            out_prefix=out_prefix,
            sample_size="10000",
            hic_bed=hic_bed,
            mock=True
        )

        genes_out = f"{out_prefix}.genes.out"
        assert os.path.exists(genes_out)
        df = pd.read_csv(genes_out, sep='\t')
        assert len(df) == 3

def test_run_magma_mock_nonexistent_gwas():
    gwas_path = "assets/test_data/non_existent.tsv"

    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "test_nonexistent")
        
        run_magma(
            gwas_path=gwas_path,
            bfile_prefix="dummy_ref",
            gene_loc="dummy.gene.loc",
            window="10,10",
            out_prefix=out_prefix,
            sample_size="10000",
            mock=True
        )

        assert os.path.exists(f"{out_prefix}.genes.out")

def test_run_magma_mock_exception_gwas():
    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "test_exception")
        # Pass tmpdir itself which is a directory, causing pd.read_csv to raise IsADirectoryError / Exception
        run_magma(
            gwas_path=tmpdir,
            bfile_prefix="dummy_ref",
            gene_loc="dummy.gene.loc",
            window="10,10",
            out_prefix=out_prefix,
            sample_size="10000",
            mock=True
        )
        assert os.path.exists(f"{out_prefix}.genes.out")

def test_run_magma_binary_not_found(monkeypatch):
    gwas_path = "assets/test_data/sample_gwas.tsv"

    def dummy_run_error(cmd, capture_output=False, text=False, check=False):
        raise FileNotFoundError("magma not found")

    monkeypatch.setattr(subprocess, "run", dummy_run_error)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "test_nofile")
        run_magma(
            gwas_path=gwas_path,
            bfile_prefix="dummy_ref",
            gene_loc="dummy.gene.loc",
            window="10,10",
            out_prefix=out_prefix,
            sample_size="10000",
            mock=False
        )
        assert os.path.exists(f"{out_prefix}.genes.out")

def test_run_magma_real_binary_simulated(monkeypatch):
    gwas_path = "assets/test_data/sample_gwas.tsv"

    class DummyCompletedProcess:
        returncode = 0
        stdout = "MAGMA v1.10"
        stderr = ""

    def dummy_run(cmd, capture_output=False, text=False, check=False):
        return DummyCompletedProcess()

    monkeypatch.setattr(subprocess, "run", dummy_run)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "test_real_sim")
        
        run_magma(
            gwas_path=gwas_path,
            bfile_prefix="dummy_ref",
            gene_loc="dummy.gene.loc",
            window="10,10",
            out_prefix=out_prefix,
            sample_size="10000",
            mock=False
        )

def test_magma_wrapper_cli(monkeypatch):
    gwas_path = "assets/test_data/sample_gwas.tsv"
    hic_bed = "assets/test_data/dummy_ref.bed"

    with tempfile.TemporaryDirectory() as tmpdir:
        out_prefix = os.path.join(tmpdir, "cli_magma")
        
        test_args = [
            "magma_wrapper.py",
            "--gwas", gwas_path,
            "--hic-bed", hic_bed,
            "--out-prefix", out_prefix,
            "--mock"
        ]
        monkeypatch.setattr(sys, "argv", test_args)
        main()

        assert os.path.exists(f"{out_prefix}.genes.out")
        assert os.path.exists(f"{out_prefix}.genes.raw")
