import os
import tempfile
import pytest
import pandas as pd
from bin import gwas_match
from bin.engine import DataEngine

def test_match_variants_pandas_fallback(monkeypatch):
    # Force pandas fallback
    monkeypatch.setattr(gwas_match, "get_engine", lambda: DataEngine(engine_type='pandas'))
    
    # Create dummy files
    with tempfile.TemporaryDirectory() as tmpdir:
        bim_file = os.path.join(tmpdir, "test.bim")
        gwas_file = os.path.join(tmpdir, "test.gwas")
        output_file = os.path.join(tmpdir, "test.matched.tsv")
        
        with open(bim_file, "w") as f:
            # CHR, ID, CM, POS, A1, A2
            f.write("1\trs1\t0\t100\tA\tT\n")
            
        with open(gwas_file, "w") as f:
            f.write("SNP\tCHR\tBP\tA1\tA2\n")
            f.write("rs1\t1\t100\tA\tT\n")
            
        gwas_match.match_variants(gwas_file, output_file, bim_file)
        
        assert os.path.exists(output_file)
        df = pd.read_csv(output_file, sep='\t')
        assert df["SNP"].iloc[0] == "rs1"

def test_match_variants_polars_if_available(monkeypatch):
    # Only run if polars is actually installed
    from bin.engine import HAS_POLARS
    if not HAS_POLARS:
        pytest.skip("Polars not available")
        
    monkeypatch.setattr(gwas_match, "get_engine", lambda: DataEngine(engine_type='polars'))
        
    with tempfile.TemporaryDirectory() as tmpdir:
        bim_file = os.path.join(tmpdir, "test.bim")
        gwas_file = os.path.join(tmpdir, "test.gwas")
        output_file = os.path.join(tmpdir, "test.matched.tsv")
        
        with open(bim_file, "w") as f:
            f.write("1\trs1\t0\t100\tA\tT\n")
            
        with open(gwas_file, "w") as f:
            f.write("SNP\tCHR\tBP\tA1\tA2\n")
            f.write("rs1\t1\t100\tA\tT\n")
            
        gwas_match.match_variants(gwas_file, output_file, bim_file)
        
        assert os.path.exists(output_file)
        df = pd.read_csv(output_file, sep='\t')
        assert df["SNP"].iloc[0] == "rs1"
