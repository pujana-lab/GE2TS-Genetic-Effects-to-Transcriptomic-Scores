import os
import sys
import tempfile
import pytest
import pandas as pd
from bin import gwas_match
from bin.engine import DataEngine

def test_match_variant_row_logic():
    id_map = {"rs100": "rs100"}
    exact_map = {"1:100:A:G": "rs101"}
    swap_map = {"1:200:G:A": "rs102"}
    strand_map = {"1:300:T:C": "rs103"}

    # 1. ID Match
    row1 = {"SNP": "rs100", "CHR": "1", "BP": "50", "A1": "A", "A2": "G"}
    res_row1, matched1 = gwas_match.match_variant_row(row1, id_map, exact_map, swap_map, strand_map)
    assert matched1 is True
    assert res_row1["SNP"] == "rs100"

    # 2. Exact Key Match
    row2 = {"SNP": "temp_2", "CHR": "1", "BP": "100", "A1": "A", "A2": "G"}
    res_row2, matched2 = gwas_match.match_variant_row(row2, id_map, exact_map, swap_map, strand_map)
    assert matched2 is True
    assert res_row2["SNP"] == "rs101"

    # 3. Swap Key Match
    row3 = {"SNP": "temp_3", "CHR": "1", "BP": "200", "A1": "G", "A2": "A"}
    res_row3, matched3 = gwas_match.match_variant_row(row3, id_map, exact_map, swap_map, strand_map)
    assert matched3 is True
    assert res_row3["SNP"] == "rs102"

    # 4. Strand Key Match
    row4 = {"SNP": "temp_4", "CHR": "1", "BP": "300", "A1": "T", "A2": "C"}
    res_row4, matched4 = gwas_match.match_variant_row(row4, id_map, exact_map, swap_map, strand_map)
    assert matched4 is True
    assert res_row4["SNP"] == "rs103"

    # 5. No Match
    row5 = {"SNP": "temp_5", "CHR": "1", "BP": "999", "A1": "C", "A2": "T"}
    res_row5, matched5 = gwas_match.match_variant_row(row5, id_map, exact_map, swap_map, strand_map)
    assert matched5 is False
    assert res_row5["SNP"] == "temp_5"

def test_stream_and_match_gwas_empty_and_gz():
    with tempfile.TemporaryDirectory() as tmpdir:
        empty_gwas = os.path.join(tmpdir, "empty.gwas")
        out_empty = os.path.join(tmpdir, "out_empty.tsv")
        with open(empty_gwas, "w") as f:
            f.write("") # Empty file
            
        matched, total = gwas_match.stream_and_match_gwas(empty_gwas, out_empty, ({}, {}, {}, {}))
        assert matched == 0
        assert total == 0

        gz_gwas = os.path.join(tmpdir, "test.gwas.gz")
        out_gz = os.path.join(tmpdir, "out_gz.tsv")
        import gzip
        with gzip.open(gz_gwas, "wt") as f:
            f.write("SNP\tCHR\tBP\tA1\tA2\n")
            f.write("rs1\t1\t100\tA\tT\n")
            
        id_map = {"rs1": "rs1"}
        matched_gz, total_gz = gwas_match.stream_and_match_gwas(gz_gwas, out_gz, (id_map, {}, {}, {}))
        assert matched_gz == 1
        assert total_gz == 1

def test_match_variants_pandas_fallback(monkeypatch):
    monkeypatch.setattr(gwas_match, "get_engine", lambda: DataEngine(engine_type='pandas'))
    
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

def test_match_variants_polars_if_available(monkeypatch):
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

def test_gwas_match_cli(monkeypatch):
    with tempfile.TemporaryDirectory() as tmpdir:
        bim_file = os.path.join(tmpdir, "cli.bim")
        gwas_file = os.path.join(tmpdir, "cli.gwas")
        output_file = os.path.join(tmpdir, "cli.matched.tsv")

        with open(bim_file, "w") as f:
            f.write("1\trs1\t0\t100\tA\tT\n")

        with open(gwas_file, "w") as f:
            f.write("SNP\tCHR\tBP\tA1\tA2\n")
            f.write("rs1\t1\t100\tA\tT\n")

        test_args = [
            "gwas_match.py",
            "--gwas", gwas_file,
            "--bim", bim_file,
            "--output-tsv", output_file
        ]
        monkeypatch.setattr(sys, "argv", test_args)
        gwas_match.main()

        assert os.path.exists(output_file)
