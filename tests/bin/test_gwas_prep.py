import os
import json
import sys
import tempfile
import pytest
from bin.gwas_prep import process_gwas, main

def test_process_gwas():
    input_path = "assets/test_data/sample_gwas.tsv"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_tsv = os.path.join(tmpdir, "out_gwas.tsv")
        qc_summary = os.path.join(tmpdir, "qc_summary.json")

        process_gwas(
            input_path=input_path,
            output_path=output_tsv,
            qc_summary_path=qc_summary,
            filter_indels=True
        )

        assert os.path.exists(output_tsv)
        assert os.path.exists(qc_summary)

        with open(qc_summary, "r") as f:
            summary = json.load(f)

        assert "TotalVariantsPreQC" in summary
        assert summary["TotalVariantsPreQC"] == 5
        assert summary["FilteredIndels"] == 1
        assert summary["TotalVariantsPostQC"] == 4

def test_process_gwas_csv_synonyms():
    csv_content = 'MarkerName,CHR,position,Allele1,Allele2,Effect,StdErr,P-value\nrs2001,1,12345,A,G,0.12,0.05,0.003\n'
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "input.csv")
        with open(csv_path, "w") as f:
            f.write(csv_content)

        output_tsv = os.path.join(tmpdir, "out_synonyms.tsv")
        qc_summary = os.path.join(tmpdir, "qc_synonyms.json")

        process_gwas(
            input_path=csv_path,
            output_path=output_tsv,
            qc_summary_path=qc_summary,
            filter_indels=True
        )

        assert os.path.exists(output_tsv)
        with open(qc_summary, "r") as f:
            summary = json.load(f)
        assert summary["TotalVariantsPostQC"] == 1

def test_process_gwas_gz_output():
    input_path = "assets/test_data/sample_gwas.tsv"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_gz = os.path.join(tmpdir, "out_gwas.tsv.gz")
        qc_summary = os.path.join(tmpdir, "qc_summary.json")

        process_gwas(
            input_path=input_path,
            output_path=output_gz,
            qc_summary_path=qc_summary,
            filter_indels=False
        )

        assert os.path.exists(output_gz)

def test_process_gwas_pandas_fallback(monkeypatch):
    input_path = "assets/test_data/sample_gwas.tsv"
    monkeypatch.setitem(sys.modules, "polars", None)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_tsv = os.path.join(tmpdir, "out_pandas.tsv")
        output_gz = os.path.join(tmpdir, "out_pandas.tsv.gz")
        qc_summary = os.path.join(tmpdir, "qc_summary_pandas.json")

        process_gwas(
            input_path=input_path,
            output_path=output_tsv,
            qc_summary_path=qc_summary,
            filter_indels=True
        )
        assert os.path.exists(output_tsv)

        process_gwas(
            input_path=input_path,
            output_path=output_gz,
            qc_summary_path=qc_summary,
            filter_indels=True
        )
        assert os.path.exists(output_gz)

def test_gwas_prep_cli(monkeypatch):
    input_path = "assets/test_data/sample_gwas.tsv"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_tsv = os.path.join(tmpdir, "cli_gwas.tsv")
        qc_summary = os.path.join(tmpdir, "cli_qc.json")

        test_args = [
            "gwas_prep.py",
            "--input", input_path,
            "--output-tsv", output_tsv,
            "--qc-summary", qc_summary
        ]
        monkeypatch.setattr(sys, "argv", test_args)
        main()

        assert os.path.exists(output_tsv)
        assert os.path.exists(qc_summary)
