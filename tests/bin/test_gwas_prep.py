import os
import json
import tempfile
import pytest
from bin.gwas_prep import process_gwas

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
        # rs1003 has A2='AG' (length 2), so it should be filtered out when filter_indels=True
        assert summary["FilteredIndels"] == 1
        assert summary["TotalVariantsPostQC"] == 4
