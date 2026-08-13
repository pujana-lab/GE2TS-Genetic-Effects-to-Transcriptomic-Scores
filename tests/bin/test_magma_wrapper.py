import os
import tempfile
import pandas as pd
from bin.magma_wrapper import run_magma

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
