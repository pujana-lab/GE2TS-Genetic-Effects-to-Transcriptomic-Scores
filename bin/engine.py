import pandas as pd
try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

class DataEngine:
    def __init__(self, engine_type=None):
        if engine_type == 'polars' and HAS_POLARS:
            self.engine = 'polars'
        else:
            self.engine = 'pandas'

    def read_csv(self, file_path, separator='\t', **kwargs):
        if self.engine == 'polars':
            return pl.read_csv(file_path, separator=separator, **kwargs)
        else:
            return pd.read_csv(file_path, sep=separator, **kwargs)

    def read_genes_out(self, genes_out_path):
        if self.engine == 'polars':
            try:
                return pl.read_csv(genes_out_path, separator='\t').to_pandas()
            except Exception:
                pass
        try:
            return pd.read_csv(genes_out_path, sep=r'\s+', engine='python')
        except Exception:
            return pd.read_csv(genes_out_path, sep='\t')

    def add_symbol_column(self, genes_out, gene_loc):
        if self.engine == 'polars':
            genes_df = pl.read_csv(genes_out, separator='\t')
            loc_df = pl.read_csv(gene_loc, separator='\t', has_header=False)
            mapping_df = loc_df.select([pl.col(loc_df.columns[0]).alias("GENE"), pl.col(loc_df.columns[-1]).alias("SYMBOL")])
            genes_df = genes_df.join(mapping_df, on="GENE", how="left")
            genes_df.write_csv(genes_out, separator='\t')
        else:
            genes_df = pd.read_csv(genes_out, sep=r'\s+')
            loc_df = pd.read_csv(gene_loc, sep=r'\s+', header=None)
            mapping = dict(zip(loc_df[0], loc_df[loc_df.columns[-1]]))
            genes_df['SYMBOL'] = genes_df['GENE'].map(mapping)
            genes_df.to_csv(genes_out, sep='\t', index=False)

    def load_bim_and_create_maps(self, bim_file):
        if self.engine == 'polars':
            ref = pl.read_csv(
                bim_file,
                separator='\t',
                has_header=False,
                new_columns=['CHR', 'ID', 'CM', 'POS', 'A1', 'A2'],
                schema_overrides={'CHR': pl.String, 'POS': pl.Int64, 'ID': pl.String, 'A1': pl.String, 'A2': pl.String}
            )
            comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            ref = ref.with_columns([
                pl.concat_str([pl.col("CHR"), pl.col("POS"), pl.col("A1"), pl.col("A2")], separator=":").alias("key_exact"),
                pl.concat_str([pl.col("CHR"), pl.col("POS"), pl.col("A2"), pl.col("A1")], separator=":").alias("key_swapped"),
                pl.col("A1").map_elements(lambda x: "".join([comp_map.get(b, b) for b in x]), return_dtype=pl.String).alias("A1_c"),
                pl.col("A2").map_elements(lambda x: "".join([comp_map.get(b, b) for b in x]), return_dtype=pl.String).alias("A2_c")
            ]).with_columns(
                pl.concat_str([pl.col("CHR"), pl.col("POS"), pl.col("A1_c"), pl.col("A2_c")], separator=":").alias("key_strand")
            )
            return dict(zip(ref['ID'], ref['ID'])), dict(zip(ref['key_exact'], ref['ID'])), dict(zip(ref['key_swapped'], ref['ID'])), dict(zip(ref['key_strand'], ref['ID']))
        else:
            ref = pd.read_csv(bim_file, sep='\t', header=None, names=['CHR', 'ID', 'CM', 'POS', 'A1', 'A2'], 
                              dtype={'CHR': str, 'POS': int, 'ID': str, 'A1': str, 'A2': str})
            comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            def complement(seq):
                return "".join([comp_map.get(b, b) for b in seq])
            ref['key_exact'] = ref['CHR'].astype(str) + ":" + ref['POS'].astype(str) + ":" + ref['A1'] + ":" + ref['A2']
            ref['key_swapped'] = ref['CHR'].astype(str) + ":" + ref['POS'].astype(str) + ":" + ref['A2'] + ":" + ref['A1']
            ref['A1_c'] = ref['A1'].apply(complement)
            ref['A2_c'] = ref['A2'].apply(complement)
            ref['key_strand'] = ref['CHR'].astype(str) + ":" + ref['POS'].astype(str) + ":" + ref['A1_c'] + ":" + ref['A2_c']
            return dict(zip(ref['ID'], ref['ID'])), dict(zip(ref['key_exact'], ref['ID'])), dict(zip(ref['key_swapped'], ref['ID'])), dict(zip(ref['key_strand'], ref['ID']))

def get_engine():
    return DataEngine(engine_type='polars' if HAS_POLARS else 'pandas')
