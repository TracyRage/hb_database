from functools import partial, reduce
import pandas as pd
import pathlib


class CreateAnnotationTable:
    """Generate annotation table (Cdd, Ko, Pfam)"""
    def __init__(self, path_to_files):
        self.cdd_files = pathlib.Path(path_to_files).glob('**/*.cdd')
        self.ko_files = pathlib.Path(path_to_files).glob('**/*.ko')
        self.pfam_files = pathlib.Path(path_to_files).glob('**/*.pfam')

    @property
    def create_file_list(self):
        """Create file lists (Cdd, Ko, Pfam) to be used"""
        cdd_files = sorted([str(file) for file in self.cdd_files])
        ko_files = sorted([str(file) for file in self.ko_files])
        pfam_files = sorted([str(file) for file in self.pfam_files])
        return cdd_files, ko_files, pfam_files

    def load_cdd_data(self, record):
        """Create Cdd dataframe"""
        cdd_df = pd.read_csv(record, sep='\t', names=['UniProt', 'Cdd'])
        return cdd_df

    def load_ko_data(self, record):
        """Create Ko dataframe"""

        ko_df = pd.read_csv(record, sep='\t', names=['UniProt', 'Ko'])
        return ko_df

    def load_pfam_data(self, record):
        """Create Pfam dataframe"""
        pfam_raw = pd.read_csv(record, sep='\s+', comment='#', header=None)
        pfam_cols = pfam_raw.iloc[:, 0:6]
        pfam_cut = pfam_cols.drop(pfam_cols.iloc[:, 1:5], axis=1)
        pfam_final = pfam_cut.rename(columns={0: 'UniProt', 5: 'Pfam'})
        return pfam_final

    def chain_df(self, cdd_files, ko_files, pfam_files):
        """Make list of dataframes"""
        df_list = [[
            self.load_cdd_data(x),
            self.load_ko_data(y),
            self.load_pfam_data(z)
        ] for x, y, z in zip(sorted(cdd_files), sorted(ko_files),
                             sorted(pfam_files))]
        return df_list

    def create_table(self, df_list):
        """Method for annotation table creation"""
        dfs = [self.merge_df(x) for x in df_list]
        return dfs

    def merge_df(self, dfs):
        """Helper #1 for merging dfs"""
        merge = partial(pd.merge, on=['UniProt'], how='outer')
        final_df = reduce(merge, dfs)
        return final_df.drop_duplicates(subset='UniProt')

    def write_dfs(self, paths, dfs):
        """Helper #2 for writing to tsv files"""
        [self.to_tsv(x, y) for x, y in zip(dfs, paths)]

    def to_tsv(self, df, path):
        """Helper #3 for writing to tsv files"""
        write_out = df.to_csv(pathlib.Path(path).with_suffix('.tsv'),
                              sep='\t',
                              na_rep='-',
                              index=False,
                              encoding='utf-8')
