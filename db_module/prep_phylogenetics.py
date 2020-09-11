import pandas as pd
from Bio import SeqIO
import pathlib
import re


class PhyloPrep:
    """Prepare fasta files for phylogenetic assessment"""
    def __init__(self, initial_table, ncbi_fa, blast_fa, query_fa):
        self.initial_table = pd.read_csv(initial_table)
        self.ncbi_fa = [str(x) for x in pathlib.Path(ncbi_fa).glob("**/*.fa")]
        self.blast_fa = [
            str(y) for y in pathlib.Path(blast_fa).glob("**/*.fa")
        ]
        self.query_fa = [
            str(z) for z in pathlib.Path(query_fa).glob("**/*.filtered")
        ]
        self.query_fa_raw = [
            str(w) for w in pathlib.Path(query_fa).glob('**/*.fa')
        ]

    @property
    def sort_files(self):
        """Sort fasta files"""
        self.ncbi_fa.sort(key=lambda x: x[re.search(r'(fasta?<=).*', x):])
        self.blast_fa.sort(key=lambda x: x[re.search(r'(hits?<=).*', x):])
        self.query_fa.sort(key=lambda x: x[re.search(r'(hits?<=).*', x):])
        self.query_fa_raw.sort(key=lambda x: x[re.search(r'(hits?<=).*', x):])
        return self.ncbi_fa, self.blast_fa, self.query_fa, self.query_fa_raw

    @property
    def get_index(self):
        """Get indexes from initial table"""
        indexes_1 = list(self.initial_table['accession'].sort_values())
        indexes_2 = list(self.initial_table['enzyme'].sort_values())
        return indexes_1, indexes_2

    def get_hits_number(self, blast_fa, query_fa, unfiltered_fa):
        """Count number of hits in each file and make a list"""
        blast_fa_count = [
            self.get_hits_number_helper(file) for file in blast_fa
        ]
        query_fa_count = [
            self.get_hits_number_helper(file) for file in query_fa
        ]
        query_fa_count_raw = [
            self.get_hits_number_helper(file) for file in unfiltered_fa
        ]
        return blast_fa_count, query_fa_count, query_fa_count_raw

    def get_hits_number_helper(self, fasta_file):
        """Parse fasta file and find out the number of sequences"""
        seq_length = [
            sequence.id for sequence in SeqIO.parse(fasta_file, 'fasta')
        ]
        return len(seq_length)

    def generate_df(self, index_1, ncbi_files, sp_merged, blast_fa_count,
                    index_2, sp_query, query_fa_count, query_fa_count_raw):
        """Generate fasta files dataframes"""
        sp_final = list(zip(index_1, ncbi_files, sp_merged, blast_fa_count))
        sp_df = pd.DataFrame(
            sp_final, columns=['accession', 'ncbi_ref', 'sp_ref', 'sp_ref#'])
        query_final = list(
            zip(index_2, sp_query, query_fa_count, query_fa_count_raw))
        query_df = pd.DataFrame(
            query_final,
            columns=['enzyme', 'sp_query', 'filtered#', 'unfiltered#'])
        return sp_df, query_df

    def create_final_table(self, sp_df, query_df):
        """Merge fasta files dataframes to initial table and write to a file"""
        first_merge = pd.merge(left=self.initial_table,
                               right=query_df,
                               on='enzyme',
                               how='outer')
        second_merge = pd.merge(left=self.initial_table,
                                right=sp_df,
                                on='accession',
                                how='outer')
        final_merge = pd.merge(left=first_merge,
                               right=second_merge,
                               on=['accession', 'enzyme', 'alpha_subunits'],
                               how='outer')
        final_merge.to_csv('analyses/annotated_table.tsv',
                           sep='\t',
                           index=False,
                           encoding='utf-8')
