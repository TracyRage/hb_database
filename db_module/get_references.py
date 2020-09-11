from Bio import SeqIO, Entrez
import pandas as pd
import re


class GetReference:
    """Get references for manually chosen
    NCBI sequences"""
    def __init__(self, email, table_input):
        self.email = email
        self.table_input = pd.read_csv(table_input, sep='\t')

    @property
    def parse_table(self):
        """Parse input csv table"""
        indexes = list(
            self.table_input.sort_values(by='accession')['accession'])
        return indexes

    def get_ref(self, accession_list):
        """Get references for each accession number"""
        references = [
            self.get_ref_helper_2(self.get_ref_helper(accession))
            for accession in accession_list
        ]
        return references

    def get_ref_helper(self, accession):
        """Helper #1 for get_ref method"""
        Entrez.email = self.email
        with Entrez.efetch(
                db='protein',
                rettype="gb",
                retmode="text",
                id=accession,
        ) as handle:
            seq_record = SeqIO.read(handle, "gb")
            return seq_record.annotations['references']

    def get_ref_helper_2(self, reference_list):
        """Helper #2 for get_ref method"""
        references = [
            self.match_line(str(reference)) for reference in reference_list
        ]
        return references

    def match_line(self, line):
        """Helper for get_title method"""
        pattern_1 = re.compile(r'title:\s+(.*)')
        if pattern_1.search(line):
            title_2 = pattern_1.findall(line)[0]
            return title_2

    def get_primary_ref(self, reference_list):
        """Get primary reference"""
        primary = [reference[0] for reference in reference_list]
        return primary

    def generate_dataframe(self, index, references):
        """Generate reference dataframe"""
        df = list(zip(index, references))
        ref_df = pd.DataFrame(df, columns=['accession', 'references'])
        return ref_df

    def merge_with_input_table(self, reference_df):
        """Merge reference dataframe with input table"""
        merge_df = pd.merge(left=self.table_input,
                            right=reference_df,
                            on='accession',
                            how='outer')
        merge_df.to_csv('analyses/annotated_table_with_ref.tsv',
                        sep='\t',
                        index=False,
                        encoding='utf-8')
