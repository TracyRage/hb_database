from Bio import SwissProt
import requests
import os
import re
import pathlib
import pandas as pd
import subprocess


class GetSwissProt:
    """Getting fasta files by keywords search"""
    def __init__(self, keyword_file):
        self.keyword_file = keyword_file

    @property
    def check_output_dir(self):
        """Check if output directory exists"""
        if os.path.exists('analyses/sp_query'):
            print('Directory already exists!')
        else:
            print('Creating output directory. ')
            os.makedirs('analyses/sp_query')

    @property
    def init_table(self):
        """Extract alpha subunits and gene names"""
        input_table = pd.read_csv(self.keyword_file)
        return input_table

    def get_dfs(self, input_table):
        """Get dataframes from table"""
        df_enzyme = input_table['enzyme']
        df_subunit = input_table['alpha_subunits']
        return df_enzyme, df_subunit

    def sp_query(self, df_subunit, df_enzyme):
        """Make a SwissProt API query"""
        url = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&reviewed=False&exact_gene="
        fasta_files = [
            self.write_fa(name, self.do_request(url + sub))
            for sub, name in zip(df_subunit, df_enzyme)
        ]
        return fasta_files

    def do_request(self, url):
        """Helper #1 request function"""
        record = requests.get(url, headers={"Accept": "text/x-fasta"})
        return record.text

    def write_fa(self, file_name, fasta_file):
        """Helper #2 write fasta files"""
        p = pathlib.Path("analyses/sp_query/" + file_name).with_suffix('.fa')
        with open(p, 'w') as f:
            f.write(fasta_file)


class FilterSequences:
    """Filter SwissProt fasta files"""

    # Filtering: uncultured and fragment sequences
    def __init__(self):
        self.path = pathlib.Path('analyses/sp_query/').glob('*.fa')

    @property
    def inquery_files(self):
        """Filter fasta files in batch"""
        [self.filter_command(fasta_file) for fasta_file in self.path]

    def filter_command(self, record):
        """Filtering command"""
        filter_subcommand = r'$4!~/[uU]ncultured/ && $4!~/[Ff]ragment/ {print ">"$name$comment"\n"$seq"\n"}'
        command = [
            'bioawk',
            '-c',
            'fastx',
            filter_subcommand,
            str(record),
        ]
        with open(str(record.with_suffix('.filtered')), 'w') as f:
            f.flush()
            subprocess.Popen(command, stdout=f)
