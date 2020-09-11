from Bio import SeqIO, Entrez
import concurrent.futures
import re
from tqdm import tqdm
import pandas as pd


class FetchAminoacid:
    """Fetch curated sequences from Entrez"""
    def __init__(self, input_table, email):
        """Assign class parameters"""
        self.input_table = input_table
        self.output = "analyses/02_alpha_subunits_aa_sequences.fa"
        self.email = email

    @property
    def parse_table(self):
        """Parse csv table with accession numbers"""
        table = pd.read_csv(self.input_table)
        ids_list = table['accession']
        return ids_list

    def record_replace(self, pattern):
        """Format fasta file, replace [:space:] with \_"""
        string_to_replace = re.compile(r' ')
        replaced_text = string_to_replace.sub(r'_', pattern)
        return replaced_text

    def entrez_inquery(self, accesion_number):
        """Make an Entrez inquery"""
        Entrez.email = self.email
        handle = Entrez.efetch(db='protein',
                               id=accesion_number,
                               rettype='gb',
                               retmode='text')
        record = SeqIO.read(handle, 'genbank')
        # Write a text file w/ aa sequences;
        # later to be used for manual KEGG BLAST inquery
        with open(self.output, 'a') as f:
            f.write('>' + str(record.id) + '_' +
                    self.record_replace(str(record.annotations['organism'])) +
                    '\n' + str(record.seq) + '\n' + '\n')
        handle.close()

    def do_multi(self):
        acc_list = self.parse_table
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as exe:
            list(
                tqdm(exe.map(self.entrez_inquery, acc_list),
                     total=len(acc_list)))
