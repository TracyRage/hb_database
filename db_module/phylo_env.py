import pandas as pd
import pathlib
from Bio import SeqIO
import re


class TreeEnvPrep:
    """Class for creating phylogenetic analyses procedures
    environment"""
    def __init__(self, input_table):
        self.input_table = pathlib.Path(input_table)
        self.table_init = pd.read_csv(self.input_table, sep='\t')

    @property
    def create_dir(self):
        """Create directory for each analyzed gene"""
        dirname = self.table_init['enzyme']
        [self.create_dir_helper(enzyme) for enzyme in dirname]

    def create_dir_helper(self, basename):
        """Helper for create_dir method"""
        path = pathlib.Path("analyses/phylo/" + basename)
        path.mkdir(parents=True, exist_ok=True)

    @property
    def init_fa_files(self):
        """Init to be processed fasta files"""
        ncbi_ref = self.table_init['ncbi_ref']
        sp_ref = self.table_init['sp_ref']
        sp_query = self.table_init['sp_query']
        dirname = self.table_init['enzyme']
        fasta_tuples = [
            (a, x, y, z)
            for a, x, y, z in zip(dirname, ncbi_ref, sp_ref, sp_query)
        ]
        return fasta_tuples

    def merge_fa(self, fasta_tuples):
        """Merge fasta files for phylogenetic analysis"""
        [self.merge_fa_helper(fasta_tuple) for fasta_tuple in fasta_tuples]

    def merge_fa_helper(self, fasta_tuple):
        """Helper for merge_fa method"""
        path = 'analyses/phylo/' + fasta_tuple[0] + '/' + fasta_tuple[0] + '.fa'
        with open(path, 'a') as f:
            for file in fasta_tuple[1:]:
                handle = [sequence for sequence in SeqIO.parse(file, 'fasta')]
                [
                    f.write('>' + self.format_header(str(seqq.id)) + '\n' +
                            str(seqq.seq) + '\n' + '\n') for seqq in handle
                ]

    def format_header(self, header):
        """Format headers"""
        pattern_blast = re.compile(r'(.*?)_(\w+?)_.*?(Name=)(.*?)_.*')
        pattern_sp = re.compile(
            r'[tr_|sp_]+([a-zA-Z0-9].*?)_.*OS=([A-Za-z0-9].*?)_.*GN=([A-Za-z0-9].*?)_.*'
        )
        pattern_ncbi = re.compile(r'([A-Za-z0-9].*\.\d)_([a-zA-Z0-9].*?)_.*')

        if pattern_blast.search(header):
            blast_header = pattern_blast.sub(r'\1-O=\2-G=\4-T=BLAST', header)
            return blast_header
        elif pattern_sp.search(header):
            sp_header = pattern_sp.sub(r'\1-O=\2-G=\3-T=Query', header)
            return sp_header
        elif pattern_ncbi.search(header):
            ncbi_header = pattern_ncbi.sub(r'\1-O=\2-NCBI', header)
            return ncbi_header
        else:
            return header
