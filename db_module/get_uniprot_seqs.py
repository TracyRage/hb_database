from Bio import ExPASy, SwissProt
import re


class GetUniProtMeta:
    """Getting UniProt metadata for Blast hits"""
    def __init__(self, file_path):
        self.file_path = file_path

    @property
    def parse_hits(self):
        """Parse Blast hits files and get accession#"""
        with open(self.file_path) as f:
            text = f.read()
            acc = self.get_accession(text)
        return acc

    def get_accession(self, text):
        """Helper function #1; getting accession#"""
        pattern = re.compile(r'(?<=\|)([OPQ]*[0-9A-Z]*)(?=\|)')
        return pattern.findall(text, re.IGNORECASE)

    def connect_expasy(self, accession_list):
        """Connect to ExPASy and get UniProt records"""
        handles = [ExPASy.get_sprot_raw(acc) for acc in accession_list]
        records = [SwissProt.read(handle) for handle in handles]
        return records

    def write_fa(self, accession_list, records):
        """Write to a fasta file"""
        # Format the fasta header
        format_fasta = [
            '>' + acc + ' ' + rec.organism + rec.gene_name + '\n' +
            rec.sequence for acc, rec in zip(accession_list, records)
        ]
        # Write to a file
        with open(re.sub(r'.txt', '.fa', self.file_path), 'w') as f:
            [f.write(str(r) + '\n\n') for r in format_fasta]
