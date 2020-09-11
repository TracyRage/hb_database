from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import concurrent.futures


class PrepBlast:
    """Split input fasta file"""
    def __init__(self, query):
        self.query = query
        self.outp = "analyses/blast_outputs/separated_fasta/"

    @property
    def parse_input(self):
        """Parse query fasta file"""
        for rec in SeqIO.parse(self.query, "fasta"):
            ids, seqs = rec.id, rec.seq
            with open(self.outp + str(ids) + ".fa", 'w') as f:
                f.write(">" + str(ids) + "\n" + str(seqs))


class DoBlast:
    """Perform a BLAST query"""
    def __init__(self, dbase, fa_list):
        self.dbase = dbase
        self.fa_list = fa_list

    def blast_input(self, in_file):
        """Method for getting BLAST hits"""
        cline = NcbiblastpCommandline(
            query="analyses/blast_outputs/separated_fasta/" + in_file,
            db=self.dbase,
            out="analyses/blast_outputs/separated_blast/" +
            in_file.rstrip('.fa') + ".xml",
            evalue=0.0001,
            outfmt=5,
            max_target_seqs=50,
        )
        # Check in case of BLAST debug
        stdout, stderr = cline()
        print(stdout, stderr)

    @property
    def do_multi(self):
        """Do a parallel BLAST query 10x"""
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as exe:
            exe.map(self.blast_input, self.fa_list)
