from Bio.Blast.Applications import NcbirpsblastCommandline
import concurrent.futures
import pathlib
import subprocess
from tqdm import tqdm


class DoRPSBlast:
    """Provide input file for RPS analysis"""
    def __init__(self, fa_input, evalue, cdd_dbase, pfam_dbase, ko_dbase,
                 ko_list, ko_config):
        self.fa_input = fa_input
        self.evalue = evalue
        self.cdd_dbase = cdd_dbase
        self.pfam_dbase = pfam_dbase
        self.ko_dbase = ko_dbase
        self.ko_list = ko_list
        self.ko_config = ko_config

    @property
    def get_files_to_blast(self):
        """Getting list of fasta files to Blast"""
        # Delete empty txt files
        txt = [f for f in pathlib.Path(self.fa_input).glob('**/*.txt')]
        [x.unlink() for x in txt if self.check_for_empty(x)]
        # List all the fasta files
        paths = set(
            [f for f in pathlib.Path(self.fa_input).glob('**/*.fa')]) | set(
                [f for f in pathlib.Path(self.fa_input).glob('**/*.filtered')])
        # Delete empty  fasta files
        [z.unlink() for z in paths if self.check_for_empty(z)]
        # Create a list of non-empty files
        files_to_blast = [str(x) for x in paths]
        return files_to_blast

    def check_for_empty(self, path):
        """Helper #1 check for empty files"""
        return path.is_file() and path.stat().st_size == 0

    def do_cdd_rps(self, record):
        """Method to get Cdd Blast hits"""
        handle = NcbirpsblastCommandline(
            query=record,
            db=self.cdd_dbase,
            out=pathlib.Path(record).with_suffix('.cdd'),
            evalue=self.evalue,
            outfmt="6 qacc sacc",
            max_target_seqs=50,
        )
        stdout, stderr = handle()
        print(stdout, stderr)

    def do_pfam_scan_pl(self, record):
        """Get Pfam hits using Pfam_scan.pl"""
        command = [
            'pfam_scan.pl',
            '-fasta',
            record,
            '-outfile',
            pathlib.Path(record).with_suffix('.pfam'),
            '-dir',
            self.pfam_dbase,
        ]
        subprocess.run(command)

    def do_kofam_scan(self, record):
        """Get KO hits using kofamscan"""
        command = [
            'exec_annotation',
            '-o',
            pathlib.Path(record).with_suffix('.ko'),
            '-p',
            self.ko_dbase,
            '-k',
            self.ko_list,
            '--config',
            self.ko_config,
            '-f',
            'mapper',
            record,
        ]
        subprocess.run(command)

    def do_kofam(self, records):
        """Run kofamscan"""
        list(tqdm(map(self.do_kofam_scan, records), total=len(records)))

    def do_multi(self, function, fa_list):
        """Parallel processing of fasta files"""
        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as exe:
            list(tqdm(exe.map(function, fa_list), total=len(fa_list)))
