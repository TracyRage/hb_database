from Bio.Blast import NCBIXML
import re
import pathlib
import pandas as pd


class PrepEnv:
    """Prepare environment for subsequent procedures"""
    def __init__(self, xml_files):
        self.xml_files = sorted(pathlib.Path(xml_files).glob("**/*.xml"))
        self.input_table = pd.read_csv('data/seqs/01_list_alpha_subunits.csv')

    @property
    def parse_table(self):
        """Get keywords"""
        accession_list = self.input_table.sort_values(
            by='accession')['alpha_subunits']
        return list(accession_list)

    @property
    def get_xml_files(self):
        xml_files = [str(x) for x in self.xml_files]
        return xml_files

    @property
    def create_env(self):
        """Create output paths"""
        p1 = pathlib.Path('analyses/blast_outputs/separated_hits/merged/')
        p2 = pathlib.Path('analyses/blast_outputs/separated_hits/relevant/')
        p3 = pathlib.Path('analyses/blast_outputs/separated_hits/outgroup')
        paths = [p1, p2, p3]
        [self.create_env_helper(x) for x in paths]

    def create_env_helper(self, path):
        """Helper for create_env method"""
        path.mkdir(parents=True, exist_ok=True)


class GetBlastHits:
    """Getting relevant and outgroup BLAST hits"""
    def __init__(self, xml_file, key_word):
        self.xml_file = xml_file
        self.key_word = key_word

    @property
    def parse_xml(self):
        """Parse xml BLAST input file"""
        result_handle = open(self.xml_file)
        blast_record = NCBIXML.read(result_handle)
        return blast_record.alignments

    def get_hits(self, input_record):
        """Get relevant and outgroup hits based on key words"""
        relevant_records = [
            rec for rec in input_record if self.match_key(rec.title)
        ]

        outgroup_records = [
            rec for rec in input_record if not self.match_key(rec.title)
        ]

        merged_records = [rec for rec in input_record]

        return relevant_records, merged_records, outgroup_records

    def match_key(self, line):
        """Helper function #2"""
        pattern_1 = re.compile(self.key_word)
        return pattern_1.search(line, re.IGNORECASE)

    def writing_to_files(self, relevant, merged, output):
        """Write records of interest to files"""
        relevant_path = pathlib.Path(
            "analyses/blast_outputs/separated_hits/relevant/")
        merged_path = pathlib.Path(
            'analyses/blast_outputs/separated_hits/merged/')
        output_path = pathlib.Path(
            'analyses/blast_outputs/separated_hits/outgroup')

        self.writing_helper(str(relevant_path), relevant)
        self.writing_helper(str(merged_path), merged)
        self.writing_helper(str(output_path), output)

    def writing_helper(self, output_path, record_to_be_written):
        """Helper function #3"""
        with open(
                output_path + '/' + pathlib.Path(self.xml_file).stem + '.txt',
                'w') as f:
            [f.write(str(rec) + "\n") for rec in record_to_be_written]
