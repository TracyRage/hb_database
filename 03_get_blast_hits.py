from db_module.get_regular_blast_hits import GetBlastHits
from db_module.get_regular_blast_hits import PrepEnv
import argparse

# Design argument parser
parser = argparse.ArgumentParser(
    description="Getting BLAST hits from input xml files")
parser.add_argument("-in",
                    "--xml_in",
                    required=True,
                    metavar="",
                    help="Provide xml BLAST file")
args = parser.parse_args()

if __name__ == "__main__":
    # Create the environment and xml files to parse
    record = PrepEnv(args.xml_in)
    record.create_env
    print('Environment has been created... ')
    xml_files = record.get_xml_files
    key_list = record.parse_table

    def get_hits(xml_file, key_word):
        """Get hits from xml files"""
        blast_hits = GetBlastHits(xml_file, key_word)
        parsed_hits = blast_hits.parse_xml
        relevant, merged, output = blast_hits.get_hits(parsed_hits)
        blast_hits.writing_to_files(relevant, merged, output)

    # Getting hits
    [get_hits(xml_file, key) for xml_file, key in zip(xml_files, key_list)]

    print("Done... ")
