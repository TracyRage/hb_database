from db_module.do_regular_blast import PrepBlast, DoBlast
import argparse
import os

# Design argument parser
parser = argparse.ArgumentParser(description="Wrapper for Blastp")
parser.add_argument("-in",
                    "--input",
                    required=True,
                    metavar="",
                    help="Provide your fasta input path")
parser.add_argument("-db",
                    "--database",
                    required=True,
                    metavar="",
                    help="Provide your local database directory path")
args = parser.parse_args()

# Dict with default output paths
default_out = {
    "default_output_path_1": "analyses/blast_outputs/separated_fasta/",
    "default_output_path_2": "analyses/blast_outputs/separated_blast/",
}


def check_for_existance(folder):
    """Check for directories designated for output files"""
    if os.path.exists(folder):
        print('Output default directory already exists!')
    else:
        print('Creating output directory')
        os.makedirs(folder)


if __name__ == "__main__":
    # Make separate fasta files from the input file
    check_for_existance(default_out['default_output_path_1'])
    blast_record = PrepBlast(args.input)
    print('Parsing input file... ')
    blast_record.parse_input
    print('Input fasta file was split!')

    # Check for directory designated for xml blast files
    check_for_existance(default_out['default_output_path_2'])
    # Make a list of files to be BLASTed"
    list_to_blast = [
        filenames
        for filenames in os.listdir(default_out['default_output_path_1'])
    ]

    # BLASTing per se
    blast_inquery = DoBlast(args.database, list_to_blast)
    print('Blast query in process... ')
    blast_inquery.do_multi
    print('Done!')
