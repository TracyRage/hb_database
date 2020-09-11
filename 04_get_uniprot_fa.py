from db_module.get_uniprot_seqs import GetUniProtMeta
import concurrent.futures
import pathlib
import argparse

# Design argument parser
parser = argparse.ArgumentParser(description="Get fasta data for Blast hits")
parser.add_argument("-in",
                    "--input",
                    required=True,
                    metavar="",
                    help="Provide input directory with all the hits files")
args = parser.parse_args()


# Helper functions
def check_files(directory_path):
    """Make a list of all blast hits files"""
    p = pathlib.Path(directory_path)
    files = [
        str(file.parent) + '/' + str(file.name) for file in p.glob("*/*.txt")
    ]
    return files


def get_meta(file):
    """Init GetUniProtMeta class and run necessary methods"""
    records = GetUniProtMeta(file)
    acc = records.parse_hits
    recs = records.connect_expasy(acc)
    records.write_fa(acc, recs)


def do_multi(files):
    """Process files in parallel"""
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as exe:
        exe.map(get_meta, files)


if __name__ == "__main__":
    files = check_files(args.input)
    print('Start! ')
    print('Working... ')
    do_multi(files)
    print('Done! ')
