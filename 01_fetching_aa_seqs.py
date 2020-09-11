from db_module.fetching_protein_seqs import FetchAminoacid
import argparse

# Design the argument parser
parser = argparse.ArgumentParser(
    description="Fetching manually curated aminoacid sequences")
parser.add_argument("-e",
                    "--email",
                    required=True,
                    metavar="",
                    help="Provide user's email (mandatory)")
parser.add_argument("-in",
                    "--input",
                    required=True,
                    metavar="",
                    help="Provide the path of the curated input table")
args = parser.parse_args()

# Get aminoacid sequences for curated alpha subunits
if __name__ == "__main__":
    print("Getting fasta file...\n ")
    data = FetchAminoacid(args.input, args.email)
    data.do_multi()
    print("\nDone\n")
