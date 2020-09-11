from db_module.do_format_seqs import FormatSequences
import argparse

if __name__ == "__main__":
    # Design argument parser
    parser = argparse.ArgumentParser(
        description="Format fasta files for other analyses")
    parser.add_argument(
        "-in",
        "--input_fa",
        required=True,
        metavar="",
        help="Provide input fasta files directory",
    )
    args = parser.parse_args()

    # Formatting per se
    rec = FormatSequences(args.input_fa)
    print('Formatting... ')
    file_list = rec.list_files
    rec.format_file(file_list)
    print('Done! ')
