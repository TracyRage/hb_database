from db_module.do_sw_query import GetSwissProt, FilterSequences
import argparse

if __name__ == "__main__":
    # Design arguments
    parser = argparse.ArgumentParser(
        description="Access UniProt API for fasta files and filter them")
    parser.add_argument(
        "-in",
        "--input_file",
        required=True,
        metavar="",
        help="Provide input csv table",
    )
    args = parser.parse_args()

    # Process per se
    rec = GetSwissProt(args.input_file)
    rec.check_output_dir
    table = rec.init_table
    df_enzyme, df_subunit = rec.get_dfs(table)
    print('Getting fasta files... ')
    fa_files = rec.sp_query(df_subunit, df_enzyme)
    print('Remove fragments and uncultured sequences... ')
    rec = FilterSequences()
    rec.inquery_files
    print('Done! ')
