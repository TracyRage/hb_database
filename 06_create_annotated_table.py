from db_module.create_annotated_table import CreateAnnotationTable
import argparse

if __name__ == "__main__":

    # Design argument parser
    parser = argparse.ArgumentParser(description='Getting annotated table')
    parser.add_argument("-in",
                        "--fasta_input",
                        required=True,
                        metavar="",
                        help="Provide directory w/ cdd, ko, pfam files")
    args = parser.parse_args()

    # Get annotated table
    print('Getting files... ')
    rec = CreateAnnotationTable(args.fasta_input)
    cdd_list, ko_list, pfam_list = rec.create_file_list
    df_list = rec.chain_df(cdd_list, ko_list, pfam_list)
    print('Creating tables... ')
    dfs = rec.create_table(df_list)
    rec.write_dfs(cdd_list, dfs)
    print('Done')
