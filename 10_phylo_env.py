from db_module.phylo_env import TreeEnvPrep
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Create phylogenetic
                                     analysis environment""")
    parser.add_argument("-in",
                        "--input_table",
                        required=True,
                        metavar="",
                        help="Provide table with all the fasta file paths")

    args = parser.parse_args()

    # Pipeline per se
    record = TreeEnvPrep(args.input_table)
    print('Creating directories... ')
    record.create_dir
    fasta_tuple = record.init_fa_files
    print('Merge files... ')
    record.merge_fa(fasta_tuple)
    print('Done')
