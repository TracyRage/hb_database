from db_module.prep_phylogenetics import PhyloPrep
import argparse

if __name__ == "__main__":
    # Design argument parser
    parser = argparse.ArgumentParser(
        description="""Generate a table with all the fasta files for
        subsequent phylogenetic analysis""")
    parser.add_argument("-it",
                        "--initial_table",
                        required=True,
                        metavar='',
                        help="Provide initial csv table")
    parser.add_argument('-nf',
                        '--ncbi_fa',
                        required=True,
                        metavar='',
                        help="Provide manually curated NCBI fasta files")
    parser.add_argument('-bsp',
                        '--sp_blast',
                        required=True,
                        metavar='',
                        help="Provide SwissProt Blast merged fasta files")
    parser.add_argument('-spq',
                        '--sp_query',
                        required=True,
                        metavar='',
                        help="Provide SwissProt query fasta files")
    args = parser.parse_args()

    # Pipeline per se
    record = PhyloPrep(args.initial_table, args.ncbi_fa, args.sp_blast,
                       args.sp_query)
    # Get sorted list of fasta files
    ncbi_fa, blast_fa, query_fa, query_fa_raw = record.sort_files
    # Get sequence number from each fasta file
    blast_fa_count, query_fa_count, query_fa_count_raw = record.get_hits_number(
        blast_fa, query_fa, query_fa_raw)
    # Get indexes for dataframes
    index_1, index_2 = record.get_index
    # Generate fasta files dataframes
    sp_df, query_df = record.generate_df(index_1, ncbi_fa, blast_fa,
                                         blast_fa_count, index_2, query_fa,
                                         query_fa_count, query_fa_count_raw)
    # Generate final table
    final_table = record.create_final_table(sp_df, query_df)
    print('Table has been generated! ')
