from db_module.do_rps_blast import DoRPSBlast
import argparse

if __name__ == "__main__":
    # Design argument parser
    parser = argparse.ArgumentParser(
        description="Getting Cdd, Pfam, Cog and Ko hits")
    parser.add_argument(
        "-in",
        "--input_fa",
        required=True,
        metavar="",
        help="Provide input fasta file",
    )
    parser.add_argument(
        "-cdd",
        "--db_cdd",
        required=True,
        metavar="",
        help="Provide Cdd database path",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        required=True,
        metavar="",
        help="Provide Evalue",
    )

    parser.add_argument(
        "-pfam",
        "--db_pfam",
        required=True,
        metavar="",
        help="Provide Pfam database path",
    )

    parser.add_argument(
        "-ko",
        "--db_ko",
        required=True,
        metavar="",
        help="Provide Ko database path",
    )

    parser.add_argument(
        "-kl",
        "--ko_list",
        required=True,
        metavar="",
        help="Provide Ko list file path",
    )

    parser.add_argument(
        "-kc",
        "--ko_config",
        required=True,
        metavar="",
        help="Provide Ko config file path",
    )

    args = parser.parse_args()

    # Do RPS Blast
    rec = DoRPSBlast(
        args.input_fa,
        args.evalue,
        args.db_cdd,
        args.db_pfam,
        args.db_ko,
        args.ko_list,
        args.ko_config,
    )
    fa_files = rec.get_files_to_blast
    print('Getting Cdd accession numbers... + \n')
    rec.do_multi(rec.do_cdd_rps, fa_files)
    print('\n Cdd annotations fetched... + \n')
    print('Getting  Pfam accession numbers... + \n')
    rec.do_multi(rec.do_pfam_scan_pl, fa_files)
    print('\n Pfam annotations fetched... \n')
    print('\n Get Ko annotation \n')
    rec.do_kofam(fa_files)
    print('\n Done!')
