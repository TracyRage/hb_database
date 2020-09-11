from db_module.get_references import GetReference
import argparse

if __name__ == "__main__":

    # Design argument parser
    parser = argparse.ArgumentParser(
        description='Get references for manually chosen sequences')
    parser.add_argument("-in",
                        "--input_table",
                        metavar='',
                        help='Provide annotated table')
    parser.add_argument("-e", "--email", metavar="", help='Provide email')

    args = parser.parse_args()

    # Pipeline per se
    print('Getting reference data... ')
    record = GetReference(args.email, args.input_table)
    # Get accession numbers
    accession_list = record.parse_table
    # Get reference list
    reference_list = record.get_ref(accession_list)
    primary_reference = record.get_primary_ref(reference_list)
    reference_dataframe = record.generate_dataframe(accession_list,
                                                    primary_reference)
    # Create annotated table
    record.merge_with_input_table(reference_dataframe)
    print('Done')
