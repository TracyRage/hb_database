from db_module.generate_best_tree import PerformMuscleAlignment
from db_module.generate_best_tree import PerformTrimal, GenerateTrees
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform MUSCLE and trimal. Generate best tree.")
    parser.add_argument("-in",
                        "--input",
                        required=True,
                        metavar="",
                        help="Provide directory with merged fasta files")

    args = parser.parse_args()

    # Pipeline per se
    record = PerformMuscleAlignment(args.input)
    record.do_alignment
    trimmed = PerformTrimal(args.input)
    print('Trimming... ')
    trimmed.do_trimming
    print('Generating best tree... ')
    rec = GenerateTrees(args.input)
    rec.generate_tree
    print('Done')
