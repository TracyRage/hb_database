import re
import pathlib


class FormatSequences:
    """Format SwissProt fasta files"""
    def __init__(self, input_directory):
        p = pathlib.Path(input_directory)
        self.path = p.glob('**/*.[fa|filtered]*')

    @property
    def list_files(self):
        list_of_files = [str(file) for file in self.path]
        return list_of_files

    def format_file(self, file_list):
        [self.read_content(file) for file in file_list]

    def read_content(self, file):
        """Helper #1 format fasta file"""
        pattern_1 = re.compile(r"\{|\}|\[|\]|\(|\)|\;|\||\'|\/|\"|\,")
        with open(file, 'r+') as f:
            content = f.read()
            replaced_text = pattern_1.sub(r'_', content)
            replaced_text_2 = re.sub(r' ', r'_', replaced_text)
            with open(file, 'w') as f:
                f.write(replaced_text_2)
