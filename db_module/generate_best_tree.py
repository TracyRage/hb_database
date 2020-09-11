import subprocess
import pathlib


class PerformMuscleAlignment:
    """Perform multiple alignment using MUSCLE"""
    def __init__(self, input_files):
        self.input_file = pathlib.Path(input_files).glob('**/*.fa')

    @property
    def do_alignment(self):
        """Perform multiple alignment"""
        [self.alignment_helper(str(file)) for file in self.input_file]

    def alignment_helper(self, file):
        """Helper for do_alignment method"""
        command = [
            'muscle',
            '-in',
            file,
            '-out',
            str(pathlib.Path(file).with_suffix('.alignment')),
        ]
        subprocess.run(command)


class PerformTrimal:
    """Alignment trimming"""
    def __init__(self, input_files):
        self.input_file = pathlib.Path(input_files).glob('**/*.alignment')

    @property
    def do_trimming(self):
        """Trim the alignment"""
        [self.trimming_helper(str(file)) for file in self.input_file]

    def trimming_helper(self, file):
        """Helper for do_trimming method"""
        command = [
            'trimal',
            '-in',
            file,
            '-out',
            str(pathlib.Path(file).with_suffix('.trimmed')),
            '-gt',
            '0.05',
        ]
        subprocess.run(command)


class GenerateTrees:
    """Class for phylogenetic tree generation with IQ-TREE"""
    def __init__(self, files):
        self.files = pathlib.Path(files).glob('**/*.trimmed')

    @property
    def generate_tree(self):
        """Method for tree generation"""
        [self.generate_tree_helper(str(file)) for file in self.files]

    def generate_tree_helper(self, file):
        """Helper for generate_tree method"""
        command = [
            'iqtree',
            '-s',
            file,
            '-m',
            'LG+R7',
            '-B',
            '1000',
            '-T',
            'AUTO',
            '-seed',
            '3333',
        ]
        subprocess.run(command)
