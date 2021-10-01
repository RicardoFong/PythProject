
import sys
import os
import argparse
import re
import shutil
from Bio.PDB import PDBIO


parser = argparse.ArgumentParser(description = """This program analyzes pdb binary chain interactions provided, 
calculate and reconstructs the polipeptide complex represented by the provided files. """)

parser.add_argument('-i', '--input',
                                        dest = "inPath",
                                        action = "store",
                                        default = os.getcwd(), # Default is current directory
                                        help = """Provide the complete path to the pdb files you
                                                  want to analyze after -i flag. If no path is provided,
                                                  the program will assume that current directory is
                                                  going to be searched for pdb files to analize.""")


parser.add_argument('-o', '--output',
                                        dest = "outfile",
                                        action = "store",
                                        default = None,
                                        help = """The output file will be named after the name provided after this flag.
                                                  This option is mandatory, if no output filename is provided, no analysis
                                                  will be run and this message will be displayed. """)

parser.add_argument('-v', '--verbose',
                                        dest = "verbose",
                                        action = "store_true",
                                        default = False,
                                        help = """This option prints the execution of the program as it is the actions are being.
                                                  performed. It is used without options as it activates the verbose mode on.
                                                  It is adviseable to run the analysis with this option in case any exception is raised.""")

parser.add_argument('-f', '--force',
                                        dest = "force",
                                        action = "store_true",
                                        default = False,
                                        help = """If this option is False and the output directory already exists before the application is
                                                   executed, exit the program execution and warn the user that the directory already exists. 
                                                   If it is True, then the program can continue running and overwrite all the contents of the
                                                   output directory""")

parser.add_argument('-c', '--core',
                                        dest = "core",
                                        action = "store",
                                        default = None,
                                        help = """If this option is None the program will try to find the core structure as the one that contain
                                                  the chain with more interactions. If this option is set, the program will use the as the core chain
                                                  the file that the user has select. The core file must be in the same directory as all the other interactions.""")

parser.add_argument('-s', '--stechiometry',
                                        dest = "stechiometry",
                                        action = "store",
                                        default = None,
                                        help = """If this option is None the program will add as many chains to the complex as chains provided in the interactions.
                                                If stechiometry is set there are two options: 
                                                    stechiometry = number: if stechiometry is a number that will be the number of chains adeded to the complex
                                                    stechiometry = file.txt: if stechiometry is a file, it must be with the format .txt and the number of chains that
                                                    each structure add to the complex is write this way: UNIPROT_id: n_chains (ex: P05412:1)""")

parser.add_argument('-a', '--analyse',
                                        dest = "analyse",
                                        action = "store_true",
                                        default = False,
                                        help = """This option allows to decide if an analysis of the reconstructed structure is to be run.
                                                  when used, the identifier of the reconstructed protein is used in order to download the original
                                                  structure from the RSCB database. This latter structure is used to construct a prosa analysis
                                                  that compares the energy profile of the reconstruction against the crystalography model.
                                                  The prosa.cmd file is built and a pdf containing the energy comparisons is outputed. Both files
                                                  are stored in the Analysis folder within the FinalComplex folder. chimera is launched so both
                                                  structures can be visualized one agaist the other."""
                                        )
                                        
parser.add_argument('-rnuc', '--remove_nucleotides',
                                        dest = "remove_nucleotides",
                                        action = "store_true",
                                        default = False,
                                        help = """If this option is False the macro-complex reconstructed will contain nucleotides in the solution pdb.
                                                  If this option is True, the program will perform a last step in which the atoms wich do not belong to 
                                                  aminoacis will be removed, and therefore no nucleotides will appear in the final solutin pdb file.""")

parser.add_argument('-rhet', '--remove_heteroatoms',
                                        dest = "remove_heteroatoms",
                                        action = "store_true",
                                        default = False,
                                        help = """If this option is False the macro-complex reconstructed may contain heteroatoms in the solution pdb.
                                                  If this option is True, the program will perform a last step in which the heteroatoms will be removed.""") 

args = parser.parse_args()

######### IncorrectInputFile Subclass ###########

class IncorrectInputFile(AttributeError):
    """ Exception when a letter that not belongs to the sequence alphabet is found"""
    __module__ = 'builtins'

    def __init__(self, pdbfile):
        self.pdbfile = pdbfile
    
    def __str__(self):
        """ Raise an exception if the file is not in the correct format"""

        return "The file %s is not appropiate" %(self.pdbfile)

######## Folder Already Exists ##############

class FolderAlreadyExists(OSError):
    """ Exception when the output folder already exists and force has not been selected"""
    __module__ = 'builtins'

    def __init__(self, path):
        self.path = path
    def __str__(self):
        """ Raise an exception if the folder already exists and force has not been selected"""
        return "The path %s already exists. Run force to override it" %(self.path)


####### Functions ########

# READ FILES

def get_files(input):
    """ 
        This method reads the directory that is provided in the command line argument.
        Default directory to read is the current directory. It searches for files with
        .pdb extension and returns them.

    """
    # Print the execution of the program if -v is set
    if args.verbose:
        print("Reading pdb files...")

    # Reading pdb files and saving into a list

    path = input

    # Regular expression to check if the input file name is correct
    input_file = re.compile(r"^(?P<name>[a-zA-Z0-9]+)(\_)(?P<chain1>[a-zA-Z])(\_)(?P<chain2>[a-zA-Z])(.pdb$)")
    input_file2 = re.compile(r"^(?P<uniprot>[a-zA-Z0-9]+)\.(?P<DNA>[a-zA-Z0-9]+)\.(?P<pdb>[a-zA-Z0-9]+)(\_)(?P<chain1>[a-zA-Z])(\_)(?P<chain1_dna>[a-zA-Z])(?P<chain2_dna>[a-zA-Z])(.pdb$)")
    stechiometry_file = re.compile(r"^(.*)(.txt$)")
    # List of pdb files
    pdb_files = []
    # If the user provides a core structure, the file can be named as he/she wants.
    if stechiometry_file.match(args.stechiometry):
        pass
    else:
        IncorrectInputFile(args.stechiometry)
        
    for f in os.listdir(path):
        if input_file.match(f):
            pdb_files.append(f)
        elif input_file2.match(f):
            pdb_files.append(f)
        elif f == args.core:
            pdb_files.append(f)
        elif stechiometry_file.match(args.stechiometry):
            stechiometry = "file"
        else:    
            raise IncorrectInputFile(f)

    return pdb_files

# OUTPUT DIRECTORY 

def save_output(final_complex):
    """
        If the output directory does not exists, create it. 
    """

    # Verbose
    if args.verbose:
        print("Writing the output...")

    subfolder_names = ["Structures", "Analysis"]

    # If the argument force is not selected create the output directory if it does not exists 

    if not args.force:
        if not os.path.exists('FinalComplex'):
            try:
                for subfolder_name in subfolder_names:
                    os.makedirs(os.path.join('FinalComplex', subfolder_name))
            except OSError as err:
                raise err
        else:
            raise FolderAlreadyExists(os.path.join(os.getcwd(),'FinalComplex'))

    # If the argument force is selected create the output directory if it does not exists and, if exists, override it
    elif args.force:
        if not os.path.exists('FinalComplex'):
            for subfolder_name in subfolder_names:
                os.makedirs(os.path.join('FinalComplex', subfolder_name))
        else:
            #Remove the directory
            shutil.rmtree('FinalComplex')
            for subfolder_name in subfolder_names:
                os.makedirs(os.path.join('FinalComplex', subfolder_name))

    io = PDBIO()
    io.set_structure(final_complex[0])

    if args.outfile:
        file_name = args.outfile

    # set for outputing the reconstructed complex in Structure folder
    output = f"FinalComplex/Structures/{file_name}.pdb"
    io.save(output)    

