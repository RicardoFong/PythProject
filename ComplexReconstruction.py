import arguments as arguments
import sys
import os
import argparse
import re
import shutil
from Bio.PDB import *
from Bio import pairwise2
from Bio.PDB import Superimposer
# from Bio.SubsMat import MatrixInfo as matlist # Used for applying subtitution matrix in sequence comparison
from Bio.pairwise2 import format_alignment # DEBUGING Used for printing the alignment. DEBUGING purposes
from statistics import mode
from Bio.PDB import NeighborSearch
from Bio.PDB import PDBIO
import string as string
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa

# Objects
io = PDBIO()
parser = PDBParser(QUIET=True)

########## GET FILE PREFIXES ###########

def get_file_prefix(pdb_files):
        """Tests for match of regex containing .pdb extension, trims it and returns
           the file prefix.
        """

        p = re.compile('(.*).pdb$')
        m = p.search(pdb_files)

        return m.group(1)


######### GET THE SEQUENCES OF EACH FILE ###########

def get_structures(pdb_files, path = arguments.args.inPath):
    """
       This method parses the provided pdb files and extracts and returns the structures
       from them as a list.
    """

    parser = PDBParser(QUIET=True)

    # Main directory
    main_dir = os.getcwd()

    # Change dir to the dir of the files / input dir
    os.chdir(path)

    interactions = []
    if arguments.args.core:
        pdb_files.append(arguments.args.core)


    for file in pdb_files:

        id = get_file_prefix(file)
        interactions.append(parser.get_structure(id, file))
        #print(structure)

        # Set the path again to main

    os.chdir(main_dir)
    return interactions

def get_structure_single_file(file, path = arguments.args.inPath):
    """
       This method parses the provided pdb files and extracts and returns the structures
       from them as a list.
    """

    parser = PDBParser(QUIET=True)

    # Main directory
    main_dir = os.getcwd()

    # Change dir to the dir of the files / input dir
    os.chdir(path)


    id = get_file_prefix(file)
    structure = parser.get_structure(id, file)
        #print(structure)

        # Set the path again to main

    os.chdir(main_dir)
    return structure


######### GET STECHIOMETRY ###########

def get_stechiometry(stechiometry, path = arguments.args.inPath):

    """ If the option stechiometry is selected it returns a dictonary with the stechiometry of its chain (the name of the
        protein as a key and its stechiometry as value), in case a file is provided, otherwise it returns the integer the user provided"""

    stechiometry_match = re.compile(r"^(?P<uniprot>[a-zA-Z0-9]+)\:(?P<chains>[0-9])")
    # Main directory
    main_dir = os.getcwd()

    # Change dir to the dir of the files / input dir
    os.chdir(path)

    stechiometry_info = {}
    if type(stechiometry) is str:
        path = path + '/' + stechiometry
        fp = open(path, 'r')
        for line in fp:
            "".join(line.split())
            m = stechiometry_match.match(line)
            stechiometry_info[m.group("uniprot")] = [m.group("chains")]

        fp.close()

        os.chdir(main_dir)
        return stechiometry_info
    if type(stechiometry) is int:
        os.chdir(main_dir)
        return stechiometry

def get_chains_structure(structure):
        """
        Returns a list of chains given a structure
        """
        chains = []

        for model in structure:
            for chain in model:
                chains.append(chain)
        return chains

def get_pdb_sequence(structure):
    """
    Retrieves the AA sequence from a PDB structure.
    """

    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
    return seq

def get_sequences_string(chain, start = 0, end = None):
    """
    Given a chain it returns its sequence as string
    """

    ppb = PPBuilder()
    pdb_seqs = ppb.build_peptides(chain)
    pp_seq = pdb_seqs[0].get_sequence()
    return pp_seq

def get_atoms_list(chain):
    """
    Given a chain it returns its CA atoms as a list if it is a Protein, or the C4 atoms if it is a nucleic acid
    """
    if molecule_type(chain) == "Protein":
        ca_atoms = []
        for residue in chain:
            if residue.get_id()[0] == " " and residue.has_id("CA"):
                atom = residue["CA"]
                ca_atoms.append(atom)
        return ca_atoms

    elif molecule_type(chain) == "Nucleotide":
        c4_atoms = []
        for residue in chain:
            if residue.get_id()[0] == " " and residue.has_id("C4\'"):
                atom = residue["C4\'"]
                c4_atoms.append(atom)
        return c4_atoms


def define_core_chain(structures):
    """ Expects a generator object containing the  list of structure in a folder.
        returns the chain object in bio.pdb format with most occurrences """

    chains = []
    for structure in structures:
        for model in structure:
            for chain in model:
                chains.append(chain)

    return mode(chains)


def define_core_structure(core_chain, structures):
    """ Expects the core chain and all the structures and returns the core structure """

    if arguments.args.core:
        for structure in get_structures([arguments.args.core], arguments.args.inPath):
            core_structure = structure
    else:
        for structure in structures:
            if core_chain in structure[0].get_chains():
                core_structure = structure
                break

    chain_ids = string.ascii_uppercase
    n = 0
    for core_chain in core_structure[0].get_chains():
        core_chain.id = chain_ids[n]
        n += 1

    return core_structure

def molecule_type(chain):
    """This function retrive if the molecule is a protein, DNA or RNA molecule"""

    DNA = ['DA','DT','DC','DG','DI']
    RNA = ['A','U','C','G','I']

    m_type = ""

    for residue in chain:
        residue_name = residue.get_resname().strip()
        break
    if residue_name in RNA or residue_name in DNA:
        m_type = "Nucleotide"
    else:
        m_type = "Protein"

    return m_type

def get_number_chains_structure(structure):
    """ Given a structure it returns the number of chains it has """

    chains_ids_in_structure = 0
    for core_chain in structure[0].get_chains():
        chains_ids_in_structure += 1


    return chains_ids_in_structure


def get_best_RMSD(core, test):
    """
    Given a core structure and the test structure. Returns the a tuple of the test model with the lower RMSD
    and the superimposer object that will be used to apply the matrix
    If the core chain is a nucleotide it passes through all the sequence to find the best passible superimposition
    with a part of the chain.
    """

    best_RMSD = -1
    best_model = None
    for test_chain in test[0].get_chains():
        test_atoms = get_atoms_list(test_chain)

        for core_chain in core[0].get_chains():
            core_atoms = get_atoms_list(core_chain)

            ini = 0
            test_atoms_length = len(test_atoms)
            if arguments.args.core and molecule_type(core_chain) == molecule_type(test_chain) == "Nucleotide":

                while test_atoms_length <= len(core_atoms):
                    core_atoms_superimpose = core_atoms[ini:test_atoms_length]


                    superimpose = Superimposer()
                    superimpose.set_atoms(core_atoms_superimpose, test_atoms)
                    RMSD = superimpose.rms

                    ini +=1
                    test_atoms_length +=1

                    if RMSD < best_RMSD or best_RMSD == -1 :
                        best_RMSD = RMSD
                        best_model = test[0]
                        superimposer_object_to_apply = superimpose

            elif len(core_atoms) == len(test_atoms) and molecule_type(test_chain) == molecule_type(core_chain):
                # Superimpose chains with same length

                superimpose = Superimposer()
                superimpose.set_atoms(core_atoms, test_atoms)
                RMSD = superimpose.rms

                if RMSD < best_RMSD or best_RMSD == -1:
                    best_RMSD = RMSD
                    best_model = test[0]
                    superimposer_object_to_apply = superimpose

    if best_model:
        return (best_model, superimposer_object_to_apply, test.id)
    else:
        return None

def get_all_RMSD(core, test):
    """ This function takes a test chain and a core chain and returns a list of tuples of the 100 best superpositions between both structures.
        If there are less than 100 it returns as many superpositions as it can find.

        If the core chain is a nucleotide it passes through all the sequence to find the best passible superimposition
        with a part of the chain.
        """
    All_Models = []

    for test_chain in test[0].get_chains():
        test_atoms = get_atoms_list(test_chain)

        for core_chain in core[0].get_chains():
            core_atoms = get_atoms_list(core_chain)

            ini = 0
            test_atoms_length = len(test_atoms)
            if arguments.args.core and molecule_type(core_chain) == molecule_type(test_chain) == "Nucleotide":

                while test_atoms_length <= len(core_atoms):
                    core_atoms_superimpose = core_atoms[ini:test_atoms_length]


                    superimpose = Superimposer()
                    superimpose.set_atoms(core_atoms_superimpose, test_atoms)

                    #test_sequence = get_sequences_string(test_chain, 0, 5)
                    #print(test_sequence)
                    ini +=1
                    test_atoms_length +=1

                    my_model = (test[0], superimpose, test.id)
                    if my_model:
                        All_Models.append(my_model)

            elif len(core_atoms) == len(test_atoms) and molecule_type(test_chain) == molecule_type(core_chain):
            # Superimpose chains with same length

                superimpose = Superimposer()
                superimpose.set_atoms(core_atoms, test_atoms)

                my_model = (test[0], superimpose, test.id)
                if my_model:
                    All_Models.append(my_model)


    All_Models = sorted(All_Models, key=lambda x: x[2])
    return All_Models[0:100]

def get_all_models_to_add(core_structure, structures):

    """ Expects the core structure and a iterator with the pdb files to add. And returns an ordered dictionary
        with:
            key : tuple(test_model, superimposer_object_to_apply)
            value : RMSD
        The dictionary is sorted by RMSD, so the lowest RMSD will be the firsts that will be added to the structure.

    """
    models_to_add = {}
    for structure in structures:
        all_models = get_all_RMSD(core_structure, structure)
        for test_model in all_models:
            if test_model:
                models_to_add[test_model] = test_model[1].rms

    models_to_add = sorted(models_to_add.items(), key=lambda x: x[1])

    return models_to_add


def get_models_to_add(core_structure, structures):

    """ Expects the core structure and a iterator with the pdb files to add. And returns an ordered dictionary
        with:
            key : tuple(best_model, superimposer_object_to_apply)
            value : RMSD
    """

    models_to_add = {}
    for structure in structures:
        test_model = get_best_RMSD(core_structure, structure)
        if test_model:
            models_to_add[test_model] = test_model[1].rms

    models_to_add = sorted(models_to_add.items(), key=lambda x: x[1])


    return models_to_add


def get_macro_complex_iterative(core_structure, interactions, n_chains_added = 0,stechiometry = None, current_stechiometry = None, n_iters = 0, iters_best_rmsd = 5, max_iterations = 50):
    """ Recursive function that builds the model. It takes as input the core structure and all the structures
        of the interactions. First, retrieve the structures to add, based on the structures that have some superimpositon
        with a chain of the model and then adds the chains that are not in the core structure"""


    n_iters += 1
    if arguments.args.verbose:
        sys.stderr.write("ITERATION %i...\n" %(n_iters))

    # Alphabet of the chains
    chain_ids = string.ascii_uppercase
    # Number of chains in the model
    chains_ids_in_core = get_number_chains_structure(core_structure)
    
    if arguments.args.verbose:
        sys.stderr.write("Your model has %i chains so far...\n" %(n_iters))

    if stechiometry:
        # Stechiometry dictionary
        if type(stechiometry) is int:
            total_chains = stechiometry
            models_to_add = get_models_to_add(core_structure, interactions)
        elif type(stechiometry) is dict:
            my_stechiometry = stechiometry
            total_chains = 0
            if arguments.args.verbose:
                sys.stderr.write("Your stechiometry is %s so far...\n" %(current_stechiometry))
            for key in my_stechiometry:
                total_chains += int(my_stechiometry[key][0])
            if n_chains_added < total_chains and n_iters > iters_best_rmsd:
                for protein in current_stechiometry:
                    if current_stechiometry[protein] == int(my_stechiometry[protein][0]):
                        for structure in interactions:
                            if protein in structure.id:
                                interactions.remove(structure)

                models_to_add = get_all_models_to_add(core_structure, interactions)
            else:
                models_to_add = get_models_to_add(core_structure, interactions)
    else:
        models_to_add = get_models_to_add(core_structure, interactions)


    stop = True


    # Iterate through all the models to be added
    for model in models_to_add:

        add_protein = True

        # Get the object model to be added and the matrix that will be applied to the chain
        best_model = model[0][0]
        superimposer_object_to_apply = model[0][1]

        # STECHIOMETRY AS A FILE/DICTIONARY
        if type(stechiometry) is dict:
            match_structure_id = re.compile(r"^(?P<uniprot>[a-zA-Z0-9]+)\.(?P<DNA>[a-zA-Z0-9]+)\.(?P<pdb>[a-zA-Z0-9]+)(\_)(?P<chain1>[a-zA-Z])(\_)(?P<chain1_dna>[a-zA-Z])(?P<chain2_dna>[a-zA-Z])")
            if match_structure_id.match(model[0][2]):
                m = match_structure_id.match(model[0][2])
                uniprot_id_protein = m.group("uniprot")
            else:
                add_protein = False
                uniprot_id_protein = None
            # Iterate through the chains of the model (interaction)
            if n_chains_added == 0:
                current_stechiometry = {}

            if uniprot_id_protein in current_stechiometry:
                if current_stechiometry[uniprot_id_protein] == int(my_stechiometry[uniprot_id_protein][0]):
                    add_protein = False
        elif type(stechiometry) is int:
            if chains_ids_in_core == total_chains:
                stop = True
                break

        # Add the chain if the stechiometry allows it, if
        if add_protein == True:

            for chain in best_model.get_chains():
                superimposer_object_to_apply.apply(chain.get_atoms())
                chain_in_core  = False

            # Iterate through the chains of the core structure. To check if the chain to add is already in the model
                for core_chain in core_structure[0].get_chains():
                    # Looks for clashes of each atom
                    clashes = NeighborSearch(get_atoms_list(core_chain))
                    test_atoms = get_atoms_list(chain)
                    total_clashes = []

                    for atom in test_atoms:
                        atom_clashes = clashes.search(atom.coord, 5)
                        if len(atom_clashes) > 0:
                            total_clashes.extend(atom_clashes) # Add the clashes of each atom to a list

                    threshold = (len(total_clashes)/len(test_atoms))*100 # Threshold as a percentage of the length of the chain

                    # If the chain has more than 30 % of its structure with clashes it is considered as added to the model
                    if threshold > 30:
                        chain_in_core = True
                        break

                # If the chain is not in the model add it
                if chain_in_core == False:
                    # Extract the structure of the chain from the interaction as a separate structure object.
                    # With this we can change its name even if in the same object there is a chain with the new name
                    chain_to_add = chain
                    io.set_structure(chain_to_add)
                    io.save("temporary_chain.pdb")
                    chain_to_add = parser.get_structure("chain", "temporary_chain.pdb")[0][chain_to_add.id]
                    os.remove("temporary_chain.pdb")

                    # If there are less than 26 chains (letters in the alphabet), just set the chain name to the next letter
                    if chains_ids_in_core < len(chain_ids):
                        chain_to_add.id = chain_ids[chains_ids_in_core]
                        chains_ids_in_core += 1

                    # If there are more than 26 chains, then the chain will have two letters
                    elif chains_ids_in_core >= len(chain_ids):
                        first_letter_chain = 0
                        second_letter_chain = 0
                        tries = 0
                        while chain_to_add in core_structure[0].get_chains():
                            chain_to_add.id = chain_ids[first_letter_chain] + chain_ids[second_letter_chain]
                            second_letter_chain += 1
                            tries += 1
                            if tries == len(chain_ids):
                                tries = 0
                                first_letter_chain += 1
                                second_letter_chain = 0
                        chains_ids_in_core += 1

                    # Update the stechiometry dictionary
                    if type(stechiometry) == dict:
                        if uniprot_id_protein not in current_stechiometry:
                            current_stechiometry[uniprot_id_protein] = 1
                        else:
                            current_stechiometry[uniprot_id_protein] += 1

                    # Add the chain to the core structure
                    core_structure[0].add(chain_to_add)
                    n_chains_added += 1
                    if arguments.args.verbose:
                        sys.stderr.write("Chain added to your structure!")

                    stop = False # If we have added a chain we do dot stop the algorithm, because maybe a new chain can be added

    if stechiometry:
        # Do not stop the program if there are still chains to be added.
        if type(stechiometry) is dict:
            if n_chains_added < total_chains:
                stop = False
            if n_chains_added == total_chains:
                stop = True
        elif type(stechiometry) is int:
            chains_ids_in_core = get_number_chains_structure(core_structure)
            if chains_ids_in_core < total_chains:
                stop = False
    if n_iters == max_iterations:
        stop == True

    if stop == True:
        if arguments.args.verbose:
            sys.stderr.write("The main algorithm has finished!")
        return core_structure
    else:
        return get_macro_complex_iterative(core_structure = core_structure, interactions = interactions, n_chains_added = n_chains_added, n_iters= n_iters, stechiometry = stechiometry, current_stechiometry = current_stechiometry)

def remove_nucleotides(output):
    """Takes the output file and remove the atoms wich do not belong to aminoacids"""

    aa = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    if arguments.args.verbose:
        sys.stderr.write("Removing the nucleotides from the output macro-model...\n")

    with open(output, 'r') as f:
        lines = f.readlines()
    with open(output, 'w') as f:
        for line in lines:
            if line.startswith("ATOM") and line[17:20] not in aa:
                pass
            else:
                f.write(line)

def remove_heteroatoms(output):
    """Takes the output file and remove the heteroatoms"""
    
    if arguments.args.verbose:
        sys.stderr.write("Removing heteroatoms from the output macro-model...\n")

    with open(output, 'r') as f:
        lines = f.readlines()
    with open(output, 'w') as f:
        for line in lines:
            if line.startswith("HETATM"):
                pass
            else:
                f.write(line)



if __name__=="__main__":

    if arguments.args.verbose:
        sys.stderr.write("Reading your input files...\n")

    files = arguments.get_files(arguments.args.inPath)
    files_dir = arguments.args.inPath

   # gettin the name of input folder for prosa analysis
    split = files_dir.split("/")
    prosa_name = split[-1]


#######################################################################

if arguments.args.verbose:
        sys.stderr.write("Defining core structure...\n")

core_structure = define_core_structure(define_core_chain(get_structures(files, files_dir)), get_structures(files, files_dir))

if arguments.args.verbose:
        sys.stderr.write("Your core structure is %s\n" %(core_structure.id))

interactions = get_structures(files, files_dir)

stechiometry_file = re.compile(r"^(.*)(.txt$)")

if arguments.args.stechiometry and stechiometry_file.match(arguments.args.stechiometry):
    defined_stechiometry = get_stechiometry(arguments.args.stechiometry)
    if arguments.args.verbose:
        sys.stderr.write("Your stechiometry dictionary is %s\n" %(defined_stechiometry))
elif arguments.args.stechiometry:
    defined_stechiometry = get_stechiometry(int(arguments.args.stechiometry))
    if arguments.args.verbose:
        sys.stderr.write("Your stechiometry is %s\n" %(defined_stechiometry))
else:
    defined_stechiometry = None

######### SAVING THE OUTPUT ##########################

if arguments.args.verbose:
    sys.stderr.write("Building you macro-complex...\n")
final_complex = get_macro_complex_iterative(core_structure, interactions, stechiometry = defined_stechiometry)

if arguments.args.verbose:
    sys.stderr.write("Saving your final structure..")

arguments.save_output(final_complex)

if arguments.args.verbose:
    sys.stderr.write("The complex is built!\n")

#######################################################

io = PDBIO()
io.set_structure(final_complex[0])

# set for outputing the reconstructed complex in Structure folder
#arguments.save_output()
output = f"FinalComplex/Structures/{arguments.args.outfile}.pdb"
#io.save(output)


# Remove the nucleotides
if arguments.args.remove_nucleotides:
    remove_nucleotides(output)

#remove the heteroatoms
if arguments.args.remove_heteroatoms:
   remove_heteroatoms(output)


########################### PROSA ANALYSIS ##################################################
if arguments.args.analyse:
    if arguments.args.verbose:
        sys.stderr.write("Analyzing your struture, please wait...")

    #parser for downloading structure
    input = f"FinalComplex/Structures/{prosa_name}.pdb"

    pdbl=PDBList()

    # Command for downloading the model PDB from RCSB database. it is saved into new PDB folder
    pdbl.retrieve_pdb_file(f'{prosa_name}',pdir='PDB', file_format='pdb')

    prosaPath = "FinalComplex/Analysis/"

    prosa = {"line1":f"read pdb ../Structures/{prosa_name}.pdb obj1 \n", "line2":f"read pdb ../../PDB/pdb{prosa_name}.ent obj2 \n", "line3":"analyse energy *\n", "line5":"shift obj2 2\n", "line6":"color comb obj1 red\n", "line7":"color comb obj2 black\n", "line8":f"graph title {prosa_name} protein reconstruction energy analysis\n", "line9":"pscolor = 1\n", "line10":f"export plot {prosa_name}Plot"}

    with open(prosaPath+f"{prosa_name}Prosa.cmd", 'w') as file:
        for key, value in prosa.items():
            file.write(value)

########### write pdb file withot DNA, this is used for prosa analysis #########
    remove_nucleotides(output)


############# Visualizing #######################################################

    os.chdir(prosaPath) # set working directory in Analysis folder
    print(os.getcwd())
    cmd = f"prosa2003 {prosa_name}Prosa.cmd &" # command for executing prosa.cmd file

    os.system(cmd) # execute command

    cmd1 = f"ps2pdf {prosa_name}Plot.ps" # command for outputing plot prosa file

    os.system(cmd1) # create prosa plot file

    os.chdir("../../") # returning to main directory

    cmd2 = f"chimera PDB/pdb{prosa_name}.ent FinalComplex/Structures/{prosa_name}.pdb &" # command for launching chimera to visualize model and reconstruction

    os.system(cmd2) # launching chimera
else:

    pass


if arguments.args.verbose:
    sys.stderr.write("PROGRAM HAS FINISHED!\n")
