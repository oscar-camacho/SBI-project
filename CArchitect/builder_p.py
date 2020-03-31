#!/usr/bin/env python
from Bio import SeqIO, PDB, pairwise2
from pdb_classes import CustomChain, CustomModel
from exception_classes import *
import utilities
import argparse, os, sys, re, random, copy


                        #EDITAR ESTO
parser = argparse.ArgumentParser(description="""ComplexArchitect is a Python application designed to generate macrocomplex structures from simple pair inetractions PDB files.""")

parser.add_argument('-i', '--input',
                        dest = "input",
                        action = "store",
                        default = None,
                        help = "Input directory where PDB files with the pair interactions are located.")

parser.add_argument('-o', '--output',      #Esto fue cambiado
                        dest = "outfile",
                        action = "store",
                        default = "complex",
                        help = "Name of the output file. No extension is needed. The output will be saved on the current working directory.")

parser.add_argument('-fa', '--fasta',
                        dest = "fasta",
                        action = "store",
                        default = None,
                        help = """FASTA file with the sequences of the chains that will conform the macrocomplex. They have to correspond to the sequences
                                of the chains from the PDB files. The file should contain unique sequences; sequences don't need to be repeated.""")

parser.add_argument('-v', '--verbose',
                        dest = "verbose",
                        action = "store_true",
                        default = False,
                        help = "Shows the progress of the program.")

parser.add_argument('-r', '--random',
                        dest = "random",
                        action = "store",
                        default = None,
                        help = """Integer defining the number of sequences to be printed in the output. If defined, a random selection of
                                the defined size has to be printed.""")

options = parser.parse_args()


def data_extraction(directory, verbose=False):
    """Takes a directory with PDB files, generates PDB objects of the models for each file and returns a list of the PDB models.

    Keyword arguments:
    directory -- path of the directory where the PDB files are
    verbose -- boolean, prints to stderr the progress of the program

    Exceptions raised:
    Directory_Not_Found -- raised if name of the directory does not exist
    No_Input_Directory_Provided -- raised if directory is no provided
    No_PDB_Found -- raised if directory does not contain PDB files
    Incorrect_Number_Chains -- raised if a PDB object contains more than 2 chains"""
    if directory is not None:
        if verbose:
            sys.stderr.write("Looking for PDB input files from %s\n" %options.input)
        if os.path.isdir(options.input):
            path = options.input + '/'
            pdb_files = [path + pdb_file for pdb_file in os.listdir(path) if pdb_file[-4:] == '.pdb']
        else:
            raise Directory_Not_Found(options.input)
    else:
        raise No_Input_Directory_Provided("""No input directory provided. You should introduce the name of the directory containing the PDB files for each interacting pair of the complex.""")

    num_files = len(pdb_files)   #Printing number of processed filenames
    if verbose:
        sys.stderr.write('%s PDB files found.\n' %(num_files))

    parser = PDB.PDBParser(PERMISSIVE=1, QUIET=True)   #Storing PDB information
    if verbose:
        sys.stderr.write("Storing PDB information...\n")
    if num_files != 0:
        pdb_models = []
        for pdb_file in pdb_files:
            pdb_models.append(parser.get_structure('pdb_model', pdb_file))
            if verbose:
                sys.stderr.write("  %s finished\n" %(pdb_file))
    else:
        raise No_PDB_Found("No PDB files found. Please make sure the given directory contains PDB files.")

    for pdb_model in pdb_models:    #Ensuring all pdb files contain 1 or 2 chains
        if len(pdb_model.child_list) > 2:
            raise Incorrect_Number_Chains("A PDB input file does not contain two chains. Please, all PDB files must only contain two chains.")
    if verbose:
        sys.stderr.write("PDB information stored.\n\n")
    return pdb_models


def new_id(id_list):
    """Checks for characters in an ASCII characters set that are not in a given list of IDs. Returns a new character that will serve as chain ID.

    Keyword arguments:
    id_list -- list of chain IDs that is progressively being apdated"""
    characters = utilities.ASCII_characters
    for character in characters:
        if character not in id_list:
            return character


def find_equivalent_chains(pdb_models, verbose=False):
    """Updates original PDB models to customized PDB models and return the list of those. Unifies the ID chains in a way that those chains with sequences
    that are very similar or identical (equivalent chains) receive the same ID.

    Keyword arguments:
    pdb_models -- list of objects of the pdb models with 2 interacting chains, extracted from the PDB files
    verbose -- boolean, prints to stderr the progress of the program

    Considerations: the function works in two different scenarios.
    FASTA file is provided -- pairwise alignments between chain sequences and FASTA sequences are performed. If they are equivalent, the chain receives the ID from the FASTA sequence.
    FASTA file is not provided -- pairwise alignments are performed between the chain sequences. Unique chains are identified and those chains that are equivalent share the same ID."""
    if verbose:
        sys.stderr.write("Identifying equivalent chains and unifying IDs.\n")
        sys.stderr.write("Performing pairwise sequence alignments...\n")
    unique_chains = []
    id_list = []
    for i in range(len(pdb_models)):
        pdb = pdb_models[i]
        pdb_model = CustomModel(str(i))
        for chain in pdb.get_chains():
            chain = CustomChain(chain)
            chain.parent = None
            if options.fasta:
                for record in SeqIO.parse(options.fasta, "fasta"):
                    equivalence = chain.has_equivalent(record.seq)
                    if equivalence:
                        m = re.search("(?<=\:)(.*?)(?=\|)", record.id)
                        id = m.group()
                        if id is not None:          #
                            chain.id = id
                        else:
                            chain.id = record.id
            else:
                if not unique_chains:
                    chain.id = new_id(id_list)
                    id_list.append(chain.id)
                    unique_chains.append(chain)
                else:
                    for unique_chain in unique_chains:
                        equivalence = chain.has_equivalent(unique_chain.get_sequence())
                        if equivalence:
                            chain.id = unique_chain.get_id()
                            break
                        if unique_chain == unique_chains[-1]:
                            chain.id = new_id(id_list)
                            id_list.append(chain.id)
                            unique_chains.append(chain)
            pdb_model.add(chain)
        pdb_models[i] = pdb_model
    if verbose:
        sys.stderr.write("Equivalent chains correctly identified and IDs unified.\n")
    return pdb_models






if options.fasta is None:
    sys.stderr.write("ATENTION! You did not provide a FASTA file. You may be required to specify the stoichometry during the course of the program.\n\n")

try:
    pdb_models = data_extraction(options.input, options.verbose)
except (Directory_Not_Found, No_Input_Directory_Provided, No_PDB_Found, Incorrect_Number_Chains) as e:
    print(e, file=sys.stderr)
    sys.exit(1)

updated_pdb_models = find_equivalent_chains(pdb_models, options.verbose)

#for pdb_model in updated_pdb_models:
#    for chain in pdb_model.get_chains():
#        print(chain.id)




def data_transformation(pdb_models, verbose):
    interaction_dict = {}
    count = 0
    for pdb_model in pdb_models:
        for chain in pdb_model.get_chains():
            count += 1
            if count%2 == 0:
                chain_index = 0
            else:
                chain_index = 1
            interaction_dict.setdefault(chain.id, [])
            if chain.id in interaction_dict:
                interacting_chain = [other_chain for other_chain in pdb_model.get_chains()][chain_index]
                interaction_dict[chain.id].append((pdb_model, chain, interacting_chain))
    return interaction_dict

interaction_dict = data_transformation(updated_pdb_models, options.verbose)
for key in interaction_dict:
    print (interaction_dict[key])


def starting_model(pdb_models, verbose):
    chain_count = {}
    max_value = 0
    possible_models = []
    for pdb_model in pdb_models:
        for chain in pdb_model.get_chains():
            if chain.id not in chain_count:
                chain_count[chain.id] = 0
            chain_count[chain.id] += 1
    for chain, value in chain_count.items():
        if value > max_value:
            max_value = value
            max_interacting_chain = chain
    for pdb_model in pdb_models:
        for chain in pdb_model.get_chains():
            if chain.id == max_interacting_chain:
                possible_models.append(pdb_model)
    starting_model = random.choice(possible_models)
    return starting_model


#starting_model = starting_model(updated_pdb_models, options.verbose)
#print (starting_model)


def generate_model_profile(model):
    """Generates a dictionary with the id chain as key and the number of repetitions of this chain as values"""
    profile = {}  # { "A":1, "B":4, ...}
    for chain in model:
        profile.setdefault(chain.id, 0)  # If id not in dic, set it to 0
        profile[chain.id] += 1 # Sum 1 to the counter of the id
    return profile


def has_clashes(move_atoms, model):
    """Compares the atoms backbone atoms of the moving chain with the backbone atoms of the model"""
    backbone = {"CA", "C1\'"}
    chain_atoms = [atom for atom in move_atoms if atom.id in backbone]  # Gets only the backbone atoms
    model_atoms = [atom for atom in model.get_atoms() if atom.id in backbone]
    ns = PDB.NeighborSearch(model_atoms)  # Generates a neigbour search tree to speed up distance calculations
    clashes = 0
    for atom in chain_atoms:
        clashes += bool(ns.search(atom.coord, 2))  # If this atom shows clashes, add 1 to the clashes counter
    if clashes/len(chain_atoms) >= 0.03:  # If more than 3% of atoms show clashes return yes
        return True
    else:  # Otherwise return no
        return False


def complex_builder(interaction_dict, pdb_models, num_models, max_chains, verbose):
    for i in range(1, num_models + 1):
        if verbose:
            sys.stderr.write("Building Macrocomplex " + str(i) + " ...\n")
        macrocomplex = starting_model(pdb_models, verbose).copy()
        model_stech = generate_model_profile(macrocomplex)
        macrocomplex.id = "Model_" + str(i)
        run = True  # While this variable is true, the program will keep trying to add chains to the macrocomplex
        num_of_chains = 2  # The model starts with 2 chains already
        num_empty_chains = 0  # NUmber of chains that have all their interactions depleted
        #while run:
        for chain in macrocomplex:  # Iterates the macrocomplex chains
            if num_of_chains < max_chains:  # If the number of chains still hasn't reached the maximum allowed
                if len(interaction_dict[chain.id]) != 0:    # If this chain still has pending interactions. Chain id is the key of the dictionary
                    random.shuffle(interaction_dict[chain.id])   # Shuffle the interactions list (to avoid repetitive behaviour)
                    for tuple in interaction_dict[chain.id]:
                        fix = tuple[1]    # Get the chain instance that corresponds to the same chain in macrocomplex
                        to_move = tuple[2]   # Get the chain instance that interacts with the other
                        sup = PDB.Superimposer()  # Generates a superimposer instance
                        chain_atoms, fix_atoms = chain.get_common_atoms(fix) # Get common atoms between the
                                                                                # macrocomplex chain and the one in the interaction dictionary
                        sup.set_atoms(chain_atoms, fix_atoms)  # Generate the superposition
                        move = to_move.copy()  # Make a copy of the chain to move
                        sup.apply(move)  # Apply superposition matrix
                        move_atoms = sorted(move.get_atoms())
                        if not has_clashes(move_atoms, macrocomplex):
                            if verbose:
                                sys.stderr.write("  Succesful superposition between " + str(chain.id) + " from macrocomplex and chain " + fix.id + " from model " + tuple[0].id + ".\n")
                                sys.stderr.write("  Chain " + str(num_of_chains) + " added: Chain " + move.id + " from model " + tuple[0].id + ".\n")
                            #move.parent = None      # Sets the parent to none to evade biopython's strict id policy
                            #macrocomplex.add(move)  # Adds the target chain to the model
                            model_stech.setdefault(move.id, 0)
                            model_stech[move.id] += 1
                            num_of_chains += 1
                            del interaction_dict[chain.id][0]   # elimino la tupla con el modelo que acabamos de añadir. faltaria eliinar la tupla con el mismo objeto.
                        else:
                            if verbose:
                                sys.stderr.write("  Unsuccesful superposition between " + str(chain.id) + " from macrocomplex and chain " + fix.id + " from model " + tuple[0].id + ".\n")
                                sys.stderr.write("  Chain " + move.id + " from model " + tuple[0].id + " NOT ADDED.\n")
                        model_id = tuple[0].get_id()    #para eliminar la otra tupla quizá hay que identificarla por la id del modelo y luego mirar como eliminarla
                        print(model_id)




# Will have to remove from interaction_dict the models with the model id
#for pdb_model in pdb_models:
#    if pdb_model.get_id() == str(1):
#        print("hi")



complex_builder(interaction_dict, pdb_models, 1, 100, options.verbose)










#def starting_model2(interacting_dict, verbose):
#    max_interactions = 0
#    for key, value in interaction_dict.items():
#        if len(value) > max_interactions:
#            len(value) = max_interactions
#            starting_model = key
#    return starting_model

#a=starting_model2(interacting_dict, options.verbose)
#print(a)





















#def build_complex(input, fasta, verbose):
#    if fasta is None:
#        sys.stderr.write("ATENTION! You did not provide a FASTA file. You may be required to specify the stoichometry during the course of the program.\n\n")

#    try:
#        pdb_models = data_extraction(input, verbose)
#    except (Directory_Not_Found, No_Input_Directory_Provided, No_PDB_Found, Incorrect_Number_Chains) as e:
#        print(e, file=sys.stderr)
#        sys.exit(1)

#    updated_pdb_models = find_equivalent_chains(pdb_models, verbose)

#    for pdb_model in updated_pdb_models:
#        for chain in pdb_model.get_chains():
#            print(chain.id)



#build_complex(options.input, options.fasta, options.verbose)





#FASTA
