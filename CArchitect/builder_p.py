#!/usr/bin/env python
from Bio import SeqIO, PDB, pairwise2
from pdb_classes import CustomChain, CustomModel
from exception_classes import *
import utilities
import argparse, os, sys, re


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



def build_complex(input, fasta, verbose):
    if fasta is None:
        sys.stderr.write("ATENTION! You did not provide a FASTA file. You may be required to specify the stoichometry during the course of the program.\n\n")

    try:
        pdb_models = data_extraction(input, verbose)
    except (Directory_Not_Found, No_Input_Directory_Provided, No_PDB_Found, Incorrect_Number_Chains) as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    updated_pdb_models = find_equivalent_chains(pdb_models, verbose)

    for pdb_model in updated_pdb_models:
        for chain in pdb_model.get_chains():
            print(chain.id)



build_complex(options.input, options.fasta, options.verbose)



#FASTA
