#!/usr/bin/env python
from Bio import SeqIO, PDB, pairwise2
import argparse
import os
import sys
from CustomPDB import CustomChain
from CustomPDB import CustomModel

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="""This program builds a macrocomplex.""")

    parser.add_argument('-i', '--input',
                        dest = "input",
                        action = "store",
                        default = None,
                        help = "Input FASTA formatted file or directory")

    parser.add_argument('-o', '--output-file',
                        dest = "outfile",
                        action = "store",
                        default = None,
                        help = "outputfile")

    parser.add_argument('-v', '--verbose',
                        dest = "verbose",
                        action = "store_true",
                        default = False,
                        help = "Print log in stderr")

    parser.add_argument('-p', '--pattern',
                        dest = "pattern",
                        action = "store",
                        default = None,
                        help = "Only sequences having the given regular expression 'pattern' will be in the output""")

    parser.add_argument('-r', '--random',
                        dest = "random",
                        action = "store",
                        default = None,
                        help = """Integer defining the number of sequences to be printed in the output. If defined, a random selection of
                                the defined size has to be printed.""")

    options = parser.parse_args()

#Storing pdb filenames into a list:
if options.input:    #specification of input
    if options.verbose:
        sys.stderr.write("Looking for PDB input files from %s\n" %options.input)
    if os.path.isdir(options.input):
        path = options.input + '/'
        pdb_files = [path + pdb_file for pdb_file in os.listdir(path) if pdb_file[-4:] == '.pdb']
    else:
        raise Exception("Directory %s does not exist. Please select a valid directory." %options.input)
else:  #no specification of directory. Look in the current directory
    raise Exception("No input provided.")

#Printing number of processed filenames
num_files = len(pdb_files)
if options.verbose:
    sys.stderr.write('%s PDB files found.\n' %(num_files))

parser = PDB.PDBParser()
if options.verbose:
    sys.stderr.write("Storing PDB information...\n")
if num_files != 0:
    pdb_models = []
    for pdb_file in pdb_files:
        try:
            pdb_models.append(parser.get_structure('pdb_model', pdb_file))
            if options.verbose:
                sys.stderr.write("  %s finished\n" %(pdb_file))
        except Exception:   #Define exception type
            sys.stderr.write("PDB files couldn't be opened. Please, revise that their format is correct.")
else:
    raise Exception("No PDB files found. Please make sure the given directory contains PDB files.")

#Ensuring all pdb files contain 1 or 2 chains
for pdb_model in pdb_models:
    if len(pdb_model.child_list) > 2:
        raise Exception("A PDB input file does not contain two chains. Please, all PDB files must only contain two chains.")
if options.verbose:
    print("PDB information stored.")


#También podemos pasarle las unique chains y que saque cada vez las ids...
def new_id(id_list):
    """Returns ID for the chain object."""
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for letter in letters:
        if letter not in id_list:
            return letter


# Si hemos encontrado homólogo, quizá no me interesa añadirlo a la lista. Me interesa tener las unique chains. Alternativa es el set/diccionario
unique_chains = []
id_list = []
for i in range(len(pdb_models)):
    pdb = pdb_models[i]
    pdb_model = CustomModel(str(i))
    for chain in pdb.get_chains():
        chain = CustomChain(chain)
        chain.parent = None
        if not unique_chains:
            chain.id = new_id(id_list)
            id_list.append(chain.id)
            unique_chains.append(chain)
        else:
            for unique_chain in unique_chains:
                homology = chain.has_homolog(unique_chain)
                if homology:
                    chain.id = unique_chain.get_id()
                    break
                if unique_chain == unique_chains[-1]:
                    chain.id = new_id(id_list)
                    id_list.append(chain.id)
                    unique_chains.append(chain)
        pdb_model.add(chain)
    pdb_models[i] = pdb_model



#for unique_chain in unique_chains:
#    print (unique_chain.id)


#for pdb_model in pdb_models:
#    for chain in pdb_model.get_chains():
#        print(chain.id)
