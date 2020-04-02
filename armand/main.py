from Bio import SeqIO, PDB, pairwise2
import argparse, os, re
from CustomPDB import CustomChain
from CustomPDB import CustomModel

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""This program builds a macrocomplex.""")

    parser.add_argument('-i', '--input',
                        dest="input",
                        action="store",
                        default=None,
                        help="Input FASTA formatted file or directory")

    parser.add_argument('-o', '--output-file',
                        dest="outfile",
                        action="store",
                        default=None,
                        help="Output file")

    parser.add_argument('-v', '--verbose',
                        dest="verbose",
                        action="store_true",
                        default=False,
                        help="Print log in standard error")

    parser.add_argument('-fa', '--fasta',
                        dest="fasta",
                        action="store",
                        default=None,
                        help="Only sequences having the given regular expression 'pattern' will be in the output")

    parser.add_argument('-r', '--random',
                        dest="random",
                        action="store",
                        default=None,
                        help="Integer defining the number of sequences to be printed in the output. "
                             "If defined, a random selection of the defined size has to be printed.")

    options = parser.parse_args()


def parse_pdb(directory):
    # READ PDBs
    if os.path.isdir(directory):
        p = PDB.PDBParser(PERMISSIVE=False, QUIET=True)
        directory_files = os.listdir(directory)
        if directory_files:
            i = 0
            list_of_objects = []
            for file in directory_files:
                if file.endswith(".pdb"):
                    i += 1
                    name = os.path.splitext(file)[0]
                    model = p.get_structure(name, os.path.join(directory, file))[0]
                    if len(model.child_list) != 2:
                        print("The file {} contains the {} chains. In concrete {} chains. "
                              "Just 2 chains are expected by the program.".format(
                                file, [element for element in model.child_list], len(model.child_list)))
                    else:
                        list_of_objects.append(model)
            if i == 0:
                print("There is no pdb files among the files inside the '{}' directory"
                      .format(os.path.split(directory)[1]))
            else:
                return list_of_objects
        else:
            print("There are no files in '{}'".format(os.path.split(directory)[1]))
    else:
        print("The directory '{}' doesn't exist".format(os.path.split(directory)[1]))


def generate_unique_id(id_list):
    """Returns ID for the chain object."""
    characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for character in characters:
        if character not in id_list:
            return character


def remove_redundant_id(list_of_objects):
    # UNIFY ID'S
    unique_chains = []
    id_list = []
    all_chains = []
    count_dict = {}
    for n in range(len(list_of_objects)):
        pdb = list_of_objects[n]
        print("Model {}".format(n))
        pdb_model = CustomModel(str(n))
        for chain in pdb.get_chains():
            print(chain)
            chain = CustomChain(chain)
            chain.parent = None
            if options.fasta:
                for record in SeqIO.parse(options.fasta, "fasta"):
                    homology = chain.has_homolog(record.seq)
                    if homology:
                        m = re.search("(?<=\:)(.*?)(?=\|)", record.id)
                        fasta_id = m.group()
                        if fasta_id is not None:  #
                            chain.id = fasta_id
                        else:
                            chain.id = record.id
            else:
                if not unique_chains:
                    chain.id = generate_unique_id(id_list)
                    id_list.append(chain.id)
                    unique_chains.append(chain)
                    all_chains.append(chain)
                    print(chain)
                    print("first chain added")
                else:
                    for unique_chain in unique_chains:
                        homology = chain.has_homolog(unique_chain.get_sequence())
                        if homology:
                            chain.id = unique_chain.get_id()
                            all_chains.append(unique_chain)
                            print(chain)
                            print("this chain is NOT unique")
                            break
                        if unique_chain == unique_chains[-1]:
                            chain.id = generate_unique_id(id_list)
                            id_list.append(chain.id)
                            unique_chains.append(chain)
                            all_chains.append(chain)
                            print(chain)
                            print("this chain IS unique")
            pdb_model.add(chain)
        list_of_objects[n] = pdb_model
    for unique_chain in unique_chains:
        unique_chain_count = all_chains.count(unique_chain)
        count_dict[unique_chain] = unique_chain_count
    return count_dict, unique_chains, all_chains


def run(directory):
    list_of_objects = parse_pdb(directory)
    count_dict = remove_redundant_id(list_of_objects)[0]
    unique_chains = remove_redundant_id(list_of_objects)[1]
    all_chains = remove_redundant_id(list_of_objects)[2]
    for model in list_of_objects:
        print(list(model.get_chains()))
    print(all_chains)
    print(count_dict)
    print(unique_chains)

run("./enterovirus")