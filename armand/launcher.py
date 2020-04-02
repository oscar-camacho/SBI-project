#!/usr/bin/env python


def read_pdb_files(directory, verbose=False):
    """Reads the input directory and generates pdb models"""
    if verbose:
        print("Reading pdb input files from %s" % directory)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    if os.path.isdir(directory) and directory.endswith("/"):
        try:
            pdbmodels = [parser.get_structure("Model_pair", directory + f)[0] for
                    f in listdir(directory) if f.endswith(".pdb")] #  Generates pdb objects for files that end with .pdb
        except:
            sys.stderr.write("PDB files couldn't be opened. Please, revise that their format is correct.")
            sys.exit(1) #personalized exception
    else:
        sys.stderr.write("Directory %s doesn't exists, please select a valid directory." % directory)
        sys.exit(1)  #personalized exception
    if not bool(pdbmodels):  # If no pdb instance is generated
        sys.stderr.write("No pdb files where read. Please make sure the given directory contains pdb files. ")
        sys.exit(1)
    for model in pdbmodels:
        if len(model.child_list) != 2:
            sys.stderr.write("A pdb input file doesn't contains two chains. Please, all input pdbs must only contain "
                             "two chains.")
            sys.exit(1)
    if verbose:
        print("Pdb objects stored")
    return pdbmodels
#include FASTA


    #Specifying behavior according to arguments regarding INTPUT:
    #if options.input:
    #    if os.path.isfile(options.input):  #Input is a file
    #        files = [options.input]
    #        path = './'
    #    else:
    #        path = options.input + '/'    #Input is a directory
    #        files = [path + fasta for fasta in os.listdir(path) if fasta[-3:] == '.fa' or fasta[-6:] == '.fasta' or fasta[-6:] == '.fa.gz' or fasta[-9:] == '.fasta.gz']
    #else:   #No specification of input
    #    files = [fasta for fasta in os.listdir('.') if fasta[-3:]=='.fa' or fasta[-6:]=='.fasta']
    #    path = './'
    #Printing number of processed files
#    num_files = len(files)
#    if options.verbose:
#        sys.stderr.write('%s FASTA files found:\n' %(num_files))
