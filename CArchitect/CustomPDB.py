#!/usr/bin/env python
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio import pairwise2
import utilities

class CustomModel(Model):
    """Custom biopython model class that allows models with more than one same id"""
    def add(self, entity):
        """Add a child to the Entity."""
        entity.set_parent(self)
        self.child_list.append(entity)


class CustomChain(Chain):      # dictionaries in utilities
    """Description"""
    protein_dict = utilities.protein_dict
    dna_dict = utilities.dna_dict
    rna_dict = utilities.rna_dict

    def __init__(self, chainObject):
        self.child_list = chainObject.child_list
        self._id = chainObject._id
        self.parent = chainObject.parent
        self.interactions = []  # Here there will be stored the chain's known interactions as tuple of residues
        self.xtra = chainObject.xtra
        self.level = chainObject.level


    def get_sequence(self):
        """Returns the chain's sequence, it can be a protein, DNA or RNA sequence"""
        sequence = ""
        first_residue = self.child_list[0].resname.strip()

        if first_residue not in self.protein_dict:
            if "D" in first_residue:
                for res in self:
                    if res.id[0] == " ":
                        sequence += self.dna_dict[res.resname.strip()]
            else:
                for res in self:
                    if res.id[0] == " ":
                        sequence += self.rna_dict[res.resname.strip()]
        else:
            for res in self:
                if res.id[0] == " ":
                    sequence += self.protein_dict[res.resname]
        return sequence


    def has_equivalent(self, other_sequence):
        cutoff = 0.95
        sequence1 = self.get_sequence()
        sequence2 = other_sequence

        alignment = pairwise2.align.globalxx(sequence1, sequence2)[0]
        align_score = alignment[2]
        align_length = (len(alignment[0]))

        cutoff_chain = align_score / align_length
        if cutoff_chain > cutoff:
            return True
        else:
            return False
