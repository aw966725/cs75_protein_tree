##
#
# FILE: Matching.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains data structures for representing individual
# peptide genes and their variants, as well as methods to determine two gene's
# relative distane apart.
#
# Blosum62 alphabet matrix code from COSC 75 Homework 1
#

import ParseFASTA
from imports import *


# Represents a phylogenetic profile tree of multiple profiles
class PhyloTree(object):

    # New instance - takes list of FASTAFile objects
    def __init__(self, FASTA_list):

        # Information about self
        self.species_names = [FASTA_file.species for FASTA_file in FASTA_list]
        self.num_species = len(FASTA_list)

        # Create list of Species objects for child classes
        self.species = self.species_list_FASTA(FASTA_list)

    # Takes a list of FASTAFile objects and returns a list of Species objects
    def species_list_FASTA(self, FASTA_list):
        temp_list = []
        for FASTA_file in FASTA_list:
            temp_list.append(Species(FASTA_file, self))
        return temp_list

    # Return number of species in tree
    def __len__(self):
        return self.num_species


# Represents a species' phylogenetic profile
class Species (object):

    # New instance - takes a FASTAFile object
    def __init__(self, FASTA_file, phylo_tree=None):

        # Information about self
        self.species = FASTA_file.species
        self.num_genes = FASTA_file.get_num_unique_genes()

        # Container class
        if phylo_tree != None and isinstance(phylo_tree, PhyloTree):
            self.phylotree = phylo_tree

        # Create list of Gene objects for child classes
        self.genes = self.genes_list_FASTA(FASTA_file)

    # Takes a FASTA file and returns a list of Gene objects
    def genes_list_FASTA(self, FASTA_file):
        temp_list = []
        for (gene_ID, variant_list) in FASTA_file.records.items():
            temp_list.append(Gene(variant_list, self))
        return temp_list

# Represents a peptide gene
class Gene (object):

    # New instance - takes a list of FASTARecord objects
    def __init__(self, variants, species=None):

        # Information about gene overall
        self.gene_ID = variants[0]
        self.num_variants = len(variants)
        self.variant_IDs = [var.variant_ID for var in variants]

        # Container class
        if species != None and isinstance(species, Species):
            self.species = species

        # Information about individual gene variants
        self.variants = self.variants_from_records(variants)
    
    # Take our FASTARecord structs and move that info to self
    def variants_from_records(self, variant_list):
        temp_list = []
        for variant in variant_list:
            temp_list.append(Variant(variant))
        return temp_list

    # Returns number of variants of this gene
    def __len__(self):
        return self.num_variants


# Represents a unique gene variant - inherits Gene superclass. Should
# only be created by internal Gene methods
# Subclass of Gene
class Variant (object):

    # New instance - takes a FASTARecord object
    def __init__(self, record, gene=None):

        # Basic information about variant
        self.variant_ID = record.variant_ID
        self.peptide_type = record.peptide_type
        self.chromosome = record.chromosome
        self.sequence = record.sequence
        self.gene_ID = record.gene_ID

        # Container class
        if gene != None and isinstance(gene, Gene):
            self.gene = gene

    # Return length of sequence
    def __len__(self):
        return len(self.sequence)


