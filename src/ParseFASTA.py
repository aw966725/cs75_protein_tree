##
#
# FILE: ParseFASTA.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains functionality for opening and reading
# a FASTA file, and then parsing the returned data into a useful
# format.
#
# Based on code authored by Sara Thiebaud (code was distributed
# with class materials for COSC 75 Homework 1
#

import os
from collections import defaultdict


## DEFINITIONS

# Maximum FASTA file size in MB that this code is tested to work with
# Allows bigger files, but prints a warning
MAX_FILE_SIZE_MB = 20
MAX_FILE_SIZE_B = MAX_FILE_SIZE_MB * 1048576


## CLASSES

# Class representing a FASTA file (contains reading and parsing methods)
# Based on code by Sara Thiebaud
# Caution: can take a while on large data sets!
class FASTAFile (object):

    # Iterator of FASTARecord objects
    def __iter__(self):
        self.read_and_parse()
        return iter(self.records)

    # New instance - takes filename to open and string name of species
    def __init__(self, filename, species):

        # Validate and open the file
        self.file = self.validate_and_open(filename)
        self.filename = filename
        self.species = species

        # Auto-generates lists for new keys and appends to list for existing keys
        # Keys are gene_IDs, values are lists of FASTARecords with that gene ID
        self.records = defaultdict(list)
        self.read_and_parse()

    # Read in a file and parse data into FASTARecord objects
    def read_and_parse(self):
        # Iterate through lines, parsing fields (they are reproducible widths)
        for line in self.file:

            # Handle trailing newlines and strip before parsing
            if line[-1] == '\n':
                line = line[:-1]
            line = line.strip()

            # Check for beginning of FASTA record (not every new line)!
            # Retrieve all information here
            if line[:1] == '>':
                self.num_records += 1

                # Split the line on spaces into fields for easier parsing
                line_fields = line.split(' ')

                # This is the unique variant ID (can be multiple variants per gene)
                variant_ID = line_fields[0][1:]

                # Get peptide type (choices are known, novel, putative)
                curr_sub_fields = line_fields[1].split(':')
                peptide_type = curr_sub_fields[1]

                # Get chromosome encoding
                curr_sub_fields = line_fields[2].split(':', 1)
                chromosome = curr_sub_fields[1]

                # Get gene ID (multiple FASTARecords may exist for this gene)
                curr_sub_fields = line_fields[3].split(':')
                gene_ID = curr_sub_fields[1]

                # Sequence begins on next line - unknown length so must while loop
                gene_sequence = ""
                seq_line = self.file.readline()
                while seq_line != '' and seq_line[:1] != '>':
                    gene_sequence += seq_line.strip()
                    seq_line = self.file.readline()

                # Create the FASTARecord with all that data
                current_record = FASTARecord(variant_ID, gene_ID, peptide_type, 
                    chromosome, gene_sequence)

                # Store it in the defaultdict under the key gene_ID
                self.records[gene_ID].append(current_record)

            # Not the beginning of a FASTA record - move on
            else:
                line = self.file.readline()

    # Return number of unique records (equal to number of variant IDs)
    def __len__(self):
        return self.num_records

     # Return number of unique gene IDs total (number of variants equal to number of records)
    def get_num_unique_genes(self):
        return len(self.records.keys())

    # Return list of unique gene IDs
    def get_unique_gene_IDs(self):
        return list(self.records.keys())

    # Return list of variants of a given ID
    def get_variants(self, gene_ID):
        return self.records.get(gene_ID)

    # Validate a filename as existent and in the correct format, then open and return it
    def validate_and_open(filename):

        # Is this a string?
        if not isinstance(filename, str):
            raise TypeError("Must pass FASTAFile constructor a string filename.")

        # Is this a .fa file?
        if not filename.lower().endswith('.fa'):
            raise IOError("Must pass FASTAFile constructor a .fa file.")

        # Warn if larger than tested
        file_size = os.path.getsize(filename)
        if file_size > MAX_FILE_SIZE_B:
            print("WARNING: File size of %d bytes is larger than tested max of %d bytes.\n", 
                file_size, MAX_FILE_SIZE_B)

        # Safely try opening it
        try:
            with open(filename, 'r') as FASTA_file:
                return FASTA_file
        except IOError:
            print("Could not open FASTA file.\n")
            return None


# Class representing an individual FASTA record
# Based on code by Sara Thiebaud
class FASTARecord:

    # New instance - takes variant ID, gene ID, peptide type, chromosome encoding, 
    # gene sequence
    def __init__(self, variant_ID, gene_ID, pep_type, chromosome, gene_seq):

        self.variant = variant_ID
        self.gene = gene_ID
        self.peptide_type = pep_type
        self.chromosome = chromosome
        self.sequence = gene_seq
        self.num_records = 0

    # Get length of gene sequence
    def __len__(self):
        return len(self.sequence)

