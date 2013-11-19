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


# Class representing an individual FASTA record
# Based on code by Sara Thiebaud
class FASTARecord:

    # New instance - takes variant ID, gene iD, chromosome type, gene sequence
    def __init__(self, var_ID, gene_ID, chrom_type, gene_seq):
        self.variant = var_ID
        self.gene = gene_ID
        self.chromosome_type = chrom_type
        self.sequence = gene_seq

    # Get length of gene sequence
    def __len__(self):
        return len(self.sequence)


# Class representing a FASTA file (contains reading and parsing methods)
# Based on code by Sara Thiebaud
# Caution: can take a while on large data sets!
class FASTAFile (object):

    # Iterator of FASTARecord objects
    def __iter__(self):
        self.read_and_parse()
        return iter(self.records)

    # New instance - takes filename to open
    def __init__(self, filename):

        # Validation
        if not isinstance(filename, str):
            raise TypeError("Must pass FASTAFile constructor a filename.")

        # Initialize storage and read in the file
        self.records = []
        self.filename = filename
        self.read_and_parse()

    # Read in a file and parse data into FASTARecord objects
    def read_and_parse(self):
        self.file = open(self.filename)

        # Iterate through lines, parsing fields (they are reproducible widths)
        # TODO: ADD ERROR CHECKING FOR FILE FORMAT IN CONSTRUCTOR
        for line in self.file:

            # Handle trailing newlines and strip
            if line[-1] == '\n':
                line = line[:-1]
            line = line.strip()


