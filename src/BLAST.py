##
#
# FILE: BLAST.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains functionality for implementing
# gapped BLAST subsequence alignment using the blosum62
# matrix.
#

from imports import *
from math import *
from collections import defaultdict
import Matching


# Default word length "w"
DEFAULT_WORD_LEN = 3
DEFAULT_SUB_MATRIX = blosum62

# Score threshold percentage - i.e. match must be in the X highest percent of scores
# evaluated thus far or it is discarded
DEFAULT_THRESHOLD = 5

# Class for holding general information about the application of a BLAST algorithm
class BLASTInfo (object):

    # Initialize - takes a list of FASTAFile objects, word length, substitution matrix, and threshold
    def __init__(self, matrix=DEFAULT_SUB_MATRIX, word_length=DEFAULT_WORD_LEN, 
        threshold=DEFAULT_THRESHOLD):

        # Info about BLAST
        self.matrix = matrix
        self.word_length = word_length
        self.threshold = threshold
        self.all_words = self.get_all_possible_words(matrix, word_length)

        # To dynamically adjust score threshold
        self.total_match_score = 0
        self.total_match_squared_score = 0
        self.total_number_matches = 0
        self.threshold_score = 0

        # Dynamically adjust matching threshold
        self.total_align_score = 0
        self.total_align_squared_score = 0
        self.total_number_aligns = 0
        self.align_threshold = 0

        self.alignments = defaultdict(list)

    # Here we get a list of all possible word_length length words in matrix and store it
    def get_all_words(self, matrix, word_length):
        all_words = []
        count = 0
        return self.get_all_words_helper(all_words, matrix, count, word_length)

    # Recursive helper method for getting all words
    # Takes a word list, sub matrix, char count so far, word length
    def get_all_words_helper(self, word_list, matrix, count, word_length):
        # Base case
        if count == word_length:
            return word_list

        # Recursive case
        words_so_far = len(word_list)
        for i in range(words_so_far):
            # Add all new possible permutations
            for j in range(len(matrix)):
                new_word = word_list[i] + matrix[j]
                word_list.append(new_word)
            # Remove old base word
            word_list.pop(i)

        # Call again with new word list and longer length
        return self.get_all_words_helper(word_list, matrix, count + 1, word_length)

    # Re-calculate the approximated threshold score after each new matching
    def recalculate_threshold_score(self):

        # Get Z-score
        mean = total_match_score / total_number_matches
        std_dev = math.pow(((pow(self.total_match_squared_score), 2) / self.total_number_matches)
            - pow(mean, 2))

        # Get approximate threshold score
        self.threshold_score = (((100 - self.threshold) / 100) * self.total_number_matches) + 0.5
        return

    # Recalculate the alignment threshold (same as above)
    def recalculate_alignment_threshold(self):

        # Get Z-score
        mean = total_align_score / total_number_aligns
        std_dev = math.pow(((pow(self.total_align_squared_score), 2) / self.total_number_aligns)
            - pow(mean, 2))

        self.align_threshold = (((100 - self.threshold) / 100) * self.total_number_aligns) + 0.5
        return

    # Score two words of length "word_length" using the given substitution matrix
    def BLAST_score(self, seq_a, seq_b, matrix, word_length):
        marg_score = 0
        for i in range(word_length):
            marg_score += matrix[seq_a[i]][seq_b[i]]

        score += marg_score
        self.total_match_score += marg_score
        self.total_match_squared_score += pow(marg_score, 2)
        self.total_number_matches += 1

        return score


# Class for implementing BLAST alignment between two Species objects
class BLASTSpeciesPair (object):

    # Initialize - takes two Species objects and an info class object
    def __init__(self, species_a, species_b, info=None):

        self.species_a = species_a
        self.species_b = species_b

        # Info struct
        if info != None and isinstance(info, BLASTInfo):
            self.info = info


# Class for containing BLAST alignment data for a given pair of Variants
class BLASTVariantPair (object):

    # Initialize - takes two Variant objects and an info class object and
    # scores their alignment against one another
    def __init__(self, variant_a, variant_b, info=None):

        self.variant_a = variant_a
        self.variant_b = variant_b

        # Info struct
        if info != None and isinstance(info, BLASTInfo):
            self.info = info

    # Get the threshold-tested words in each of the subsequences
    def get_threshold_words(self):

        # Iterate over and check all word combinations, keeping only those above threshold
        threshold_words = defaultdict(list)
        for i in range(len(self.variant_a.sequence) - self.info.word_length + 1):
            for j in range(len(self.variant_b.sequence) - self.info.word_length + 1):
                if BLAST_score(self.variant_a.sequence[i:i+self.info.word_length],
                    self.variant_b.sequence[j:j+word_length]) > self.info.threshold_score:
                
                    threshold_words[self.variant_a.sequence[i:i+self.info.word_length]].append(self.variant_b.sequence[j:j+self.infoword_length])

        # Recalculate threshold
        self.info.recalculate_threshold_score()

        return threshold_words

    # Get the best alignment for seq_a and seq_b
    def get_best_alignment_BLAST(self, query, word):

        # Get initial indices in full sequence and initial score
        i = self.variant_a.sequence.find(query)
        j = self.variant_b.sequence.find(word)
        align_length = 0
        score = BLAST_score(query, word, self.info.matrix, self.info.word_length)

        while i > 0 and j > 0:
            while self.info.matrix[self.variant_a.sequence[i-1]][self.variant_b.sequence[j-1]] > 0:
                i -= 1
                j -= 1
                align_length += 1
                score += self.info.matrix[self.variant_a.sequence[i-1]][self.variant_b.sequence[j-1]]

        while i <= len(self.variant_a.sequence) and j <= len(self.variant_b.sequence):
            while self.info.matrix[self.variant_a.sequence[i+self.info.word_length]][self.variant_b.sequence[j+self.info.word_length]]:
                align_length += 1
                score += self.info.matrix[self.variant_a.sequence[i+self.info.word_length]][self.variant_b.sequence[j+self.info.word_length]]

        # Return substrings and score
        return (self.variant_a.sequence[i:i+self.info.word_length],
            self.variant_b.sequence[j:j+self.info.word_length], score)

    # BLAST alignment function - takes two sequences and an info object
    def align_BLAST(self, seq_a, seq_b, info):

        # Get threshold-tested words
        threshold_words = get_threshold_words()
        for (query, word_list) in threshold_words.items():
            for word in word_list:
                best = get_best_alignment_BLAST(query, word)
                if best[2] > self.info.align_threshold:
                    self.info.alignments.append(best)

        # Recalculate score threshold
        self.info.recalculate_alignment_threshold()

        return
