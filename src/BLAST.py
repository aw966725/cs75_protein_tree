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






# TODO: add minimum sample to functions






from imports import *
import math
from collections import defaultdict
from random import *

import Matching
import NormalLookup

## DEFINITIONS

# Default word length "w"
DEFAULT_WORD_LEN = 3
DEFAULT_SUB_MATRIX = SubstitutionMatrix.blosum62

# Score threshold percentage - i.e. match must be in the X highest percent of scores
# evaluated thus far or it is discarded
DEFAULT_THRESHOLD = 5

# Minimum sample - the number of samples at which we begin to test the threshold
MINIMUM_SAMPLE = 10

## CLASSES

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

        # Dynamically adjust matching threshold
        self.total_match_score = 0
        self.total_match_squared_score = 0
        self.total_number_matches = 0
        self.threshold_score = 0

        # Dynamically adjust alignment threshold
        self.total_align_score = 0
        self.total_align_squared_score = 0
        self.total_number_aligns = 0
        self.align_threshold_score = 0

        # Dynamically adjust relation threshold
        self.total_relation_score = 0
        self.total_relation_squared_score = 0
        self.total_number_relations = 0
        self.relation_threshold_score = 0

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
        z = lookup_z_score(100 - self.threshold)

        # Get approximate threshold score
        self.threshold_score = (z * std_dev) + mean
        return

    # Recalculate the alignment threshold (same as above)
    def recalculate_alignment_threshold_score(self):

        # Get Z-score
        mean = total_align_score / total_number_aligns
        std_dev = math.pow(((pow(self.total_align_squared_score), 2) / self.total_number_aligns)
            - pow(mean, 2))
        z = lookup_z_score(100 - self.threshold)

        # Get approximate threshold score
        self.align_threshold_score = (z * std_dev) + mean
        return

    # Recalculate the relation threshold (same as above)
    def recalculate_relation_threshold_score(self):

        # Get Z-score
        mean = total_relation_score / total_number_relations
        std_dev = math.pow(((pow(self.total_relation_squared_score), 2) / self.total_number_relations)
            - pow(mean, 2))
        z = lookup_z_score(100 - self.threshold)

        # Get approximate threshold score
        self.align_threshold_score = (z * std_dev) + mean

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

        # Relations defined between genes as families        
        self.relations = []

        # BLASTVariantPair objects
        self.variant_pairs = self.get_alignments()

        # "Families" of related proteins
        self.protein_families = self.get_protein_families()

        # Info struct
        if info != None and isinstance(info, BLASTInfo):
            self.info = info

    # Process alignments for the two species (may take a long time)
    def get_alignments(self):

        variant_pairs = []

        # Iterate over first species' variants
        for gene_a in self.species_a.genes:
            for variant_a in gene_a:

                # Iterate over second species' variants
                for gene_b in self.species_b.genes:
                    for variant_b in gene_b:

                        # Create a BLASTVariantPair object and get the alignments
                        this_pair = BLASTVariantPair(variant_a, variant_b, self.info)
                        variant_pairs.append(this_pair)

                        # Test for relation between these variants
                        #if
        return

    # Categorize the proteins into "families" to decrease processing time
    def get_protein_families(self, num_families=3):

        clusters = {}

        # Randomly select initial variants for clustering
        for i in range(num_families):
            initCentroid = variant_pairs[random.randint(0, len(variant_pairs))].variant_a
            clusters[initCentroid] = []

        prev_clusters = {}

        # Keep iterating until clusters don't change
        while True:
                
            # Assign variants to closest clusters
            for gene in self.species_a.genes:
                for variant in gene.variants:
                    clusters = assign_to_cluster(variant, clusters)
                    
            for gene in self.species_b.genes:
                for variant in gene.variants:
                     clusters = assign_to_cluster(variant, clusters)

            #TODO: This might not test for equality properly, have to check
            if clusters == prev_clusters:
                break
                
            prev_clusters = clusters
            clusters = recompute_centroids(prev_clusters)

        return clusters
            
            
                    
    # Assign a given variant to the closest cluster
    def assign_to_cluster(self, variant, clusters):
        closest_centroid = None
        score = 0
        
        #skip variant if it's one of the centroids
        if variant in clusters:
            return clusters
        
        # find all sets of pairs with the current variant present and
        # keep track of the highest scoring centroid
        for variant_pair in self.variant_pairs:
            if variant == variant_pair.variant_a and variant_pair.variant_b in clusters:
                if variant_pair.score > score:
                    closest_centroid = variant_pair.variant_b
                    score = variant_pair.score 
                elif variant == variant_pair.variant_b and variant_pair.variant_a in clusters:
                    if variant_pair.score > score:
                        closest_centroid = variant_pair.variant_a
                        score = variant_pair.score
                        

                        # add variant to highest scoring centroid's cluster
                        if closest_centroid == None:
                            print("Something's wrong")
                        else:
                            clusters[closest_centroid].append(variant)
        
        return clusters

    # ASSUMES THERE WILL BE VARIANT PAIRS FOR EVERY SET OF VARIANTS
    def recompute_centroids(self, clusters):
        new_clusters = {}
        
        for centroid in clusters:
            avg_scores = {}
            
            # Get list of variants in the cluster
            variants = clusters[centroid]
            variants.append(centroid)

            # For each variant, compute avg distance between other variants
            for variant1 in variants:
                avg_scores[variant1] = 0
                for variant2 in variants:
                    if variant1 == variant2:
                        continue

                    # Find variant_pair of variant1corresponding to each variant2 to compute average score
                    for variant_pair in self.variant_pairs:
                        if (variant_pair.variant_a == variant1 and variant_pair.variant_b == variant2) or (variant_pair.variant_b == variant1 and variant_pair.variant_a == variant2):
                            avg_scores[variant1] += variant_pair.score
                avg_scores[variant1] /= len(variants)

            new_centroid = max(avg_scores, key=avg_scores.get)
            new_clusters[new_centroid] = []

        return new_clusters

            

    # Score species pair by blasting each variant against each other
    def score_species_pair(self):
        return


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

        # BLAST these protein sequences against each other and obtain score
        self.score = score_variant_pair()

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

    # Get the best alignment for a query and a similar word
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
            while self.info.matrix[self.variant_a.sequence[i+align_length]][self.variant_b.sequence[j+align_length]]:
                align_length += 1
                score += self.info.matrix[self.variant_a.sequence[i+align_length]][self.variant_b.sequence[j+align_length]]

        # Return substrings and score
        return (self.variant_a.sequence[i:i+align_length],
            self.variant_b.sequence[j:j+align_length], score)

    # BLAST alignment function - takes two sequences
    def align_BLAST(self, seq_a, seq_b):

        # Get threshold-tested words
        threshold_words = get_threshold_words()
        for (query, word_list) in threshold_words.items():
            for word in word_list:
                best = get_best_alignment_BLAST(query, word)
                if best[2] > self.info.align_threshold_score:
                    # Alignment dictionary keys are tuples of the two variant IDs
                    self.info.alignments[(self.variant_a.variant_ID, self.variant_b.variant_ID)].append(best)

        # Recalculate score threshold
        self.info.recalculate_alignment_threshold_score()
        return

    # Scores an alignment of a pair of variants by blasting them and summing acceptable alignments
    def score_variant_pair(self):
        align_BLAST(self.variant_a.sequence, self.variant_b.sequence)
        return sum(self.info.alignments[(self.variant_a.variant_ID, self.variant_b.variant_ID)])
