##
#
# FILE: Profiling.py
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains functionality for implementing
# phylogentic profiling of species based on protein families
#
#

class SpeciesClusters(object):

    # New instance - takes a list of species and a dictionary of protein families
    def __init__(self, species_list, protein_families):
        self.species_list = species_list
        self.protein_families = protein_families
        self.profiles = Profiles(species_list, protein_families)
        self.clusters = self.compute_clusters()


    # Computes species clusters based on profiles
    def compute_clusters(self):
        profile_list = self.profiles.profile_list
        clusters = []
        
        # Level to assess depth of cluster formation
        level = 0
        for species in profile_list:
            clusters.append([])
            cluster = [species]
            clusters[level].append(cluster)
        
        while len(clusters[level]) > 1:
            # go to next level
            clusters.append([])
            level += 1
            
            #compute min distance and cluster
            (cluster1, cluster2) = self.compute_most_similar_pair(clusters)
            
            clusters[level] = []
            for cluster in clusters:
                
                # Merge cluster1 and cluster2 into a new cluster at cur level
                if cluster == cluster1:
                    #copy cluster into new cluster
                    new_cluster = list(cluster)
                    new_cluster.append(cluster2)
                    clusters[level].append(new_cluster)
                    
                # Ignore second cluster because merged with first cluster at cur level
                elif cluster == cluster2:
                    continue
                    
                # If not part of group, just readd as is to next level
                else:
                    clusters[level].append(cluster)

        return clusters
        

    # Computes the most similar cluster pair
    def compute_most_similar_pair(self,clusters):
        max_pair = (None, None)
        max_similarity = 0

        # Check each set of profiles against all others
        for test_cluster in clusters:

            for other_cluster in clusters:
                if test_cluster == other_cluster:
                    continue
                
                cur_best = self.compute_cluster_similarity(test_cluster, other_cluster)

                if cur_best > max_similarity:
                    max_similarity = curbest
                    max_pair = (test_cluster, other_cluster)

        return max_pair
            

    # Currently rudimentary: minimum distance between any two profiles
    # Computes the similarity measure between two profile sets corresponding to clusters
    def compute_cluster_similarity(self, cluster1, cluster2):
        max_similarity = 0
        
        species_list = fold(cluster1)
        

        for species in cluster1:
            for others in cluster2:
                        cur_similarity = self.compute_profile_similarity(species, other)
                        if cur_similarity > max_similarity:
                            max_similarity = cur_similarity

        return max_similarity

    # Computes similarity measured by summing presence of proteins in common families
    def compute_profile_similarity(self, species1, species2):
        v1 = self.profiles.profile_list[species1].profile_vector
        v1 = self.profiles.profile_list[species1].profile_vector
        count = 0
        
        for family in self.protein_families:
            if v1[family] == 1 and v2[family] == 1:
                count+=1
                
        return count
        


class Profiles(object):

    # New instance - takes a list of species and a dictionary of protein families 
    def __init__(self, species_list, protein_families):

        self.species_list = species_list
        self.protein_families = protein_families
        self.profile_list = self.create_profiles()

    def create_profiles(self):

        profiles = {}
        
        for species in self.species_list:
            profiles[species] = Profile(species, self.protein_families)
                
        return profiles
            
class Profile(object):

    # New instance - takes a species and a dictionary of protein families
    def __init__(self, species, protein_families):

        self.species = species

        genes = species.genes
        species_proteins = []
        for gene in genes:
            species_proteins.append(gene.variants)
            
        self.species_proteins = species_proteins
        self.protein_families = protein_families
        self.profile_vector = self.initialize_profile()
        self.compute_profile()

    # Initialize profile to all zeroes (not present)
    def initialize_profile(self):
        profile_vector = {}

        # Initialize all protein families to not present
        for family in self.protein_families:

            profile_vector[family] = 0

        return profile_vector

    # Iterate over each family and see if any family members
    # are present in the species's protein list. 
    def compute_profile(self):
        for family in self.protein_families:
            members = self.protein_families[family]

            for member in members:
                if member in self.species_proteins:
                    self.profile_vector[family] = 1
                    break
