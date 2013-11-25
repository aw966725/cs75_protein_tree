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
        self.clusters = compute_clusters()


    # computes species clusters based on profiles
    # Heirarchical?
    def compute_clusters(self):
        profile_list = self.profiles.profile_list

        # how do you determine distance function once they're clustered?
        # average scores of checks against each in cluster?
        clusters = []

        # level to assess depth of cluster formation
        level = 0
        for profile in profile_list:
            cluster = { count : ([profile], level)}
            clusters.append(cluster)

        #compute min distance and cluster
        (x, y) = compute_min_pair(clusters)


    def compute_min_pair(clusters):
        for cluster in clusters:
            p

        
        


class Profiles(object):

    # New instance - takes a list of species and a dictionary of protein families 
    def __init__(self, species_list, protein_families):

        self.species_list = species_list
        self.protein_families = protein_families
        self.profile_list = create_profiles()

    def create_profiles(self):

        profiles = {}
        
        for species in self.species_list:
            profile_list[species] = new Profile(species, protein_families)
                
        return profiles
            
class Profile(object):

    # New instance - takes a species and a dictionary of protein families
    def __init___(self, species, protein_families):

        self.species = species
        self.species_proteins = species.genes.variants
        self.protein_families = protein_families
        self.profile = initialize_profile(self)
        compute_profile()

    # Initialize profile to all zeroes (not present)
    def initialize_profile(self):
        profile = {}

        # Initialize all protein families to not present
        for family in self.protein_families:

            # Are we doing partial points for variants now that we're
            # doing families?
            profile[family] = 0

        return profile

    # Iterate over each family and see if any family members
    # are present in the species's protein list. 
    def compute_profile(self):
        for family in self.protein_families:
            members = self.protein_families[family]

            for member in members:
                if member in self.species_proteins:
                    self.profile[family] = 1
                    break


    
        
