##
# 
# FILE: GUI.py 
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains graphical user interface functionality
# for Protein Tree project.
#

WINDOW_WIDTH = 640
WINDOW_HEIGHT = 640
MAX_WORD = 80
CANVAS_HEIGHT = 640
CANVAS_WIDTH = 640
from tkinter import *
import math

class TreeDraw(Frame):
    
    def __init__(self, window, clusters):
        Frame.__init__(self, window)

        self.clusters = clusters
        self.window = window 
        self.canvas = Canvas(self, width=CANVAS_WIDTH, height=CANVAS_HEIGHT)

        # Get list of species names
        species_list = []
        for cluster in clusters[0]:
            species_list.append(cluster[0])

        self.species_list = species_list
        
        self.initUI()
        self.drawTree(clusters, 10, CANVAS_HEIGHT/2, 0, 100)
        

    def initUI(self):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)

        self.canvas.pack(fill=BOTH, expand=1)

    def drawTree(self, clusters, x, y, angle, length):
        canvas = self.canvas
        curx = x
        deltay = length
        
        # Draw the clustering backwards
        for level in range(len(clusters) - 1, 0, -1):
            cluster_list = clusters[level]
            next_clusters = clusters[level - 1]
        
            flat_clusters = [species for cluster in cluster_list for species in cluster]
            flat_next = [species for cluster in next_clusters for species in cluster]
            
            # list of old clusters that remained as is
            pair1 = list(set(flat_next) & set(flat_clusters))
            
            # list of new clusters
            pair2 = list(set(flat_next) - set(flat_clusters))

            same_clusters = []
            changed_clusters = []

            #generate subclusters based on presence in current cluster
            for cluster in next_clusters:
                if cluster[0] in pair1:
                    same_clusters.append(cluster)
                elif cluster[0] in pair2:
                    changed_clusters.append(cluster)
                else:
                    print("uh oh")
            
            # reorder species list
            self.species_list = same_clusters + changed_clusters

            # draw current level
            canvas.create_line(curx, y, curx, y + deltay)
            canvas.create_line(curx, y, curx, y - deltay)
            canvas.create_line(curx, y + deltay, curx + length, y + deltay)
            canvas.create_line(curx, y - deltay, curx + length, y - deltay)
            curx += length
            deltay += length

        # level now zero, so finish up
        # fix this by drawing same as above
        for i in range(len(self.species_list)):
            canvas.create_line(curx, y + deltay, curx + length, y + deltay)
            canvas.create_text(curx + length + MAX_WORD, y - deltay, text=self.species_list[i], width=80)
            deltay -= length*2
            
        
        
        
        
        
            
def main():
    
    root = Tk()
    root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT)+ "+450+300")

    clusters = [[] for x in range(3)]
    
    clusters[0].append(["Dog"])
    clusters[0].append(["Cat"])
    clusters[0].append(["Snake"])
    clusters[1].append(["Dog", "Cat"])
    clusters[1].append(["Snake"])
    clusters[2].append(["Dog", "Cat", "Snake"])
            
         
    tree = TreeDraw(root, clusters)
    root.mainloop()

if __name__ == '__main__':
    main()
