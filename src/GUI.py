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
STARTX = 10
STARTY = 50
from tkinter import *
import math

class TreeNode(object):

    def __init__(self, species, name, startx, starty):
        self.species = species
        self.name = name
        self.x = startx
        self.y = starty

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
        self.drawTree(clusters, STARTX, STARTY, 0, 100)
        

    def initUI(self):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)

        self.canvas.pack(fill=BOTH, expand=1)

    def drawTree(self, clusters, x, y, angle, length):
        canvas = self.canvas

        nodes = []
        spacing = CANVAS_HEIGHT / len(self.species_list)
        count = 0
        for species in self.species_list:
            nodes.append(TreeNode([species], species,  x, y + spacing*count))
            count += 1
            
        for node in nodes:
            canvas.create_text(node.x + (MAX_WORD/2), node.y, text=node.species, width=80)
            canvas.create_line(node.x + MAX_WORD, node.y, node.x + MAX_WORD + length, node.y)
            

        for level in range(len(clusters) - 1):
            cur_clusters = clusters[level]
            next_clusters = clusters[level + 1]
            
            for cluster in next_clusters:
                found = False
                for node in nodes:
                    if node.species == cluster:
                        found = True
                        node.x += length
                if not found:
                    rig = rig
                    #nodes.append(TreeNode(cluster, 
            
                
            
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
