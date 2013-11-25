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

class PhyloTree(Frame):
    
    def __init__(self, window, clusters):
        Frame.__init__(self, window)

        self.clusters = clusters
        self.window = window 
        self.canvas = Canvas(self, width=CANVAS_WIDTH, height=CANVAS_HEIGHT)
        
        self.initUI()
        self.drawTree(clusters, 200)
        

    def initUI(self):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)

        self.canvas.pack(fill=BOTH, expand=1)

    def drawTree(self, clusters, length):
        canvas = self.canvas
        startx = 10
        cury = CANVAS_HEIGHT/2

        
        
        
            
def main():
    
    root = Tk()
    root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT)+ "+450+300")

    clusters = [[] for x in range(2)]
    
    clusters[0].append(["Dog"])
    clusters[0].append(["Cat"])
            
         
    tree = PhyloTree(root, clusters)
    root.mainloop()

if __name__ == '__main__':
    main()
