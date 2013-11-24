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
WINDOW_HEIGHT = 480

from tkinter import *

class PhyloTree(Frame):
    
    def __init__(self, window, clusters):
        Frame.__init__(self, window)

        self.window = window 

        self.initUI(clusters)

    def initUI(self, clusters):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)
        
        canvas = Canvas(self)
        startx = 50
        starty = WINDOW_HEIGHT / 2
        curx = startx
        ydist = 0
        colsize = 180
        rowsize = 100

        #determine this later
        name_width = 70
        
        #width should be text size divided by depth of clusters
        canvas.create_text(startx, starty, text=clusters, width=name_width)

        #width should be dependent on strength of relationship
        curx = startx + colsize
        canvas.create_line(startx + name_width, starty, curx, starty, width=2)

        #compute ydist properly
        canvas.create_line(curx, starty, curx, starty + rowsize, width=2)
        canvas.create_line(curx, starty, curx, starty - rowsize, width=2)        
    
        
        
        
        canvas.pack(fill=BOTH, expand=1)

def main():
    
    root = Tk()
    root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT)+ "+450+300")
    tree = PhyloTree(root, "Ciona intestinalis")
    root.mainloop()

if __name__ == '__main__':
    main()
