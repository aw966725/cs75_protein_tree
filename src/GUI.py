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
MAX_WORD = 20

from tkinter import *
import math

class PhyloTree(Frame):
    
    def __init__(self, window, clusters):
        Frame.__init__(self, window)

        self.window = window 
        self.canvas = Canvas(self)
        
        self.initUI(clusters)
        self.drawTree(0, WINDOW_HEIGHT/2, 0, 150, 4)

    def initUI(self, clusters):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)

        self.canvas.pack(fill=BOTH, expand=1)

    def drawTree(self, x, y, angle, length, level):
        canvas = Canvas(self)
        if level == 0:
            self.canvas.create_line(x, y, x + length, y)
            self.canvas.create_text((x + length + MAX_WORD), y, text="Rig", width=80)

        else:
            
            #Change length for next iteration based on relationship strength
            self.canvas.create_line(x, y, x + length, y)
            self.canvas.create_line(x + length, y, x + length, y + length)
            self.canvas.create_line(x + length, y, x + length, y - length)
            self.drawTree(x + length, 
                     y - length, 
                     angle + math.pi/4, 
                     length * (1.0/2.0), level - 1)

            self.drawTree(x + length, 
                     y + length, 
                     angle - math.pi/4, 
                     length * (1.0/2.0), level - 1)

def main():
    
    root = Tk()
    root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT)+ "+450+300")
    tree = PhyloTree(root, "Ciona intestinalis")
    root.mainloop()

if __name__ == '__main__':
    main()
