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
import math

class PhyloTree(Frame):
    
    def __init__(self, window, clusters):
        Frame.__init__(self, window)

        self.window = window 
        self.canvas = Canvas(self)
        
        self.initUI(clusters)
        self.drawTree(300, 500, math.pi/2, 150, 5)

    def initUI(self, clusters):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)

        self.canvas.pack(fill=BOTH, expand=1)

    def drawTree(self, x, y, angle, length, level):
        canvas = Canvas(self)
        if level == 0:
            self.canvas.create_line(x, y, x + length * math.cos(angle), y - length * math.sin(angle))

        else:
            self.canvas.create_line(x, y, x + length * math.cos(angle), y - length * math.sin(angle))
            self.drawTree(x + length * math.cos(angle), 
                     y - length * math.sin(angle), 
                     angle + math.pi/4, 
                     length * (2.0/3.0), level - 1)

            self.drawTree(x + length * math.cos(angle), 
                     y - length * math.sin(angle), 
                     angle - math.pi/4, 
                     length * (2.0/3.0), level - 1)

def main():
    
    root = Tk()
    root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT)+ "+450+300")
    tree = PhyloTree(root, "Ciona intestinalis")
    root.mainloop()

if __name__ == '__main__':
    main()
