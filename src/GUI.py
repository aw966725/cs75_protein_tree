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
CANVAS_HEIGHT = 500
CANVAS_WIDTH = 500
from tkinter import *
import math

class PhyloTree(Frame):
    
    def __init__(self, window):
        Frame.__init__(self, window)

        self.window = window 
        self.canvas = Canvas(self, width=CANVAS_WIDTH, height=CANVAS_HEIGHT)
        
        self.initUI()
        

    def initUI(self):
        
        self.window.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)

        #uses initial values
        level_slider = Scale(self.window, from_=1, to=5, command = self.updateTree, length=50)
        level_slider.pack(fill=BOTH, side=RIGHT)

        self.canvas.pack(fill=BOTH, expand=1)

    def drawTree(self, x, y, angle, length, level):
        canvas = self.canvas
        if level == 0:
            line = self.canvas.create_line(x, y, x + length, y)
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

    def updateTree(self, new_level):
        self.canvas.delete(ALL)
        self.drawTree(0, WINDOW_HEIGHT/2, 0, 150, int(new_level))
            
def main():
    
    root = Tk()
    root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT)+ "+450+300")
    tree = PhyloTree(root)
    root.mainloop()

if __name__ == '__main__':
    main()
