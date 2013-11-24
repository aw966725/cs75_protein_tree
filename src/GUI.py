##
# 
# FILE: GUI.py 
#
# AUTHORS: Alex Welton and Varun Ravishanker
#
# DESCRIPTION: Contains graphical user interface functionality
# for Protein Tree project.
#

from tkinter import *

class PhyloTree(Frame):
    
    def __init__(self, name):
        Frame.__init__(self, name)

        self.name = name

        self.initUI()

    def initUI(self):
        
        self.name.title("Phylogenetic Tree")
        self.pack(fill=BOTH, expand=1)
        
        canvas = Canvas(self)
        canvas.create_oval(10, 10, 80, 80, outline="blue", fill = "red", width=2)

        canvas.pack(fill=BOTH, expand=1)

def main():
    
    root = Tk()
    root.geometry("640x480+300+300")
    tree = PhyloTree(root)
    root.mainloop()


if __name__ == '__main__':
    main()
