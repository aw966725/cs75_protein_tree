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


def gui_init():
    width, height = 1024, 768

    root = tkinter.Tk()
    window = tkinter.Toplevel()
    canvas = tkinter.Canvas(window, width=width, height=height)
    canvas.pack()

def get_file_to_open():
    return tkinter.filedialog.askopenfilename()
