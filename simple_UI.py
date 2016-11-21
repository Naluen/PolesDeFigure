from __future__ import print_function
from lib import figuredePoles as fpf


try:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename
except ImportError:
    print("Python3 Detected...")

try:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
except ImportError:
    print("Python2 Detected...")


Tk().withdraw()
# we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename()
# show an "Open" dialog box and return the path to the selected file
gridder = fpf.main(filename)
