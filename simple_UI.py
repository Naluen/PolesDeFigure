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
filename = askopenfilename(
    title='Choose Poles Figure File...',
    filetypes=[("Raw files","*.raw")]
)
gridder = fpf.main(filename)
