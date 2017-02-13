from __future__ import print_function, unicode_literals
import os
import sys
from lib.bruker3 import convert_raw_to_uxd
import codecs
import tempfile
import logging
import struct
from bruker3 import DatasetDiffractPlusV3

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

try:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename
except ImportError:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    print("Python3 Detected...")
else:
    print("Python2 Detected...")

logging.config.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,  # this fixes the problem

    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        },
    },
    'handlers': {
        'default': {
            'filename': os.path.join(
                os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'formatter': 'standard',
        },
    },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'DEBUG',
            'propagate': True
        }
    }
})


class polesFIGURE(object):
    def __init__(self, raw_file):

        self.raw_file = raw_file
        uxd_file = os.path.join(raw_file.replace("raw", "uxd"))

        self.directory = os.path.dirname(self.raw_file)
        config = configparser.ConfigParser()
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(self.directory, 'config.ini')
        ])
        self.preference = dict(config._sections['db'])

    @staticmethod
    def raw_file_reader(input_raw_file):
        ds = DatasetDiffractPlusV3(open(input_raw_file, 'rb'))
        data = (ds.pretty_format(print_header=True)).encode('ISO-8859-1')
        return data


    def _read_binary_data(self):


    def plot_2d_image(self, is_show_image=False):
        """Plot the 2D Image"""


if __name__ == '__main__':
    Tk().withdraw()
    raw_file = askopenfilename(
        title='Choose Poles Figure File...',
        filetypes=[("Raw files", "*.raw")]
    )

    sample = polesFIGURE(raw_file)
