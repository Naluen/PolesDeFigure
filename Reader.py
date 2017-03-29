from H5Reader import H5Reader as H5Reader
from RawReader import RawReader as RawReader
import os


class Reader(object):
    def __init__(self, file):
        self.file = file

    def matrix_data(self):
        file_type = os.path.basename(self.file).split('.')[1]
        file_type_dict = {
            'h5': H5Reader,
            'raw': RawReader
        }
        data_dict, instance_dict = file_type_dict[file_type](
            self.file).matrix_data()
        return data_dict, instance_dict
