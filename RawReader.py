import configparser
import logging
import os
import re

import numpy as np

from bruker3 import DatasetDiffractPlusV3


class DetectorScanRaw(object):
    def __init__(self, data_list):
        self.data_list = data_list

    def process(self):
        int_data_list = [float(i.split()[1]) for i in self.data_list[1]
                         if (re.match('\d', i) or i.startswith("-"))]
        int_data_matrix = np.asanyarray(int_data_list)
        return {'int_data': int_data_matrix}


class TwoDScanRaw(DetectorScanRaw):
    def process(self):
        import re
        phi_data_list = [float(i.split()[0]) for i in self.data_list[1]
                         if (re.match('\d', i) or i.startswith("-"))]
        int_data_list = [
            [float(j.split('\t')[1]) for j in value
             if (re.match('\d', j) or j.startswith("-"))]
            for value in self.data_list[1:]
            ]
        step_time_list = [
            [float(j.split('=')[1]) for j in value if
             j.startswith('_STEPTIME')]
            for value in self.data_list[1:]
            ]
        khi_data_list = [
            [float(j.split('=')[1]) for j in value if j.startswith('_KHI')]
            for value in self.data_list[1:]
            ]
        phi_data_matrix = np.asanyarray(phi_data_list)
        int_data_matrix = (
            np.asanyarray(int_data_list) / np.asanyarray(step_time_list)
        )
        khi_data_matrix = np.asanyarray(khi_data_list)
        return {
            'int_data': int_data_matrix,
            'phi_data': phi_data_matrix,
            'khi_data': khi_data_matrix
        }


class RawReader(object):
    def __init__(self, raw_file):
        logging.info("Reading raw file {0}".format(raw_file))
        self.raw_file = raw_file

    @staticmethod
    def get_head(data_list):
        scan_dict = {i.split('=')[0].strip(): i.split('=')[1].strip()
                     for i in data_list[0] if i.startswith("_")}
        return scan_dict

    def str_data(self):
        ds = DatasetDiffractPlusV3(open(self.raw_file, 'rb'))
        data_string = ds.pretty_format(print_header=True)
        data_list = data_string.split('\n')
        from itertools import groupby
        data_list = [list(g) for k, g in
                     groupby((line.strip() for line in data_list), bool) if k]

        return data_list

    @staticmethod
    def read_config_file(raw_file):
        def guess_sample_name(raw_file_path):
            try:
                sample_name = (
                    re.findall(
                        r'S\d\d\d\d',
                        raw_file_path,
                        flags=re.IGNORECASE
                    )
                )
            except FileNotFoundError:
                print("Error, no sample name.")
                return ''
            else:
                if sample_name:
                    sample_name = sample_name[-1]
                return sample_name

        config = configparser.ConfigParser()
        import sys
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(os.path.dirname(raw_file), 'config.ini')
        ])
        config['db']['sample'] = guess_sample_name(raw_file)
        config['db']['directory'] = os.path.dirname(raw_file)
        config['db']['raw_file'] = os.path.basename(raw_file)

        return config

    def matrix_data(self):
        data_list = self.str_data()
        scan_dict = self.get_head(data_list)
        process_dict = {
            'SingleScanPlot': DetectorScanRaw,
            'TwoDPlot': TwoDScanRaw
        }
        data_dict = process_dict[scan_dict['_TYPE']](data_list).process()
        plot_dict = self.read_config_file(self.raw_file)

        return data_dict, plot_dict


if __name__ == '__main__':
    logging.basicConfig(
        # filename=os.path.join(
        #     os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename

    Tk().withdraw()
    raw_file_name = askopenfilename(
        title='Choose Poles Figure File...',
        filetypes=[("Raw files", "*.raw")]
    )
    logging.info("File {0} was chosen.".format(raw_file_name))
    if raw_file_name:
        sample = RawReader(raw_file_name)
        data_dict, plot_dict = sample.matrix_data()
        print(data_dict)
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
