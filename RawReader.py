import abc
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


class RsmScanRaw(DetectorScanRaw):
    def process(self):
        phi_data_matrix = np.asanyarray(
            [
                [float(j.split('=')[1]) for j in value if j.startswith('_PHI')]
                for value in self.data_list[1:]
                ]
        )
        omega_data_matrix = np.asanyarray(
            [
                [
                    float(j.split('=')[1]) for j in value
                    if j.startswith('_OMEGA')
                    ]
                for value in self.data_list[1:]
                ]
        )
        int_data_matrix = np.asanyarray(
            [
                [float(j.split('\t')[1]) for j in
                 value
                 if (re.match('\d', j) or j.startswith("-"))]
                for value in self.data_list[1:]
                ]
        )
        two_theta_data_matrix = np.asanyarray(
            [
                [float(j.split('\t')[0]) for j in
                 value
                 if (re.match('\d', j) or j.startswith("-"))]
                for value in self.data_list[1:]
                ]
        )

        return {
            'int_data': int_data_matrix,
            'phi_data': phi_data_matrix,
            'omega_data': omega_data_matrix,
            'two_theta_data': two_theta_data_matrix
        }


class RawReader(object):
    def __init__(self, raw_file):
        logging.info("Reading raw file {0}".format(raw_file))
        self.file = raw_file

    def get_head(self):
        data_list = self.__str_data(self.file)
        scan_dict = {i.split('=')[0].strip(): i.split('=')[1].strip()
                     for i in data_list[0] if i.startswith("_")}
        return scan_dict

    @staticmethod
    def __str_data(raw_file):
        ds = DatasetDiffractPlusV3(open(raw_file, 'rb'))
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

        def beam_intensity(directory, default_beam_intensity):
            logging.info("Starts calculation beam intensity.")

            source_file_list = [
                os.path.join(directory, 'beam1mm.raw'),
                os.path.join(directory, 'beam8mm.raw')]
            beam_int_list = [
                (max(RawReader(i).matrix_data()[0]['int_data']) * 8940)
                for i in source_file_list if os.path.isfile(i)]

            if not beam_int_list:
                logging.warning("Could not found default beam intensity files")
                beam_int_list = [default_beam_intensity]
                logging.info(
                    "Got Source Beam Intensity from config file.\n",
                    "Source Beam Intensity = {0}.".format(
                        beam_int_list[0])
                )
            beam_int_list = np.asarray(beam_int_list)
            beam_int = np.mean(beam_int_list)

            return beam_int

        config = configparser.ConfigParser()
        import sys
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(os.path.dirname(raw_file), 'config.ini')
        ])
        config['db']['sample'] = guess_sample_name(raw_file)
        config['db']['directory'] = os.path.dirname(raw_file)
        config['db']['file'] = os.path.basename(raw_file)

        config['db']['beam_intensity'] = str(
            beam_intensity(
                config['db']['directory'],
                float(config['db']['beam_intensity'])
            )
        )

        return config

    @abc.abstractmethod
    def matrix_data(self):
        data_list = self.__str_data(self.file)
        scan_dict = self.get_head()
        process_dict = {
            'SingleScanPlot': DetectorScanRaw,
            'TwoDPlot': TwoDScanRaw,
            'RSMPlot': RsmScanRaw
        }
        data_dict = process_dict[scan_dict['_TYPE']](data_list).process()
        if scan_dict['_TYPE'] == 'TwoDPlot':
            instance_dict = self.read_config_file(self.file)
        else:
            instance_dict = None

        return data_dict, instance_dict


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
        sample = RawReader()
        data_dict, plot_dict = sample.matrix_data(raw_file_name)
        print(data_dict)
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
