from __future__ import print_function, unicode_literals

import logging
import logging.config
import logging.handlers
import os
import re
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure
from scipy.ndimage.filters import gaussian_filter

import numpy as np
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


class BeamIntensityFile(object):
    def __init__(self, raw_file):
        self.raw_file = raw_file

    def get_data(self):
        logging.info("Reading raw file {0}".format(self.raw_file))
        ds = DatasetDiffractPlusV3(open(self.raw_file, 'rb'))
        data_string = ds.pretty_format(print_header=True)
        data_list = data_string.split('\n')

        return data_list

    def raw_file_reader(self):
        data_list = self.get_data()
        int_data_matrix = [float(i.split()[1]) for i in data_list if re.match('\d', i)]
        max_intensity = max(int_data_matrix) * 8940

        return max_intensity


class PolesFigureFile(BeamIntensityFile):
    def __init__(self, raw_file):
        super(PolesFigureFile, self).__init__(raw_file)

        self.directory = os.path.dirname(self.raw_file)
        config = configparser.ConfigParser()
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(self.directory, 'config.ini')
        ])
        self.preference_dict = dict(config._sections['db'])

    def raw_file_reader(self):
        data_list = self.get_data()
        khi_data_matrix = np.asarray(
            [float(i.split("=")[1].strip('\n')) for i in data_list if i.startswith("_KHI")]
        )
        step_time_set = {float(i.split("=")[1].strip('\n')) for i in data_list if i.startswith("_STEPTIME")}
        scan_data_matrix = np.asarray(
            [list(map(float, i.split())) for i in data_list if (re.match('\d', i) or i.startswith("-"))]
        )
        matrix_width_int = len(khi_data_matrix)
        logging.debug("Matrix width is {0}".format(matrix_width_int))
        phi_data_matrix = scan_data_matrix[:, 0][:matrix_width_int]
        int_data_matrix = scan_data_matrix[:, 1] / next(iter(step_time_set))
        int_data_matrix = np.reshape(
            int_data_matrix,
            (matrix_width_int, int(len(int_data_matrix) / matrix_width_int)))
        if len(step_time_set) is not 1:
            logging.warning("Step time is not identical.")
        del data_list

        return {
            'int_data': int_data_matrix,
            'phi_data': phi_data_matrix,
            'khi_data': khi_data_matrix
        }

    @staticmethod
    def _shift_angel(intensity_matrix, alpha):
        """Shift the intensity_matrix fig.

        shift the intensity_matrix fig  towards left for alpha degree and add the part
        which exceeded to the right.

        Args:
        intensity_matrix: input array.
        alpha: shift degree.

        Returns:
        intensity_new_matrix: array shifted.
        """
        _, x_length = intensity_matrix.shape
        intensity_new_matrix = np.zeros(shape=intensity_matrix.shape)
        if alpha is not 0:
            intensity_new_matrix[:(x_length - alpha) % x_length, :] = intensity_matrix[
                                                                      (x_length + alpha) % x_length:, :]
            intensity_new_matrix[(x_length - alpha) % x_length:, :] = intensity_matrix[
                                                                      :(x_length + alpha) % x_length, :]
        else:
            intensity_new_matrix = intensity_matrix

        return intensity_new_matrix

    @staticmethod
    def peak_search(int_data):
        neighborhood = generate_binary_structure(2, 1)
        for i in range(3):
            int_data = gaussian_filter(int_data, 4, mode='nearest')
        local_max = maximum_filter(int_data, footprint=neighborhood) == int_data
        index = np.asarray(np.where(local_max))
        x_threshold = [0, 70]
        y_threshold = [0, 360]
        filtered_index_list = [[i, j] for (i, j) in zip(index[0, :], index[1, :])]
        while len(filtered_index_list) > 8:
            if x_threshold[1] <= x_threshold[0]:
                logging.warning("Could not find the peak automatically.")
                break
            filtered_index_list = [
                [i, j] for (i, j, v, u, w, x) in zip(
                    index[0, :],
                    index[1, :],
                    index[1, :] > y_threshold[0],
                    index[1, :] < y_threshold[1],
                    index[0, :] > x_threshold[0],
                    index[0, :] < x_threshold[1]
                ) if all([v, u, w, x])]

            x_threshold = [x_threshold[0] + 1, x_threshold[1] - 1]
            y_threshold = [y_threshold[0] + 1, y_threshold[1] - 1]
        filtered_index_list = [i for i in filtered_index_list if i[0] < 40]
        return filtered_index_list

    @staticmethod
    def square(intensity_matrix, index_list, size_list, axes):
        logging.debug("size_list: {0}".format(size_list))
        peak_intensity_list = []
        peak_matrix_points = []
        logging.debug("The square size is {0}".format(size_list))
        for i in index_list:
            x_min = int(np.floor(i[1] - size_list[1] / 2.0))
            x_max = int(np.floor(i[1] + size_list[1] / 2.0))
            y_min = int(np.floor(i[0] - size_list[0] / 2.0))
            y_max = int(np.floor(i[0] + size_list[0] / 2.0))
            x = np.asarray([x_min, x_max, x_max, x_min, x_min])
            y = np.asarray([y_min, y_min, y_max, y_max, y_min])
            axes.plot(x, y, linewidth=0.5)

            logging.debug("x_min\tx_max\ty_min\ty_max:\n{0}\t{1}\t{2}\t{3}".format(x_min, x_max, y_min, y_max))
            intensity_result_matrix = intensity_matrix[y_min:y_max, x_min:x_max]
            peak_intensity_list.append(np.sum(intensity_result_matrix))
            b, p = intensity_result_matrix.shape
            peak_matrix_points.append(b * p)

        peak_intensity_matrix = np.asanyarray(peak_intensity_list)
        peak_matrix_points_matrix = np.asanyarray(peak_matrix_points)
        return peak_intensity_matrix, peak_matrix_points_matrix

    def plot_2d_image(self, is_show_image=True, is_log_scale=False):
        """Plot the 2D Image"""
        data_dict = self.raw_file_reader()

        # data_dict['int_data'] = self._shift_angel(data_dict['int_data'], 360 - int(data_dict['phi_data'][-1]))
        data_dict['phi_data'] = np.radians(
            data_dict['phi_data'] +
            [float(self.preference_dict['phi_offset'])]
        )
        logging.debug("Add {0} degree offset to phi".format(self.preference_dict['phi_offset']))
        r, theta = np.meshgrid(data_dict['khi_data'], data_dict['phi_data'])

        plt.figure(figsize=(25, 5))
        ax2d = plt.subplot(111)
        im = ax2d.imshow(
            data_dict['int_data'],
            origin="lower",
            norm=LogNorm(
                vmin=int(self.preference_dict['v_min']),
                vmax=int(self.preference_dict['v_max'])
            )
        )
        plt.colorbar(im, ax=ax2d, extend='max', fraction=0.046, pad=0.04)
        plt.title(self.preference_dict['sample'] + "\n")

        plt.figure(figsize=(25, 5))
        ax_2d_calculation = plt.subplot(111)
        index_list = self.peak_search(data_dict['int_data'])
        size_list = list(map(int, self.preference_dict['square_size'].split(',')))
        ax_2d_calculation.imshow(
            data_dict['int_data'],
            origin="lower",
            norm=LogNorm(
                vmin=int(self.preference_dict['v_min']),
                vmax=int(self.preference_dict['v_max'])
            )
        )
        inner_peak_intensity_matrix, inner_peak_matrix_points_matrix = self.square(
            data_dict['int_data'], index_list, size_list, ax_2d_calculation)
        outer_peak_intensity_matrix, outer_peak_matrix_points_matrix = self.square(
            data_dict['int_data'], index_list, [i + 4 for i in size_list], ax_2d_calculation)

        background_noise_intensity_float = np.average(
            (outer_peak_intensity_matrix - inner_peak_intensity_matrix) /
            (outer_peak_matrix_points_matrix - inner_peak_matrix_points_matrix)
        )
        peak_net_intensity_matrix = (
            inner_peak_intensity_matrix -
            background_noise_intensity_float * inner_peak_matrix_points_matrix)
        if is_show_image:
            plt.show(all)


if __name__ == '__main__':
    logging.basicConfig(
        filename=os.path.join(os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )

    Tk().withdraw()
    raw_file_name = askopenfilename(
        title='Choose Poles Figure File...',
        filetypes=[("Raw files", "*.raw")]
    )
    logging.info("File {0} was chosen.".format(raw_file_name))
    sample = PolesFigureFile(raw_file_name)
    sample.plot_2d_image()
