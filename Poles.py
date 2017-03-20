from __future__ import print_function, unicode_literals

import csv
import logging
import logging.config
import logging.handlers
import os
import re
import sys
from functools import wraps

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure

from bruker3 import DatasetDiffractPlusV3

try:
    import configparser
except ImportError:
    import ConfigParser as configparser
try:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename, askopenfilenames
except ImportError:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename, askopenfilenames

    print("Python3 Detected...")
else:
    print("Python2 Detected...")


class Square(object):
    def __init__(self, central_point_list, size_list, **param):
        self.central_point_list = central_point_list
        self.size_list = size_list
        self.param = param

        if ('limitation' in param) and param['limitation']:
            (y_limit, x_limit) = param['limitation']
        else:
            x_limit = 500000
            y_limit = 500000

        x_min = max(0, int(
            np.floor(self.central_point_list[1] - self.size_list[1] / 2.0)))
        x_max = min(
            int(np.floor(self.central_point_list[1] + self.size_list[1] / 2.0)),
            x_limit)
        y_min = max(0, int(
            np.floor(self.central_point_list[0] - self.size_list[0] / 2.0)))
        y_max = min(
            int(np.floor(self.central_point_list[0] + self.size_list[0] / 2.0)),
            y_limit)
        self.x_list = np.asarray([x_min, x_max, x_max, x_min, x_min])
        self.y_list = np.asarray([y_min, y_min, y_max, y_max, y_min])

        self.figure_handle = None

    def draw(self, **param):
        self.__init__(self.central_point_list, self.size_list, **self.param)
        self.figure_handle, = plt.plot(self.x_list, self.y_list, linewidth=0.5)

    def sum(self, intensity_matrix):
        intensity_result_matrix = intensity_matrix[
                                  self.y_list[0]:self.y_list[2],
                                  self.x_list[0]:self.x_list[1]
                                  ]
        peak_intensity_int = np.sum(intensity_result_matrix)
        b, p = intensity_result_matrix.shape
        peak_matrix_points = b * p

        return peak_intensity_int, peak_matrix_points

    def remove(self):
        axes = plt.gca()
        if self.figure_handle is not None and axes is not None:
            axes.lines.remove(self.figure_handle)

    def is_contained(self, point_list):
        if (
                            self.x_list[0] < point_list[0] < self.x_list[1] and
                            self.y_list[0] < point_list[1] < self.y_list[2]
        ):
            return True
        else:
            return False


class BeamIntensityFile(object):
    def __init__(self, raw_file):
        self.raw_file = raw_file
        self.data_dict = None

    def get_data(self):
        logging.info("Reading raw file {0}".format(self.raw_file))
        ds = DatasetDiffractPlusV3(open(self.raw_file, 'rb'))
        data_string = ds.pretty_format(print_header=True)
        data_list = data_string.split('\n')

        return data_list

    def raw_file_reader(self):
        data_list = self.get_data()
        int_data_matrix = [
            float(i.split()[1]) for i in data_list if re.match('\d', i)
            ]
        max_intensity = max(int_data_matrix) * 8940

        return max_intensity


class TwoDFigure(BeamIntensityFile):
    def __init__(self, raw_file):
        super(TwoDFigure, self).__init__(raw_file)

        config = configparser.ConfigParser()
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(os.path.dirname(raw_file), 'config.ini')
        ])
        self.preference_dict = dict(config._sections['db'])
        self.guess_sample_name()
        self.label = '2D'

    def beam_intensity(self):
        logging.info("Starts calculation beam intensity.")
        directory = os.path.dirname(self.raw_file)
        beam_intensity_float_list = []

        default_beam_intensity_file_list = [
            os.path.join(directory, 'beam1mm.raw'),
            os.path.join(directory, 'beam8mm.raw')]
        for i in default_beam_intensity_file_list:
            if os.path.isfile(i):
                j = BeamIntensityFile(i)
                beam_intensity_float_list.append(j.raw_file_reader())
                logging.info(
                    "Successfully load beam intensity file {0}".format(i)
                )

        if beam_intensity_float_list == []:
            logging.info("Could not found default beam intensity files")
            try:
                beam_intensity_float_list.append(
                    float(self.preference_dict['beam_intensity'])
                )
            except AttributeError:
                try:
                    from PyQt5.QtWidgets import QFileDialog as QFileDialog
                except ImportError:
                    root = Tk().withdraw()
                    beam_intensity_file_list = askopenfilenames(
                        title='Choose Intensity File...',
                        initialdir=directory,
                        filetypes=[("Raw files", "*.raw")]
                    )
                    beam_intensity_file_list = root.tk.splitlist(
                        beam_intensity_file_list
                    )
                else:
                    from PyQt5.QtWidgets import QApplication

                    app = QApplication([])
                    beam_intensity_file_list = QFileDialog.getOpenFileNames(
                        caption='Choose Intensity File...',
                        directory=directory,
                        filter="RAW files (*.raw)"
                    )
                    beam_intensity_file_list = beam_intensity_file_list[0]
                finally:
                    if beam_intensity_file_list is not []:
                        for i in beam_intensity_file_list:
                            j = BeamIntensityFile(i)
                            beam_intensity_float_list.append(
                                j.raw_file_reader())
                            logging.debug(
                                ("Successfully load beam"
                                 " intensity file {0}").format(i)
                            )
                    else:
                        print("Intensity Error!")
                        logging.warning(
                            "Failed to load beam intensity, use "
                            "16000*8940 instead."
                        )
                        beam_intensity_float_list.append(16000 * 8940)
            else:
                logging.debug(
                    "Successfully load beam intensity from config"
                )

        beam_intensity_float_list = np.mean(beam_intensity_float_list)

        return beam_intensity_float_list

    @staticmethod
    def correction(chi, thickness):
        theta = 14.22  # for GaP
        omega = 14.22  # for GaP
        if thickness is not None:
            thickness /= 10000.
        else:
            thickness = 900. / 10000.
        e_angle = 90. - np.rad2deg(
            np.arccos(
                np.cos(np.deg2rad(chi)) * np.sin(np.deg2rad(theta))))
        i_angle = 90. - np.rad2deg(
            np.arccos(
                np.cos(np.deg2rad(chi)) * np.sin(np.deg2rad(omega))))
        offset = e_angle - i_angle
        eta = 1. / 37.6152  # 1/um at 8.05keV (CXRO)
        p = np.sin(np.deg2rad(e_angle + offset))
        q = np.sin(np.deg2rad(e_angle - offset))
        coefficient_b = p / (p + q)
        coefficient_c = 1. - np.exp(-eta * thickness * (1. / p + 1. / q))
        coefficient = coefficient_b * (1. / eta) * coefficient_c
        logging.debug(
            "eta:{0}, thickness:{1}, p:{2}, q:{3}, c_c:{4}".format(
                eta, thickness, p, q, coefficient_c
            )
        )

        return coefficient

    def guess_sample_name(self):
        try:
            sample_name = (
                re.findall(
                    r'S\d\d\d\d',
                    os.path.abspath(self.raw_file),
                    flags=re.IGNORECASE
                )
            )
        except FileNotFoundError:
            print("Error, no sample name.")
        else:
            if sample_name is not []:
                self.preference_dict['sample'] = sample_name[-1]

    def raw_file_reader(self):
        data_list = self.get_data()
        khi_data_matrix = np.asarray(
            [
                float(i.split("=")[1].strip('\n')) for i in data_list
                if i.startswith("_KHI")
                ]
        )
        step_time_set = {
            float(i.split("=")[1].strip('\n')) for i in data_list
            if i.startswith("_STEPTIME")
            }
        scan_data_matrix = np.asarray(
            [
                list(map(float, i.split())) for i in data_list
                if (re.match('\d', i) or i.startswith("-"))
                ]
        )
        matrix_width_int = len(khi_data_matrix)
        logging.debug("Matrix width is {0}".format(matrix_width_int))
        int_data_matrix = scan_data_matrix[:, 1] / next(iter(step_time_set))
        phi_data_matrix = scan_data_matrix[:, 0]
        phi_data_matrix = phi_data_matrix[
                          0:int(len(int_data_matrix) / matrix_width_int)
                          ]
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

        shift the intensity_matrix fig  towards left for alpha degree and
        add the part which exceeded to the right.

        Args:
        intensity_matrix: input array.
        alpha: shift degree.

        Returns:
        intensity_new_matrix: array shifted.
        """
        _, x_length = intensity_matrix.shape
        intensity_new_matrix = np.zeros(shape=intensity_matrix.shape)
        if alpha is not 0:
            intensity_new_matrix[
            :(x_length - alpha) % x_length, :
            ] = intensity_matrix[(x_length + alpha) % x_length:, :]
            intensity_new_matrix[
            (x_length - alpha) % x_length:, :
            ] = intensity_matrix[:(x_length + alpha) % x_length, :]
        else:
            intensity_new_matrix = intensity_matrix

        return intensity_new_matrix

    def save_config(self):
        config = configparser.ConfigParser()
        directory = os.path.dirname(self.raw_file)
        self.preference_dict['directory'] = os.path.abspath(directory)
        config['db'] = self.preference_dict
        with open(os.path.join(directory, 'config.ini'), 'w') as configfile:
            config.write(configfile)

    def _print_log(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            logging.info("Start to %s." % func.__name__)
            result = func(*args, **kwargs)
            logging.info(
                "{0} Finished!\n----------------------------------"
                "------------------------------".format(func.__name__)
            )
            return result

        return wrapper

    @_print_log
    def plot(self, **param):
        """Plot the 2D Image.
        Args:
            is_show_image: Boolean value control if show the image
        Returns:
            None
        """
        if self.data_dict is None:
            data_dict = self.raw_file_reader()
            self.data_dict = data_dict
        else:
            data_dict = self.data_dict

        ax2d = plt.gca()
        im = ax2d.imshow(
            data_dict['int_data'],
            origin="lower",
            norm=LogNorm(
                vmin=int(self.preference_dict['v_min']),
                vmax=int(self.preference_dict['v_max'])
            )
        )
        ax2d.tick_params(axis='both', which='major', labelsize=16)
        plt.colorbar(
            im, ax=ax2d, extend='max', fraction=0.03, pad=0.04, shrink=0.6
        )
        self.display_save_image(**param)

    def display_save_image(self, **param):
        sample_name_str = self.preference_dict['sample']

        plt.title(sample_name_str + self.label + "\n")

        if 'is_show_image' in param and param['is_show_image']:
            plt.show()

        figure_file_name = os.path.join(
            os.path.dirname(self.raw_file),
            sample_name_str + '_' + self.label + '.png'
        )

        if 'is_save_image' in param and not param['is_save_image']:
            pass
        else:
            logging.info(
                "Saving figure file to {0}".format(figure_file_name))
            plt.savefig(
                figure_file_name,
                dpi=200,
                bbox_inches='tight')


class PolarFigure(TwoDFigure):
    def __init__(self, raw_file):
        super(PolarFigure, self).__init__(raw_file)
        self.label = 'Polar'

    @TwoDFigure._print_log
    def plot(self, is_log_scale=1, **param):
        """Plot the Polar coordinate system Image.
            Args:
                is_show_image: Boolean value control if show the image
                is_log_scale: Boolean value control if plot in log scale
            Returns:
                None
        """
        if self.data_dict is None:
            data_dict = self.raw_file_reader()
            self.data_dict = data_dict
        else:
            data_dict = self.data_dict

        ctn_number = int(self.preference_dict['ctn_number'])
        v_max = int(self.preference_dict['v_max'])
        v_min = int(self.preference_dict['v_min'])
        step = (v_max - v_min) / ctn_number
        contour_levels = np.arange(v_min, v_max, step)
        int_data_matrix = data_dict['int_data'].T

        phi_data_matrix = np.radians(
            data_dict['phi_data'] +
            [float(self.preference_dict['phi_offset'])]
        )
        logging.debug(
            "Add {0} degree offset to phi".format(
                self.preference_dict['phi_offset']
            )
        )
        r, theta = np.meshgrid(data_dict['khi_data'], phi_data_matrix)

        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.grid(color="white")
        ax.tick_params(axis="y", labelsize=15, labelcolor="white")
        ax.tick_params(axis="x", labelsize=15)
        ax.set_rmin(0.)

        if is_log_scale:
            cax = ax.contourf(
                theta,
                r,
                int_data_matrix,
                contour_levels,
                norm=LogNorm(
                    vmin=v_min,
                    vmax=v_max)
            )

        else:
            cax = ax.contourf(
                theta,
                r,
                int_data_matrix,
                contour_levels
            )

        color_bar = fig.colorbar(cax, pad=0.1, format="%.f", extend='max')
        color_bar.ax.tick_params(labelsize=16)
        color_bar.set_label(r'$Counts\ per\ second $', size=15)
        self.display_save_image(**param)


class MeasureFigure(TwoDFigure):
    def __init__(self, raw_file):
        super(MeasureFigure, self).__init__(raw_file)
        self.label = 'measure'

    @classmethod
    def peak_search(cls, int_data):
        def sort_index_list(index_list):
            logging.debug("index list before sort:{0}".format(index_list))
            khi_sorted_list = sorted(index_list, key=lambda pair: pair[1])
            chi_index_list = [l[0] for l in khi_sorted_list]
            shifted_index_int = chi_index_list.index(max(chi_index_list))
            sorted_index_list = [None] * len(index_list)
            for l in range(len(index_list)):
                sorted_index_list[l - shifted_index_int] = khi_sorted_list[l]
            logging.debug(
                "index list after sort:{0}".format(sorted_index_list)
            )
            return sorted_index_list

        neighborhood = generate_binary_structure(2, 2)
        for i in range(3):
            int_data = gaussian_filter(int_data, 4, mode='nearest')

        local_max = (
            maximum_filter(int_data, footprint=neighborhood) == int_data
        )
        index = np.asarray(np.where(local_max))
        filtered_index_list = [
            [i, j] for (i, j) in zip(index[0, :], index[1, :])
            ]

        chi_threshold = 40
        filtered_index_list = [
            i for i in filtered_index_list if i[0] < chi_threshold
            ]

        peak_intensity_matrix, points, _ = cls.square(
            int_data, filtered_index_list, [10, 10], is_plot_square=0
        )
        peak_intensity_matrix_bigger, points_more, _ = cls.square(
            int_data, filtered_index_list, [20, 20], is_plot_square=0
        )
        peak_intensity_matrix = (
            peak_intensity_matrix -
            (
                (peak_intensity_matrix_bigger - peak_intensity_matrix) /
                (points_more - points) *
                points
            )
        )
        peak_intensity_matrix = list(peak_intensity_matrix)
        peak_intensity_matrix_copy = peak_intensity_matrix.copy()
        new_filtered_index_list = []
        for i in range(4):
            index = peak_intensity_matrix.index(
                max(peak_intensity_matrix_copy)
            )
            new_filtered_index_list.append(
                filtered_index_list[index]
            )
            peak_intensity_matrix_copy.pop(
                peak_intensity_matrix_copy.index(
                    max(peak_intensity_matrix_copy)
                )
            )
        filtered_index_list = new_filtered_index_list

        filtered_index_list = sort_index_list(filtered_index_list)
        return filtered_index_list

    @staticmethod
    def square(intensity_matrix, index_list, size_list, is_plot_square=1):
        logging.debug("size_list: {0}".format(size_list))
        peak_intensity_list = []
        peak_matrix_points = []
        square_instances = []
        logging.info("The square size is {0}".format(size_list))

        for i in index_list:
            square_instance = Square(
                i, size_list, limitation=intensity_matrix.shape
            )
            square_instances.append(square_instance)

        for i in square_instances:
            if is_plot_square:
                i.draw()

            peak_intensity_int, peak_matrix_point_int = i.sum(
                intensity_matrix
            )
            peak_intensity_list.append(peak_intensity_int)
            peak_matrix_points.append(peak_matrix_point_int)

        peak_intensity_matrix = np.asanyarray(peak_intensity_list)
        peak_matrix_points_matrix = np.asanyarray(peak_matrix_points)

        return (
            peak_intensity_matrix,
            peak_matrix_points_matrix,
            square_instances
        )

    @TwoDFigure._print_log
    def plot(self, is_calculation=0, **param):
        result = self.int_points(is_plot_square=(not is_calculation), **param)
        if not is_calculation:
            super(MeasureFigure, self).plot(**param)
        return result

    def int_points(self, is_plot_square=1, **param):
        if self.data_dict is None:
            data_dict = self.raw_file_reader()
            self.data_dict = data_dict
        else:
            data_dict = self.data_dict

        if "outer_index_list" in param and param['outer_index_list']:
            index_list = param['outer_index_list']
        else:
            index_list = self.peak_search(data_dict['int_data'])
            logging.debug(index_list)

        size_list = list(
            map(int, self.preference_dict['square_size'].split(',')))

        (inner_peak_intensity_matrix,
         inner_peak_matrix_points_matrix,
         inner_square_instances) = self.square(
            data_dict['int_data'],
            index_list,
            size_list,
            is_plot_square=is_plot_square
        )
        (outer_peak_intensity_matrix,
         outer_peak_matrix_points_matrix,
         outer_square_instances) = self.square(
            data_dict['int_data'],
            index_list,
            [i + 4 for i in size_list],
            is_plot_square=is_plot_square
        )

        square_instances = inner_square_instances + outer_square_instances

        background_noise_intensity_float = np.average(
            (outer_peak_intensity_matrix - inner_peak_intensity_matrix) /
            (outer_peak_matrix_points_matrix - inner_peak_matrix_points_matrix)
        )
        peak_net_intensity_matrix = (
            inner_peak_intensity_matrix -
            background_noise_intensity_float * inner_peak_matrix_points_matrix)

        del data_dict

        result = {
            'peak_intensity_matrix': peak_net_intensity_matrix,
            'index': index_list,
            'square_instances': square_instances
        }

        return result

    def mt_intensity_to_fraction(self, result):
        index_list = result['index']
        peak_net_intensity_matrix = result['peak_intensity_matrix']
        logging.info("Peak intensity is {0}".format(peak_net_intensity_matrix))
        try:
            thickness = float(self.preference_dict['thickness'])
        except KeyError:
            thickness = None
            logging.warning("Thickness was not provided, use 90nm instead")
        eta = [self.correction(x[0], thickness) for x in index_list]
        logging.info("Intensity correction index eta is {0}".format(eta))
        coefficient_list = np.asarray([
            0.939691064,
            0.666790711 / (eta[1] / eta[0]),
            0.426843274 / (eta[2] / eta[0]),
            0.72278158 / (eta[3] / eta[0])
        ])
        beam_intensity_float = self.beam_intensity()
        peak_net_intensity_matrix = (
            peak_net_intensity_matrix *
            10000 *
            coefficient_list /
            beam_intensity_float
        )
        result = {
            'peak_intensity_matrix': peak_net_intensity_matrix,
            'index': index_list
        }
        return result

    @TwoDFigure._print_log
    def print_result_csv(self):
        result = self.plot()
        result = self.mt_intensity_to_fraction(result)

        mt_table_file = os.path.join(
            os.path.dirname(self.raw_file),
            '{0}_result.csv'.format(self.preference_dict['sample'])
        )

        with open(mt_table_file, 'w') as tableTeX:
            spam_writer = csv.writer(tableTeX, dialect='excel')
            spam_writer.writerow(['MT-A', 'MT-D', 'MT-C', 'MT-B'])
            spam_writer.writerow(result['peak_intensity_matrix'])


if __name__ == '__main__':
    logging.basicConfig(
        # filename=os.path.join(
        #     os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )

    Tk().withdraw()
    raw_file_name = askopenfilename(
        title='Choose Poles Figure File...',
        filetypes=[("Raw files", "*.raw")]
    )
    logging.info("File {0} was chosen.".format(raw_file_name))
    sample = MeasureFigure(raw_file_name)
    # sample.plot_polar_image()
    sample.print_result_csv()
    sample.save_config()
    # sample.plot_2d_measurement(is_show_image=True)
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
