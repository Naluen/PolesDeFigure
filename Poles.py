import configparser
import logging
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure

from bruker3 import DatasetDiffractPlusV3
from collections import deque


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


class Square(object):
    def __init__(
            self, central_point_list, size_list, limitation=(500000, 500000)):
        self.central_point_list = central_point_list
        self.size_list = size_list

        (y_limit, x_limit) = limitation

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

    def draw(self, limitation=(500000, 500000)):
        self.__init__(self.central_point_list, self.size_list, limitation)
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
        """
        Remove the the square lines.
        :return: None
        """
        axes = plt.gca()
        try:
            axes.lines.remove(self.figure_handle)
        except ValueError:
            logging.warning("Could not find the square.")

    def __contains__(self, item):
        """
        Check if the point is in the square.
        :param item: the position of point [x,y].
        :return: The boolean value.
        """
        if (
                            self.x_list[0] < item[0] < self.x_list[1] and
                            self.y_list[0] < item[1] < self.y_list[2]
        ):
            return True
        else:
            return False


class PrintLogDecorator(object):
    def __init__(self, *args, **kwargs):
        # store arguments passed to the decorator
        self.args = args
        self.kwargs = kwargs

    def __call__(self, func):
        def new_func(*args, **kwargs):
            logging.info("Start to %s." % func.__name__)
            # call the method
            ret = func(*args, **kwargs)
            logging.info(
                "{0} Finished!\n----------------------------------"
                "------------------------------".format(func.__name__)
            )

            return ret

        new_func.__doc__ = func.__doc__
        return new_func


class TwoDFigureRaw(object):
    def __init__(self, raw_file):
        self.plot_dict = None
        self.label = '2D'
        self.raw_file = raw_file

    @staticmethod
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

    @classmethod
    def raw_reader(cls, raw_file):
        # Get Config.
        config = configparser.ConfigParser()
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(os.path.dirname(raw_file), 'config.ini')
        ])
        sample_name = cls.guess_sample_name(os.path.abspath(raw_file))
        if sample_name:
            config['db']['sample'] = sample_name
        config['db']['raw_file'] = os.path.basename(raw_file)
        config['db']['directory'] = os.path.dirname(raw_file)
        # Get Data.
        logging.info("Reading raw file {0}".format(raw_file))
        ds = DatasetDiffractPlusV3(open(raw_file, 'rb'))
        data_string = ds.pretty_format(print_header=True)
        data_list = data_string.split('\n')
        from itertools import groupby
        data_list = [list(g) for k, g in
                     groupby((line.strip() for line in data_list), bool) if k]
        scan_dict = {i.split('=')[0].strip(): i.split('=')[1].strip()
                     for i in data_list[0] if i.startswith("_")}

        process_dict = {
            'SingleScanPlot': DetectorScanRaw,
            'TwoDPlot': TwoDScanRaw
        }
        data_dict = process_dict[scan_dict['_TYPE']](data_list).process()
        data_dict['plot_dict'] = config

        return data_dict

    @staticmethod
    def save_config(plot_dict):
        directory = plot_dict['db']['directory']
        with open(os.path.join(directory, 'config.ini'), 'w') as configfile:
            plot_dict.write(configfile)

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0):
        data_dict = self.raw_reader(self.raw_file)
        self.plot_dict = data_dict['plot_dict']

        v_max = int(self.plot_dict['db']['v_max'])
        v_min = int(self.plot_dict['db']['v_min'])

        ax2d = plt.gca()
        im = ax2d.imshow(
            data_dict['int_data'],
            origin="lower",
            norm=LogNorm(vmin=v_min, vmax=v_max)
        )
        ax2d.tick_params(axis='both', which='major', labelsize=16)

        ticks = np.logspace(0, np.log10(v_max), np.log10(v_max) + 1)
        plt.colorbar(im, pad=0.1, format="%.e", extend='max', ticks=ticks)
        self.display_save_image(is_show, is_save)

    def display_save_image(self, is_show, is_save):
        sample_name_str = self.plot_dict['db']['sample']
        plt.title(sample_name_str + " " + self.label + "\n")

        if is_show:
            plt.show()

        figure_file_name = os.path.join(
            self.plot_dict['db']['directory'],
            sample_name_str + '_' + self.label + '.png'
        )

        if is_save:
            logging.info(
                "Saving figure file to {0}".format(figure_file_name))
            plt.savefig(
                figure_file_name,
                dpi=200,
                bbox_inches='tight')


class PolarFigureRaw(TwoDFigureRaw):
    def __init__(self, raw_file):
        super(PolarFigureRaw, self).__init__(raw_file)
        self.label = 'Polar'

    def raw_reader(self, raw_file):
        data_dict = TwoDFigureRaw.raw_reader(raw_file)
        plot_dict = data_dict['plot_dict']
        print(data_dict['phi_data'])
        data_dict['phi_data'] = np.radians(
            data_dict['phi_data'] +
            [float(plot_dict['db']['phi_offset'])]
        )
        print(data_dict['phi_data'])
        return data_dict

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0):
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.grid(color="white")
        ax.tick_params(axis="y", labelsize=15, labelcolor="white")
        ax.tick_params(axis="x", labelsize=15)
        ax.set_rmin(0.)
        data_dict = self.raw_reader(self.raw_file)
        self.plot_dict = data_dict['plot_dict']

        v_max = int(self.plot_dict['db']['v_max'])
        v_min = int(self.plot_dict['db']['v_min'])
        [xx, yy] = np.meshgrid(
            data_dict['khi_data'][:,0].T, data_dict['phi_data'])
        plt.coutourf(xx, yy, data_dict['int_data'])

        ax=plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=16)

        ticks = np.logspace(0, np.log10(v_max), np.log10(v_max) + 1)
        plt.colorbar(im, pad=0.1, format="%.e", extend='max', ticks=ticks)
        self.display_save_image(is_show, is_save)

class MeasureFigureRaw(TwoDFigureRaw):
    def __init__(self, raw_file):
        self.data_dict = self.raw_reader(raw_file)
        super(MeasureFigureRaw, self).__init__(raw_file)
        self.label = 'measure'

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

    @classmethod
    def peak_search(cls, int_data):
        def sort_index_list(index_list):
            logging.debug("index list before sort:{0}".format(index_list))
            khi_sorted_list = sorted(index_list, key=lambda pair: pair[1])
            chi_index_list = [l[0] for l in khi_sorted_list]
            shifted_index_int = chi_index_list.index(max(chi_index_list))
            khi_deque= deque(khi_sorted_list)
            khi_deque.rotate(-shifted_index_int)
            sorted_index_list = list(khi_deque)
            logging.info(
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
            int_data, filtered_index_list, [10, 10], is_plot=0
        )
        peak_intensity_matrix_bigger, points_more, _ = cls.square(
            int_data, filtered_index_list, [20, 20], is_plot=0
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

        filtered_index_list = [
            x for (y,x)
            in sorted(zip(peak_intensity_matrix,filtered_index_list),
                      key=lambda pair: pair[0])][-4:]

        filtered_index_list = sort_index_list(filtered_index_list)
        return filtered_index_list

    @staticmethod
    def square(intensity_matrix, index_list, size_list, is_plot):
        logging.debug("size_list: {0}".format(size_list))
        peak_intensity_list = []
        peak_matrix_points = []
        square_instances = []
        logging.debug("The square size is {0}".format(size_list))

        for i in index_list:
            square_instance = Square(
                i, size_list, limitation=intensity_matrix.shape
            )
            square_instances.append(square_instance)

        for i in square_instances:
            if is_plot:
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

    def beam_intensity(self):
        logging.info("Starts calculation beam intensity.")
        directory = self.plot_dict['db']['directory']

        default_beam_intensity_file_list = [
            os.path.join(directory, 'beam1mm.raw'),
            os.path.join(directory, 'beam8mm.raw')]
        beam_intensity_float_list = [max(self.raw_reader(i)['int_data']) * 8940
                                     for i in default_beam_intensity_file_list
                                     if os.path.isfile(i)]
        if not beam_intensity_float_list:
            logging.warning("Could not found default beam intensity files")
            beam_intensity_float_list = [
                float(self.plot_dict['db']['beam_intensity'])]
            logging.info(
                "Got Source Beam Intensity from config file.\n",
                "Source Beam Intensity = {0}.".format(
                    beam_intensity_float_list[0])
            )
        beam_intensity_float_list = np.asarray(beam_intensity_float_list)
        beam_intensity_float_list = np.mean(beam_intensity_float_list)

        return beam_intensity_float_list

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0, **param):
        data_dict = self.data_dict
        plot_dict = data_dict['plot_dict']
        self.plot_dict = data_dict['plot_dict']

        if "outer_index_list" in param and param['outer_index_list']:
            index_list = param['outer_index_list']
        else:
            index_list = self.peak_search(data_dict['int_data'])
            logging.debug(index_list)

        size_list = list(
            map(int, plot_dict['db']['square_size'].split(',')))
        (inner_peak_intensity_matrix,
         inner_peak_matrix_points_matrix,
         inner_square_instances) = self.square(
            data_dict['int_data'],
            index_list,
            size_list,
            is_plot=is_show
        )
        (outer_peak_intensity_matrix,
         outer_peak_matrix_points_matrix,
         outer_square_instances) = self.square(
            data_dict['int_data'],
            index_list,
            [i + 4 for i in size_list],
            is_plot=is_show
        )
        square_instances = inner_square_instances + outer_square_instances

        background_noise_intensity_float = np.average(
            (outer_peak_intensity_matrix - inner_peak_intensity_matrix) /
            (outer_peak_matrix_points_matrix - inner_peak_matrix_points_matrix)
        )
        peak_net_intensity_matrix = (
            inner_peak_intensity_matrix -
            background_noise_intensity_float * inner_peak_matrix_points_matrix)

        self.display_save_image(is_show=0, is_save=is_save)

        result = {
            'peak_intensity_matrix': peak_net_intensity_matrix,
            'index': index_list,
            'square_instances': square_instances,
            'dict': plot_dict
        }

        return result

    def mt_intensity_to_fraction(self, result):
        index_list = result['index']
        peak_net_intensity_matrix = result['peak_intensity_matrix']
        logging.info("Peak intensity is {0}".format(peak_net_intensity_matrix))

        thickness= float(result['dict']['db']['thickness'])
        logging.info("Sample thickness is {0}\n".format(thickness))

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

        result['peak_intensity_matrix']= peak_net_intensity_matrix
        return result

    @PrintLogDecorator()
    def print_result_csv(self):
        result = self.plot()
        result = self.mt_intensity_to_fraction(result)

        directory = result['dict']['db']['directory']
        sample = result['dict']['db']['sample']

        mt_table_file = os.path.join(
            directory,
            '{0}_result.csv'.format(sample)
        )
        import csv
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
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename

    Tk().withdraw()
    raw_file_name = askopenfilename(
        title='Choose Poles Figure File...',
        filetypes=[("Raw files", "*.raw")]
    )
    logging.info("File {0} was chosen.".format(raw_file_name))
    if raw_file_name:
        sample = PolarFigureRaw(raw_file_name)
        sample.plot(is_show=1)
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
