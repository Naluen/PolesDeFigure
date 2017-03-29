import logging
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import deque
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure

from Reader import Reader as Reader
from Square import Square as Square


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
    def __init__(self, file):
        self.label = '2D'
        self.file = file

    @staticmethod
    def save_config(plot_dict):
        directory = plot_dict['db']['directory']
        with open(os.path.join(directory, 'config.ini'), 'w') as configfile:
            plot_dict.write(configfile)
        import sys
        with open(
                os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'), 'w'
        ) as configfile:
            plot_dict.write(configfile)

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0):
        data_dict, plot_dict = Reader(self.file).matrix_data()

        plot_dict['db']['directory'] = os.path.dirname(self.file)
        plot_dict['db']['file'] = os.path.basename(self.file)

        v_max = int(plot_dict['db']['v_max'])
        v_min = int(plot_dict['db']['v_min'])

        ax2d = plt.gca()
        im = ax2d.imshow(
            data_dict['int_data'],
            origin="lower",
            norm=LogNorm(vmin=v_min, vmax=v_max)
        )
        ax2d.tick_params(axis='both', which='major', labelsize=16)

        ticks = np.logspace(0, np.log10(v_max), np.log10(v_max) + 1)
        plt.colorbar(im, pad=0.1, format="%.e", extend='max', ticks=ticks)

        self.display_save_image(plot_dict, is_show, is_save)
        self.save_config(plot_dict)

    def display_save_image(self, plot_dict, is_show=0, is_save=1):
        sample_name_str = plot_dict['db']['sample']
        plt.title(sample_name_str + " " + self.label + "\n")

        if is_show:
            plt.show()

        figure_file_name = os.path.join(
            plot_dict['db']['directory'],
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
    def __init__(self, file):
        super(PolarFigureRaw, self).__init__(file)
        self.label = 'Polar'

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0):
        data_dict, plot_dict = Reader(self.file).matrix_data()

        data_dict['phi_data'] = np.radians(
            data_dict['phi_data'] +
            [float(plot_dict['db']['phi_offset'])]
        )

        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

        ax.tick_params(axis="y", labelsize=15, labelcolor="white")
        ax.tick_params(axis="x", labelsize=15)
        ax.set_rmin(0.)

        [xx, yy] = np.meshgrid(data_dict['phi_data'],
                               data_dict['khi_data'][:, 0].T)
        v_max = int(plot_dict['db']['v_max'])
        v_min = int(plot_dict['db']['v_min'])
        im = ax.pcolormesh(xx, yy, data_dict['int_data'],
                           norm=LogNorm(vmin=v_min, vmax=v_max))
        ax.grid(color="white")

        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=16)

        ticks = np.logspace(0, np.log10(v_max), np.log10(v_max) + 1)
        plt.colorbar(im, pad=0.1, format="%.e", extend='max', ticks=ticks)
        self.display_save_image(plot_dict, is_show, is_save)
        self.save_config(plot_dict)


class MeasureFigureRaw(TwoDFigureRaw):
    def __init__(self, file):
        super(MeasureFigureRaw, self).__init__(file)
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
            khi_deque = deque(khi_sorted_list)
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

        filtered_index_list = [x for (y, x) in sorted(
            zip(peak_intensity_matrix, filtered_index_list),
            key=lambda pair: pair[0])][-4:]

        filtered_index_list = sort_index_list(filtered_index_list)
        return filtered_index_list

    @staticmethod
    def square(intensity_matrix, index_list, size_list, is_plot):
        logging.debug("size_list: {0}".format(size_list))
        square_instances = []
        logging.debug("The square size is {0}".format(size_list))

        for i in index_list:
            square_instance = Square(
                i, size_list, limitation=intensity_matrix.shape
            )
            square_instances.append(square_instance)

        peak = np.asanyarray(
            [i.sum(intensity_matrix) for i in square_instances]
        )
        [i.draw() for i in square_instances if is_plot]

        peak_intensity_matrix = peak[:, 0]
        peak_matrix_points_matrix = peak[:, 1]

        return (
            peak_intensity_matrix,
            peak_matrix_points_matrix,
            square_instances
        )

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0, **param):
        data_dict, plot_dict = Reader(self.file).matrix_data()

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

        bg_noise_float = np.average(
            (outer_peak_intensity_matrix - inner_peak_intensity_matrix) /
            (outer_peak_matrix_points_matrix - inner_peak_matrix_points_matrix)
        )
        peak_net_intensity_matrix = (
            inner_peak_intensity_matrix -
            bg_noise_float * inner_peak_matrix_points_matrix)

        self.display_save_image(plot_dict, is_save=is_save, is_show=0)
        self.save_config(plot_dict)

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

        thickness = float(result['dict']['db']['thickness'])

        logging.info("Sample thickness is {0}\n".format(thickness))

        eta = [self.correction(x[0], thickness) for x in index_list]

        logging.info("Intensity correction index eta is {0}".format(eta))

        coefficient_list = np.asarray([
            0.939691064,
            0.666790711 / (eta[1] / eta[0]),
            0.426843274 / (eta[2] / eta[0]),
            0.72278158 / (eta[3] / eta[0])
        ])
        beam_intensity_float = float(result['dict']['db']['beam_intensity'])
        peak_net_intensity_matrix = (
            peak_net_intensity_matrix *
            10000 *
            coefficient_list /
            beam_intensity_float
        )

        result['peak_intensity_matrix'] = peak_net_intensity_matrix
        return result

    @PrintLogDecorator()
    def print_result_csv(self):
        result = self.plot(is_save=1, is_show=0)
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
        sample.plot(is_show=0, is_save=1)
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
