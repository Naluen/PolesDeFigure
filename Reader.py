import abc
import configparser
import logging
import os
import re
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure

import Square
from bruker3 import DatasetDiffractPlusV3

__author__ = 'Ang ZHOU (azhou@insa-rennes.fr)'
__project__ = 'XrdAnalysis'


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


class XrdScan(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.data_dict = {}
        self.scan_dict = {}

    def save_config(self):
        directory = self.scan_dict['directory']
        config = configparser.ConfigParser()
        config.add_section('db')
        [config.set('db', key, self.scan_dict[key]) for key
         in self.scan_dict.keys()]
        path = os.path.join(directory, 'config.ini')
        with open(path, 'w') as configfile:
            config.write(configfile)
        del config

    @staticmethod
    def two_d_data(data_list, index):
        data = np.asanyarray(
            [
                [float(j.split('\t')[index]) for j in
                 value
                 if (re.match('\d', j) or j.startswith("-"))]
                for value in data_list[1:]
                ]
        )
        return data

    @staticmethod
    def one_d_data(data_list, key_word):
        data = np.asanyarray(
            [
                [float(j.split('=')[1]) for j
                 in value if j.startswith(key_word)]
                for value in data_list[1:]
                ]
        )
        return data

    @abc.abstractmethod
    def stack_raw_data(cls, data_list):
        pass


class DetectorScan(XrdScan):
    def __init__(self):
        super(DetectorScan, self).__init__()

    def stack_raw_data(self, data_list):
        int_data_matrix = self.two_d_data(data_list, 1)
        two_theta_data_matrix = self.two_d_data(data_list, 0)
        self.data_dict = {'int_data': int_data_matrix,
                          'two_theta_data': two_theta_data_matrix
                          }
        return self.data_dict

    def get_max(self):
        return np.max(self.data_dict['int_data'])


class TwoDScan(XrdScan):
    def __init__(self):
        super(TwoDScan, self).__init__()

    def stack_raw_data(self, data_list):
        phi_data_matrix = self.two_d_data(data_list, 0)
        int_data_matrix = (self.two_d_data(data_list, 1) /
                           self.one_d_data(data_list, '_STEPTIME'))
        khi_data_matrix = self.one_d_data(data_list, '_KHI')
        self.data_dict = {
            'int_data': int_data_matrix,
            'phi_data': phi_data_matrix,
            'khi_data': khi_data_matrix
        }

        return self.data_dict

    def display_save_image(self, label, is_show, is_save):
        sample_name_str = self.scan_dict['sample']
        plt.title(sample_name_str + " " + label + "\n")

        if is_show:
            plt.show()

        figure_file_name = os.path.join(
            self.scan_dict['directory'],
            sample_name_str + '_' + label + '.png'
        )

        if is_save:
            logging.info(
                "Saving figure file to {0}".format(figure_file_name))
            plt.savefig(
                figure_file_name,
                dpi=200,
                bbox_inches='tight')

        self.save_config()

    # @staticmethod
    # def add_color_bar(im, aspect=20, pad_fraction=0.5, **kwargs):
    #     """Add a vertical color bar to an image plot."""
    #     from mpl_toolkits import axes_grid1
    #
    #     divider = axes_grid1.make_axes_locatable(im.axes)
    #     width = axes_grid1.axes_size.AxesY(im.axes, aspect=1 / aspect)
    #     pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    #     current_ax = plt.gca()
    #     cax = divider.append_axes("right", size=width, pad=pad)
    #     plt.sca(current_ax)
    #     return im.axes.figure.colorbar(im, cax=cax, **kwargs)

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0):
        plot_data_dict = self.data_dict.copy()
        scan_dict = self.scan_dict.copy()
        v_max = int(scan_dict['v_max'])
        v_min = int(scan_dict['v_min'])
        ver_min = plot_data_dict['phi_data'][0, :][0]
        ver_max = plot_data_dict['phi_data'][0, :][-1]
        hor_min = plot_data_dict['khi_data'][0]
        hor_max = plot_data_dict['khi_data'][-1]

        ax2d = plt.gca()
        im = ax2d.imshow(
            plot_data_dict['int_data'],
            origin="lower",
            norm=LogNorm(vmin=v_min, vmax=v_max),
            extent=[ver_min, ver_max, hor_min, hor_max]
        )
        ax2d.tick_params(axis='both', which='major', labelsize=10)

        ticks = np.logspace(1, np.log10(v_max), np.log10(v_max))
        plt.colorbar(im, fraction=0.012, pad=0.04,
                     format="%.e", extend='max', ticks=ticks)

        self.display_save_image('2D', is_show, is_save)

    @PrintLogDecorator()
    def polar_plot(self, is_save=1, is_show=0):
        """
        Plot polar figure
        :param is_save: if save the image.
        :param is_show: if show the image.
        :return:
        """
        scan_dict = self.scan_dict.copy()
        plot_data_dict = self.data_dict.copy()
        plot_data_dict['phi_data'] = np.radians(
            plot_data_dict['phi_data'] +
            [float(scan_dict['phi_offset'])]
        )

        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.tick_params(axis="y", labelsize=15, labelcolor="white")
        ax.tick_params(axis="x", labelsize=15)
        ax.set_rmin(0.)

        [xx, yy] = np.meshgrid(plot_data_dict['phi_data'][0, :].T,
                               plot_data_dict['khi_data'])
        v_max = int(scan_dict['v_max'])
        v_min = int(scan_dict['v_min'])
        im = ax.pcolormesh(xx, yy, plot_data_dict['int_data'],
                           norm=LogNorm(vmin=v_min, vmax=v_max))
        ax.grid(color="white")

        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=16)

        ticks = np.logspace(0, np.log10(v_max), np.log10(v_max) + 1)
        plt.colorbar(im, pad=0.1, format="%.e", extend='max', ticks=ticks)
        self.display_save_image('Polar', is_show, is_save)

    @staticmethod
    def correction(chi, thickness):
        theta = 14.22  # for GaP
        omega = 14.22  # for GaP
        thickness /= 10000.
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
    def peak_search(cls, bg_data_dict):
        def sort_index_list(index_list):
            """
            Sort index list to fit ABCD micro-Twins, where chi of A is max.
            :param index_list: list for each point with form [chi, khi]
            :return: sorted list.
            """
            khi_sorted_list = sorted(index_list, key=lambda pair: pair[1])
            chi_index_list = [l[0] for l in khi_sorted_list]
            shifted_index_int = chi_index_list.index(max(chi_index_list))
            from collections import deque
            khi_deque = deque(khi_sorted_list)
            khi_deque.rotate(-shifted_index_int)
            sorted_index_list = list(khi_deque)

            logging.debug("index list before sort:{0}".format(index_list))
            logging.info(
                "index list after sort:{0}".format(sorted_index_list)
            )
            return sorted_index_list

        int_data_matrix = bg_data_dict['int_data']
        ver_min = bg_data_dict['phi_data'][0, :][0]
        ver_max = bg_data_dict['phi_data'][0, :][-1]
        hor_min = bg_data_dict['khi_data'][0]
        hor_max = bg_data_dict['khi_data'][-1]

        neighborhood = generate_binary_structure(2, 2)
        for i in range(3):
            int_data_matrix = gaussian_filter(int_data_matrix, 4,
                                              mode='nearest')
        local_max = (
            maximum_filter(int_data_matrix,
                           footprint=neighborhood) == int_data_matrix
        )
        index = np.asarray(np.where(local_max))
        ft_index_list = [
            [i, j] for (i, j) in zip(index[0, :], index[1, :])
            ]

        chi_threshold = 40
        ft_index_list = [
            i for i in ft_index_list if i[0] < chi_threshold
            ]

        in_sq_instance_list = [
            Square.Square(
                i, [10, 10], int_data_matrix,
                limitation=([hor_min, hor_max], [ver_min, ver_max])
            ) for i in ft_index_list
            ]
        ot_sq_instance_list = [
            Square.Square(
                i, [20, 20], int_data_matrix,
                limitation=([hor_min, hor_max], [ver_min, ver_max])
            ) for i in ft_index_list
            ]
        int_list = [i - k for (i, k)
                    in zip(in_sq_instance_list, ot_sq_instance_list)]
        ft_index_list = [x for (y, x) in sorted(
            zip(int_list, ft_index_list),
            key=lambda pair: pair[0])][-4:]
        ft_index_list = sort_index_list(ft_index_list)

        return ft_index_list

    @PrintLogDecorator()
    def plot_square(self, is_save=1, is_show=0, is_plot=1,**param):
        int_data_matrix = self.data_dict['int_data']
        scan_dict = self.scan_dict.copy()

        if 'outer_index_list' in param:
            index_list = param['outer_index_list']
        else:
            index_list = self.peak_search(self.data_dict)
        size_list = list(map(int, scan_dict['square_size'].split(',')))

        ver_min = self.data_dict['phi_data'][0, :][0]
        ver_max = self.data_dict['phi_data'][0, :][-1]
        hor_min = self.data_dict['khi_data'][0]
        hor_max = self.data_dict['khi_data'][-1]
        in_sq_instance_list = [
            Square.Square(
                i, size_list, int_data_matrix,
                limitation=([hor_min, hor_max], [ver_min, ver_max]))
            for i in index_list
            ]
        ot_sq_instance_list = [
            Square.Square(
                i, [i + 4 for i in size_list], int_data_matrix,
                limitation=([hor_min, hor_max], [ver_min, ver_max]))
            for i in index_list
            ]
        square_instances = in_sq_instance_list + ot_sq_instance_list
        [i.draw() for i in square_instances if is_plot]

        peak_net_intensity_matrix = np.asarray([
            i - k for (i, k) in zip(in_sq_instance_list, ot_sq_instance_list)
            ])

        self.display_save_image('Measurement', is_show, is_save)

        logging.info("Net Peak intensity for each peak is {0}".format(
            peak_net_intensity_matrix
        ))
        logging.debug("The square size is {0}".format(size_list))

        result_dict = {
            'peak_intensity_matrix': peak_net_intensity_matrix,
            'index': index_list,
            'square_instances': square_instances}

        return result_dict

    def mt_intensity_to_fraction(self, result):
        scan_dict = self.scan_dict.copy()
        index_list = result['index']
        peak_net_intensity_matrix = result['peak_intensity_matrix']
        thickness = float(scan_dict['thickness'] or 900)
        beam_intensity_float = float(scan_dict['beam_intensity'])

        eta = [self.correction(x[0], thickness) for x in index_list]

        logging.info("Intensity correction index eta is {0}".format(eta))

        coefficient_list = np.asarray([
            0.939691064,
            0.666790711 / (eta[1] / eta[0]),
            0.426843274 / (eta[2] / eta[0]),
            0.72278158 / (eta[3] / eta[0])
        ])
        peak_net_intensity_matrix = (
            peak_net_intensity_matrix *
            10000 *
            coefficient_list /
            beam_intensity_float
        )

        logging.info("Sample thickness is {0}\n".format(thickness))
        logging.info("Peak intensity is {0}".format(peak_net_intensity_matrix))

        return peak_net_intensity_matrix

    @PrintLogDecorator()
    def print_result_csv(self, peak_intensity_matrix):
        scan_dict = self.scan_dict.copy()

        directory = scan_dict['directory']
        sample_name = scan_dict['sample']

        mt_table_file = os.path.join(
            directory,
            '{0}_result.csv'.format(sample_name)
        )
        import csv
        with open(mt_table_file, 'w') as tableTeX:
            spam_writer = csv.writer(tableTeX, dialect='excel')
            spam_writer.writerow(['MT-A', 'MT-D', 'MT-C', 'MT-B'])
            spam_writer.writerow(peak_intensity_matrix)

        return  peak_intensity_matrix


class RsmScan(XrdScan):
    def __init__(self):
        super(RsmScan, self).__init__()

    def stack_raw_data(self, data_list):
        phi_data_matrix = self.one_d_data(data_list, '_PHI')
        omega_data_matrix = self.one_d_data(data_list, '_OMEGA')
        int_data_matrix = self.two_d_data(data_list, 1)
        two_theta_data_matrix = self.two_d_data(data_list, 0)
        self.data_dict = {
            'int_data': int_data_matrix,
            'phi_data': phi_data_matrix,
            'omega_data': omega_data_matrix,
            'two_theta_data': two_theta_data_matrix
        }

        return self.data_dict


class XrdFile(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, file):
        self.file = file

    @abc.abstractmethod
    def read_data(self):
        pass

    @abc.abstractmethod
    def get_scan_type(self):
        pass

    def create_scan_instance(self):
        scan_type = str(self.get_scan_type())
        scan_type_dict = {
            'SingleScanPlot': DetectorScan,
            'TwoDPlot': TwoDScan,
            'RSMPlot': RsmScan
        }
        scan_instance = scan_type_dict[scan_type]()

        logging.info("Scan Type is {0}".format(scan_type))

        return scan_instance


class RawFile(XrdFile):
    def __init__(self, raw_file):
        if not raw_file.endswith('.raw'):
            raise TypeError("Illegal raw file name.")
        super(RawFile, self).__init__(raw_file)

        logging.info("Reading raw file {0}".format(raw_file))

    def get_scan_dict(self):
        data_list = self.__str_data(self.file)
        scan_dict = {i.split('=')[0].strip(): i.split('=')[1].strip()
                     for i in data_list[0] if i.startswith("_")}
        return scan_dict

    def get_scan_type(self):
        scan_dict = self.get_scan_dict()
        scan_type = scan_dict['_TYPE']
        return scan_type

    @staticmethod
    def __str_data(raw_file):
        with open(raw_file, 'rb') as fp:
            ds = DatasetDiffractPlusV3(fp)
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

        def beam_intensity(directory):
            logging.info("Starts calculation beam intensity.")

            source_file_list = [
                os.path.join(directory, 'beam1mm.raw'),
                os.path.join(directory, 'beam8mm.raw')]
            beam_int_list = [
                (RawFile(i).read_data().get_max() * 8940)
                for i in source_file_list if os.path.isfile(i)]

            if not beam_int_list:
                logging.warning("Could not found default beam intensity files")
                beam_df_int = 6000*8940
                beam_int_list = [beam_df_int]
                logging.info(
                    "Use default Source Beam Intensity = {0}.".format(
                        beam_int_list[0])
                )
            beam_int_list = np.asarray(beam_int_list)
            beam_int = np.mean(beam_int_list)

            logging.info("Beam Int is {0}".format(beam_int))

            return beam_int

        config = configparser.ConfigParser()
        config.read([
            os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
            os.path.join(os.path.dirname(raw_file), 'config.ini')
        ])

        logging.debug(dict(config['db']))

        config_db = config['db']
        config_db['sample'] = str(
            config_db['sample'] or
            guess_sample_name(raw_file)
        )
        config_db['directory'] = str(
            config_db['directory'] or
            os.path.dirname(raw_file)
        )
        config_db['raw_file'] = str(
            config_db['raw_file'] or
            os.path.basename(raw_file)
        )
        config_db['beam_intensity'] = str(
            config_db['beam_intensity'] or
            beam_intensity(config_db['directory'])
        )

        logging.debug(dict(config_db))

        return config_db

    def read_data(self):
        data_list = self.__str_data(self.file)
        scan_instance = self.create_scan_instance()
        scan_instance.stack_raw_data(data_list)
        instance_dict = self.get_scan_dict()
        if scan_instance.__class__.__name__ == 'TwoDScan':
            instance_dict.update(self.read_config_file(self.file))
        scan_instance.scan_dict = instance_dict

        logging.debug("Scan dict: {0}".format(instance_dict))

        return scan_instance


class H5File(XrdFile):
    def __init__(self, h5_file):
        if not h5_file.endswith('.h5'):
            raise TypeError("Illegal h5 file name.")
        super(H5File, self).__init__(h5_file)
        file_handle = h5py.File(self.file, 'a')

    def get_scan_type(self):
        scan_type = self.file[1].split('/')[-1]
        logging.info("Scan type is {0}".format(scan_type))
        return scan_type

    def create_file(self):
        with open(self.file, 'a'):
            file_handle = h5py.File(self.file)
        file_handle.attrs['HDF5_Version'] = h5py.version.hdf5_version
        file_handle.attrs['h5py_version'] = h5py.version.version

    def read_raw(self, raw_file):
        # Require data.
        instance = RawFile(raw_file)
        raw_data_dict, raw_plot_dict = instance.matrix_data()
        scan_dict = instance.get_head()
        # Require data set.
        file_handle = h5py.File(self.file, 'a')
        sample_name_str = raw_plot_dict['db']['sample']
        sample_handle = file_handle.require_group(sample_name_str)
        lr_handle = sample_handle.require_group(scan_dict['_TYPE'])
        # Record data.
        for key, value in raw_data_dict.items():
            if key in lr_handle:
                del lr_handle[key]

            lr_handle.create_dataset(
                key,
                data=value
            )
        for key, value in scan_dict.items():
            try:
                lr_handle.attrs.modify(key, value)
            except TypeError:
                pass
        if raw_plot_dict:
            for key, value in raw_plot_dict['db'].items():
                try:
                    lr_handle.attrs.modify(key, value)
                except TypeError:
                    pass
        # Close data set.
        file_handle.close()

    def matrix_data(self):
        pass


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
        sample = RawFile()
        data_dict, plot_dict = sample.matrix_data(raw_file_name)
        print(data_dict)
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
