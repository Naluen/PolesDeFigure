import abc
import configparser
import logging
import os
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure

import XrdAnalysis
from XrdAnalysis import Square
from XrdAnalysis.bruker3 import DatasetDiffractPlusV3

__author__ = 'Ang ZHOU (azhou@insa-rennes.fr)'
__project__ = 'XrdAnalysis'


class Material(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass


class GaP(Material):
    def __init__(self):
        super(GaP, self).__init__()
        self.hkl = {
            "0 0 2": 32,
            "0 0 4": 68,
            "0 0 6": 115,
            "-2 -2 4": 87
        }


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
        config = configparser.ConfigParser()
        config.add_section('db')
        [config.set('db', key, self.scan_dict[key]) for key
         in self.scan_dict.keys()]
        path = 'config.ini'
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
        data = np.asarray(
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

    @staticmethod
    def render(template_path, context):
        try:
            import jinja2
        except ImportError:
            jinja2 = None
            return ''
        else:
            path, filename = os.path.split(template_path)
            result_str = jinja2.Environment(
                loader=jinja2.FileSystemLoader(path or './')
            ).get_template(filename).render(context)
            return result_str

    def dp_sv_con_fig(self, label, is_show, is_save, **param):
        sample_name_str = self.scan_dict['sample']
        plt.title(sample_name_str + " " + label + "\n")
        self.save_config()

        if is_show:
            plt.show()

        if is_save:
            if "save_name" not in param:
                fig_file_name = os.path.join(
                    sample_name_str + '_' + label + '.png')
            else:
                fig_file_name = param['save_name']
            logging.info("Saving figure file to {0}".format(fig_file_name))
            plt.savefig(
                fig_file_name,
                dpi=200,
                bbox_inches='tight')
            return fig_file_name

    @abc.abstractmethod
    def plot(self):
        return NotImplementedError


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

    def plot(self):
        plt.plot(self.data_dict['two_theta_data'],
                 self.data_dict['int_data'])

        self.dp_sv_con_fig('Detector', is_save=0, is_show=1)


class TwoDScan(XrdScan):
    def __init__(self):
        super(TwoDScan, self).__init__()

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

    @PrintLogDecorator()
    def plot(self, is_save=1, is_show=0, **param):
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

        fig_file_name = self.dp_sv_con_fig('2D', is_show, is_save, **param)
        return fig_file_name

    @PrintLogDecorator()
    def polar_plot(self, is_save=1, is_show=0, **param):
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
        fig_file_name = self.dp_sv_con_fig('Polar', is_show, is_save, **param)

        return fig_file_name

    @PrintLogDecorator()
    def plot_square(self, is_save=1, is_show=0, is_plot=1, **param):
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
            i - k for (i, k) in
            zip(in_sq_instance_list,
                ot_sq_instance_list)
        ])

        fig_name = self.dp_sv_con_fig('Measurement', is_show, is_save, **param)

        logging.info("Net Peak intensity for each peak is {0}".format(
            peak_net_intensity_matrix
        ))
        logging.debug("The square size is {0}".format(size_list))

        result_dict = {
            'peak_intensity_matrix': peak_net_intensity_matrix,
            'index': index_list,
            'square_instances': square_instances,
            'fig_file_name': fig_name
        }

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
    def print_result_csv(self, peak_int_matrix):
        scan_dict = self.scan_dict.copy()

        sample_name = scan_dict['sample']

        mt_table_file = os.path.join(
            '{0}_result.csv'.format(sample_name)
        )
        import csv
        with open(mt_table_file, 'w') as tableTeX:
            spam_writer = csv.writer(tableTeX, dialect='excel')
            spam_writer.writerow(['MT-A', 'MT-D', 'MT-C', 'MT-B'])
            spam_writer.writerow(peak_int_matrix)

        return mt_table_file

    @PrintLogDecorator()
    def print_result_tex(self):
        plt.clf()

        class_path = os.path.dirname(XrdAnalysis.__file__)
        self.plot(is_save=0)
        res_dict = self.plot_square()
        fig_name = os.path.normcase(res_dict['fig_file_name'])
        peak_int_matrix = self.mt_intensity_to_fraction(res_dict)
        scan_dict = self.scan_dict.copy()
        rounded_int_list = np.around(peak_int_matrix, decimals=2).tolist()
        peak_int_list = (
            ['Intensity'] +
            [str(i) for i in rounded_int_list] +
            [str(np.round(sum(rounded_int_list)))]
        )

        tmp_file = os.path.join(class_path, 'templates', 'xrdtmp.tex')
        mt_dict = {'mt': peak_int_list,
                   'sr_int': float(scan_dict['beam_intensity']),
                   'mt_fig': fig_name}
        mt_v_fraction_data_str = self.render(tmp_file, mt_dict)

        return mt_v_fraction_data_str


class RsmScan(XrdScan):
    def __init__(self):
        super(RsmScan, self).__init__()

    @staticmethod
    def fill_array(array):
        array = np.asanyarray([i for i in array if i is not []])
        return array

    def stack_raw_data(self, data_list):
        phi_data_matrix = self.one_d_data(data_list, '_PHI')
        omega_data_matrix = self.one_d_data(data_list, '_OMEGA')
        int_data_matrix = self.two_d_data(data_list, 1)
        two_theta_data_matrix = self.two_d_data(data_list, 0)

        int_data_matrix = self.fill_array(int_data_matrix)
        two_theta_data_matrix = self.fill_array(two_theta_data_matrix)

        self.data_dict = {
            'int_data': int_data_matrix,
            'phi_data': phi_data_matrix,
            'omega_data': omega_data_matrix,
            'two_theta_data': two_theta_data_matrix
        }

        return self.data_dict

    @staticmethod
    def set_plt_style_plane():
        ax = plt.gca()

        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        ax.tick_params(
            axis="both", which="both", bottom="off", top="off",
            labelbottom="on",
            left="off", right="off", labelleft="on")

        xlabel = r'$S_x (nm^{-1})$'
        ylabel = r'$S_z (nm^{-1})$'
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        # ax.set_title('RSM', fontsize=14)

        ax.set_aspect('auto')

    def plot(self, is_save=1, is_show=0, **param):
        tth = self.data_dict['two_theta_data']
        phi = self.data_dict['phi_data'][0]
        omega = self.data_dict['omega_data']
        int_data = self.data_dict['int_data']
        wavelength_dict = {'Cu-1': 0.154055911278}
        material_class_dict = {'GaP': GaP}

        lam = wavelength_dict['Cu-1']
        mt_ins = material_class_dict['GaP']()
        hkl = [i for i in list(mt_ins.hkl.keys())
               if mt_ins.hkl[i] - 3 <= tth[0][0] <= mt_ins.hkl[i] + 3]
        if len(hkl) is not 1:
            logging.error('HKL Value Error')
            hkl = "0 0 0"
        else:
            hkl = hkl[0]
            self.scan_dict['hkl'] = hkl

        s_mod = 2. / lam * np.sin(np.radians(tth / 2.))
        psi = omega - tth / 2.
        s_x = s_mod * np.sin(np.radians(psi))
        s_z = s_mod * np.cos(np.radians(psi))
        w, l = int_data.shape

        if (phi > -2) and (phi < 2):
            s_x = -s_x

        fig, ax = plt.subplots()

        xi = np.linspace(s_x.min(), s_x.max(), w)
        yi = np.linspace(s_z.min(), s_z.max(), l)
        xi, yi = np.meshgrid(xi, yi)
        zi = griddata((s_x.flatten(), s_z.flatten()), int_data.flatten(),
                      (xi, yi), method='linear')

        im = plt.imshow(
            zi,
            origin='lower',
            norm=LogNorm(vmin=int_data.min() + 1, vmax=int_data.max()),
            extent=[s_x.min(), s_x.max(), s_z.min(), s_z.max()]
        )

        self.set_plt_style_plane()

        cb = fig.colorbar(im, ax=ax, extend='max')
        cb.set_label(
            r'$Intensity\ (Counts\ per\ second)$',
            fontsize=14)

        label = hkl + ' RSM'
        fig_name = self.dp_sv_con_fig(label, is_show, is_save, **param)

        return fig_name

    def print_result_tex(self):
        class_path = os.path.dirname(XrdAnalysis.__file__)
        fig_name = self.plot()
        fig_name = os.path.normcase(fig_name)
        scan_dict = self.scan_dict.copy()

        tmp_file = os.path.join(class_path, 'templates', 'rsmtmp.tex')
        rsm_dict = {'fig_name': scan_dict['hkl'],
                    'rsm_fig': fig_name}
        txt_str = self.render(tmp_file, rsm_dict)

        return txt_str


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
    def read_config_file(raw_file, is_ignore_intensity=0):
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
                return 'sample'
            else:
                if sample_name:
                    sample_name = sample_name[-1]
                    return sample_name
                else:
                    return 'sample'

        def beam_intensity():
            logging.info("Starts calculation beam intensity.")

            source_file_list = ['beam1mm.raw', 'beam8mm.raw']
            beam_int_list = [
                (RawFile(i).read_data().get_max() * 8940)
                for i in source_file_list if os.path.isfile(i)]

            if not beam_int_list:
                logging.warning("Could not found default beam intensity files")
                beam_df_int = 6000 * 8940
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
            os.path.join(os.path.dirname(XrdAnalysis.__file__), 'config.ini'),
            os.path.join(os.path.dirname(raw_file), 'config.ini')
        ])

        logging.debug(dict(config['db']))

        config_db = config['db']
        config_db['sample'] = str(
            guess_sample_name(raw_file) or
            config_db['sample']
        )
        os.chdir(os.path.dirname(raw_file))
        if not is_ignore_intensity:
            config_db['beam_intensity'] = str(
                config_db['beam_intensity'] or
                beam_intensity()
            )

        logging.debug(dict(config_db))

        return config_db

    def read_data(self):
        data_list = self.__str_data(self.file)
        scan_instance = self.create_scan_instance()
        instance_dict = self.get_scan_dict()
        if instance_dict['_TYPE'] == 'TwoDPlot':
            instance_dict.update(
                self.read_config_file(self.file)
            )
        else:
            instance_dict.update(
                self.read_config_file(self.file, is_ignore_intensity=1)
            )

        scan_instance.stack_raw_data(data_list)
        scan_instance.scan_dict = instance_dict

        logging.debug("Scan dict: {0}".format(instance_dict))

        os.chdir(os.path.dirname(self.file))

        return scan_instance


class H5File(XrdFile):
    def __init__(self, h5_file):
        if isinstance(h5_file, str):
            h5_file = [h5_file, '']
            logging.warning("Sub file name is required by access_file")

        if not h5_file[0].endswith('.h5'):
            raise TypeError("Illegal h5 file name.")

        super(H5File, self).__init__(h5_file)

        self.file_handle, self.h5_file_handle = self.access_file()

    def access_file(self, **param):
        h5_file_name = self.file[0]
        if 'sub_file_name' not in param:
            sub_file_name = self.file[1]
        else:
            sub_file_name = param['sub_file_name']

        h5_file_handle = h5py.File(h5_file_name, 'a')
        try:
            file_handle = h5_file_handle[sub_file_name]
        except KeyError:
            logging.error("No such sub file.")

        if file_handle is None:
            logging.warning("No such file in lib.")

        return file_handle, h5_file_handle

    def get_scan_dict(self):
        scan_dict = {i: self.file_handle.attrs[i]
                     for i in self.file_handle.attrs.keys()}

        return scan_dict

    def get_scan_type(self):
        scan_type = self.file_handle.attrs.get('_TYPE')

        return scan_type

    def read_raw(self, raw_file):
        # Require data.
        instance = RawFile(raw_file).read_data()
        sub_file_name = (
            instance.scan_dict['sample'] +
            '/' + os.path.basename(raw_file).split('.')[0])
        self.h5_file_handle.close()
        file_handle = self.h5_file_handle.require_group(sub_file_name)

        # Record data.
        data_dict = instance.data_dict.copy()
        for key in data_dict.keys():
            try:
                del file_handle[key]
            except (TypeError, KeyError):
                pass
            file_handle.create_dataset(
                key,
                data=data_dict[key]
            )
        scan_dict = instance.scan_dict.copy()
        for key in scan_dict.keys():
            try:
                file_handle.attrs.modify(key, scan_dict[key])
            except TypeError:
                pass

                # Close data set.

    def read_data(self):
        scan_instance = self.create_scan_instance()
        data_set = self.file_handle.items()

        scan_instance.data_dict = {i[0]: np.asarray(i[1])for i in data_set}
        scan_instance.scan_dict = self.get_scan_dict()

        logging.debug("Scan dict: {0}".format(scan_instance.scan_dict))

        return scan_instance

    def __del__(self):
        self.h5_file_handle.close()


def reader(file):
    if isinstance(file, str) and file.endswith('.raw'):
        return RawFile(file)
    else:
        return H5File(file)
