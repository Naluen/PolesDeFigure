import logging

import matplotlib.pyplot as plt
import numpy as np


class Square(object):
    def __init__(
            self,
            cr_list,
            size_tuple,
            int_matrix,
            limitation=([0, 70], [-30, 330])
    ):
        self.cr_list = cr_list
        self.size_tuple = size_tuple
        self.int_matrix = int_matrix
        self.limitation = limitation
        self.figure_handle = None

    def lim(self):
        (y_limit, x_limit) = self.limitation
        x_min = max(0, int(
            np.floor(self.cr_list[1] - self.size_tuple[1] / 2.0)))
        x_max = min(
            int(np.floor(self.cr_list[1] + self.size_tuple[1] / 2.0)),
            int(x_limit[1] - x_limit[0]))
        y_min = max(0, int(
            np.floor(self.cr_list[0] - self.size_tuple[0] / 2.0)))
        y_max = min(
            int(np.floor(self.cr_list[0] + self.size_tuple[0] / 2.0)),
            int(y_limit[1] - y_limit[0]))

        return x_min, x_max, y_min, y_max

    def draw(self):
        x_min, x_max, y_min, y_max = self.lim()
        (y_limit, x_limit) = self.limitation
        x_list = np.asarray([x_min, x_max, x_max, x_min, x_min]) + x_limit[0]
        y_list = np.asarray([y_min, y_min, y_max, y_max, y_min]) + y_limit[0]
        self.figure_handle, = plt.plot(x_list, y_list, linewidth=0.5)

    def move(self, direct_tuple):
        self.remove()
        self.cr_list = [i + j for (i, j) in
                        zip(self.cr_list, list(direct_tuple))]
        self.draw()

    def sum(self):
        x_min, x_max, y_min, y_max = self.lim()

        intensity_result_matrix = self.int_matrix[y_min:y_max, x_min:x_max]
        peak_intensity_int = np.sum(intensity_result_matrix)
        b, p = intensity_result_matrix.shape
        peak_matrix_points = b * p

        return peak_intensity_int, peak_matrix_points

    def remove(self):
        """
        Remove the the square _lines.
        :return: None
        """
        axes = plt.gca()
        try:
            axes.lines.remove(self.figure_handle)
        except ValueError:
            logging.warning("Could not find the square.")

    def __sub__(self, x):
        if isinstance(x, Square):
            pk_int, pk_pt = self.sum()
            x_pk_int, x_pk_pt = x.sum()
            bg_noise_float = (pk_int - x_pk_int) / (pk_pt - x_pk_pt)

            return pk_int - pk_pt * bg_noise_float
        else:
            return NotImplemented

    def __contains__(self, item):
        """
        Check if the point is in the square.
        :param item: the position of point [x,y].
        :return: The boolean value.
        """
        x_min, x_max, y_min, y_max = self.lim()
        (y_limit, x_limit) = self.limitation
        if (
                                x_min + x_limit[0] < item[0] < x_max + x_limit[
                        0] and
                                y_min + y_limit[0] < item[1] < y_max + y_limit[
                        0]):
            return True
        else:
            return False
