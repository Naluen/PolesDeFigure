import logging

import matplotlib.pyplot as plt
import numpy as np


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
