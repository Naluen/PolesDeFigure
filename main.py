import logging
import os
import sys

import matplotlib.pyplot as plt
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)

from Poles import TwoDFigureRaw, MeasureFigureRaw
from gui.main_window import UiMainWindow

try:
    import configparser
except ImportError:
    import ConfigParser as configparser


class PfGui(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = UiMainWindow(self)
        self.ui.push_button.setEnabled(False)

        self.results = {
            'peak_intensity_matrix': None,
            'index': None,
            'square_instances': None
        }
        self.sample = None
        self.press = None
        self.selected_square = []
        self.cid_press = None
        self.cid_release = None
        self.cid_motion = None

        self.ui.actionOpen_2.triggered.connect(self.open_file)
        self.ui.actionSave.triggered.connect(self.save_image)

        self.figure = plt.figure()
        self.figure.set_size_inches((30, 6))
        self.canvas = FigureCanvas(self.figure)
        self.ui.verticalLayout.addWidget(self.canvas)
        self.connect_canvas()

        self.ui.push_button.clicked.connect(self.automatically_detect_peak)

    def open_file(self):

        self.results = {
            'peak_intensity_matrix': None,
            'index': None,
            'square_instances': None
        }

        config = configparser.ConfigParser()
        config_file_str = os.path.join(
            os.path.dirname(sys.argv[0]),
            'config.ini'
        )
        config.read(
            config_file_str
        )

        if os.path.exists(config.get('db', 'directory')):
            dir_item = config.get('db', 'directory')
        else:
            dir_item = '\home'

        raw_file_name = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file',
            dir_item,
            "RAW files (*.raw)")
        raw_file_name = str(raw_file_name[0])

        if raw_file_name:
            sample = TwoDFigureRaw(raw_file_name)
            self.figure.clear()
            plt.gca()

            sample.plot()
            self.canvas.draw()

            self.ui.push_button.setEnabled(True)

            self.sample = raw_file_name

        else:
            pass

    def save_image(self):
        config = configparser.ConfigParser()
        config_file_str = os.path.join(
            os.path.dirname(sys.argv[0]),
            'config.ini'
        )
        config.read(
            config_file_str
        )

        if os.path.exists(config.get('db', 'directory')):
            dir_item = config.get('db', 'directory')
        else:
            dir_item = '\home'

        supported_file_dict = self.canvas.get_supported_filetypes()
        supported_file_list = [
            "{0} (*.{1})".format(supported_file_dict[i], i)
            for i in supported_file_dict
            ]
        supported_file_str = ";;".join(supported_file_list)
        image_file_name = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Open file',
            dir_item,
            supported_file_str)
        image_file_name = str(image_file_name[0])

        self.figure.savefig(image_file_name, bbox_inches='tight')

        config.set('db', 'directory', os.path.dirname(image_file_name))
        with open(config_file_str, 'w') as config_file:
            config.write(config_file)

    def automatically_detect_peak(self):
        sample = MeasureFigureRaw(self.sample)
        self.instance = sample

        if self.results['square_instances'] is not None:
            for i in self.results['square_instances']:
                i.remove()

        results = sample.plot(is_show=1)
        self.results = results
        self.canvas.draw()

        results = sample.mt_intensity_to_fraction(results)
        text_label_list = [
            self.ui.mt_a_value,
            self.ui.mt_d_value,
            self.ui.mt_c_value,
            self.ui.mt_b_value
        ]
        for (i, j) in zip(text_label_list, results['peak_intensity_matrix']):
            i.setText(str(j))

    def connect_canvas(self):
        self.cid_press = self.canvas.mpl_connect(
            'button_press_event', self.on_press
        )
        self.cid_release = self.canvas.mpl_connect(
            'button_release_event', self.on_release
        )
        self.cid_motion = self.canvas.mpl_connect(
            'motion_notify_event', self.on_motion
        )

    def on_press(self, event):
        ax_list = self.figure.axes
        if self.results['square_instances'] is None:
            return
        else:
            for square in self.results['square_instances']:
                if event.inaxes != ax_list[0]:
                    pass
                else:
                    if [event.xdata, event.ydata] in square:
                        [x0, y0] = square.central_point_list
                        self.selected_square.append(square)
                        self.press = x0, y0, event.xdata, event.ydata

    def on_motion(self, event):
        if self.results['square_instances'] is None:
            return
        if self.results['square_instances'] is not []:
            for square in self.selected_square:
                if event.inaxes != self.figure.axes[0]:
                    pass
                else:
                    x0, y0, init_mouse_x, init_mouse_y = self.press
                    dy = event.xdata - init_mouse_x
                    dx = event.ydata - init_mouse_y
                    square.remove()
                    square.central_point_list = [x0 + dx, y0 + dy]
                    square.draw()

        self.canvas.draw()

    def on_release(self, event):
        if self.results['square_instances'] is None:
            return
        self.press = None
        self.selected_square = []
        self.canvas.draw()

        outer_index_list = [
            i.central_point_list for i in self.results['square_instances'][:4]
            ]

        results = self.instance.plot(
            is_save_image=0,
            is_calculation=1,
            outer_index_list=outer_index_list
        )
        results = self.instance.mt_intensity_to_fraction(results)
        text_label_list = [
            self.ui.mt_a_value,
            self.ui.mt_d_value,
            self.ui.mt_c_value,
            self.ui.mt_b_value
        ]
        for (i, j) in zip(text_label_list, results['peak_intensity_matrix']):
            i.setText(str(j))

    def disconnect(self):
        self.canvas.mpl_disconnect(self.cid_press)
        self.canvas.mpl_disconnect(self.cid_release)
        self.canvas.mpl_disconnect(self.cid_motion)


if __name__ == '__main__':
    logging.basicConfig(
        # filename=os.path.join(
        #     os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )

    Program = QtWidgets.QApplication(sys.argv)
    MyProgram = PfGui()
    MyProgram.show()
    Program.aboutToQuit.connect(MyProgram.disconnect)
    sys.exit(Program.exec_())
