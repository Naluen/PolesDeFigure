import logging
import os
import sys

import matplotlib.pyplot as plt
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)

from XrdAnalysis import Reader
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
        self.instance = Reader.TwoDScan()
        self.selected_square = []
        self.cid_press = None
        self.cid_release = None
        self.cid_motion = None

        os.chdir(os.path.dirname(sys.argv[0]))
        self.gui_cfg = configparser.ConfigParser()
        self.gui_cfg.read('GUI.cfg')

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

        dir_item = self.gui_cfg.get('DEFAULT', 'directory') or '\home'

        raw_file_name = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Open file',
            dir_item,
            "RAW files (*.raw)")
        raw_file_name = str(raw_file_name[0])

        if raw_file_name:
            raw_instance = Reader.RawFile(raw_file_name).read_data()
            self.figure.clear()
            raw_instance.plot(is_show=0, is_save=0)
            self.canvas.draw()
            self.ui.push_button.setEnabled(True)
            self.instance = raw_instance
            self.gui_cfg.set('DEFAULT', 'directory',
                             os.path.dirname(raw_file_name))
        else:
            pass

    def save_image(self):
        supported_file_dict = self.canvas.get_supported_filetypes()
        supported_file_list = [
            "{0} (*.{1})".format(supported_file_dict[i], i)
            for i in supported_file_dict
            ]
        supported_file_str = ";;".join(supported_file_list)

        dir_item = self.gui_cfg.get('DEFAULT', 'directory') or '\home'

        image_file_name = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Open file',
            dir_item,
            supported_file_str)
        image_file_name = str(image_file_name[0])

        if image_file_name:
            self.figure.savefig(image_file_name, bbox_inches='tight')
        else:
            return

    def save_config(self):
        script_path = os.path.dirname(sys.argv[0])
        with open(os.path.join(script_path, 'GUI.cfg'), 'w') as config_file:
            self.gui_cfg.write(config_file)

    def automatically_detect_peak(self):
        if self.results['square_instances']:
            for i in self.results['square_instances']:
                i.remove()

        results = self.instance.plot_square(is_show=0, is_save=0)
        self.results = results
        self.canvas.draw()

        results = self.instance.mt_intensity_to_fraction(results)
        self.update_fraction_list(results)

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
        if event.inaxes != ax_list[0]:
            return
        for square in self.results['square_instances']:
            if [event.xdata, event.ydata] in square:
                self.selected_square.append(square)

    def on_motion(self, event):
        if not self.selected_square:
            return
        if event.inaxes != self.figure.axes[0]:
            return
        [init_mouse_y, init_mouse_x] = self.selected_square[0].cr_list
        for square in self.selected_square:
            dy = event.xdata - init_mouse_x + 30
            dx = event.ydata - init_mouse_y
            square.move((dx, dy))

        self.canvas.draw()

    def update_fraction_list(self, results):
        text_label_list = [
            self.ui.mt_a_value,
            self.ui.mt_d_value,
            self.ui.mt_c_value,
            self.ui.mt_b_value
        ]
        for (i, j) in zip(text_label_list, results):
            i.setText(str(j))

    def on_release(self, event):
        if self.results['square_instances'] is None:
            return
        self.selected_square = []

        outer_index_list = [
            i.cr_list for i in self.results['square_instances'][:4]
            ]
        results_dict = self.instance.plot_square(
            is_show=0,
            is_save=0,
            is_plot=0,
            outer_index_list=outer_index_list
        )
        results = self.instance.mt_intensity_to_fraction(results_dict)
        self.instance.print_result_csv(results)
        self.update_fraction_list(results)

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
    Program.aboutToQuit.connect(MyProgram.save_config)
    sys.exit(Program.exec_())
