# -*- coding: utf-8 -*-

# Created by: PyQt5 UI code generator 5.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class UiMainWindow(object):
    def __init__(self, main_window):
        main_window.setObjectName("main_window")
        main_window.resize(1188, 371)
        self.central_widget = QtWidgets.QWidget(main_window)
        self.central_widget.setObjectName("central_widget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.central_widget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.line = QtWidgets.QFrame(self.central_widget)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.horizontalLayout.addWidget(self.line)

        self.control_layout = QtWidgets.QVBoxLayout()
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()

        self.push_button = QtWidgets.QPushButton(self.central_widget)
        self.horizontalLayout_2.addWidget(self.push_button)

        spacerItem = QtWidgets.QSpacerItem(20, 40,
                                           QtWidgets.QSizePolicy.Minimum,
                                           QtWidgets.QSizePolicy.Expanding)
        self.horizontalLayout_2.addItem(spacerItem)
        self.horizontalLayout_2.addWidget(self.push_button)
        self.control_layout.addLayout(self.horizontalLayout_2)
        spacerItem1 = QtWidgets.QSpacerItem(40, 150,
                                            QtWidgets.QSizePolicy.Expanding,
                                            QtWidgets.QSizePolicy.Minimum)
        self.control_layout.addItem(spacerItem1)
        # Create FormLayout
        self.form_group_box = QtWidgets.QGroupBox()
        self.formLayout = QtWidgets.QFormLayout()

        self.mt_a_label = QtWidgets.QLabel(self.central_widget)
        self.mt_a_label.setMinimumSize(QtCore.QSize(50, 50))
        self.mt_a_value = QtWidgets.QLabel(self.central_widget)
        self.mt_a_value.setMinimumSize(QtCore.QSize(100, 50))
        self.formLayout.addRow(self.mt_a_label, self.mt_a_value)

        self.mt_d_label = QtWidgets.QLabel(self.central_widget)
        self.mt_d_label.setMinimumSize(QtCore.QSize(50, 50))
        self.mt_d_value = QtWidgets.QLabel(self.central_widget)
        self.mt_d_value.setMinimumSize(QtCore.QSize(100, 50))
        self.formLayout.addRow(self.mt_d_label, self.mt_d_value)

        self.mt_c_label = QtWidgets.QLabel(self.central_widget)
        self.mt_c_label.setMinimumSize(QtCore.QSize(50, 50))
        self.mt_c_value = QtWidgets.QLabel(self.central_widget)
        self.mt_c_value.setMinimumSize(QtCore.QSize(100, 50))
        self.formLayout.addRow(self.mt_c_label, self.mt_c_value)

        self.mt_b_label = QtWidgets.QLabel(self.central_widget)
        self.mt_b_label.setMinimumSize(QtCore.QSize(50, 50))
        self.mt_b_value = QtWidgets.QLabel(self.central_widget)
        self.mt_b_value.setMinimumSize(QtCore.QSize(100, 50))
        self.formLayout.addRow(self.mt_b_label, self.mt_b_value)

        self.form_group_box.setLayout(self.formLayout)
        self.control_layout.addWidget(self.form_group_box)

        self.horizontalLayout.addLayout(self.control_layout)
        self.line.raise_()
        main_window.setCentralWidget(self.central_widget)
        self.menubar = QtWidgets.QMenuBar(main_window)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1188, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        main_window.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(main_window)
        self.statusbar.setObjectName("statusbar")
        main_window.setStatusBar(self.statusbar)
        self.actionOpen_2 = QtWidgets.QAction(main_window)
        self.actionOpen_2.setObjectName("actionOpen_2")
        self.actionSave = QtWidgets.QAction(main_window)
        self.actionSave.setObjectName("actionSave")
        self.menuFile.addAction(self.actionOpen_2)
        self.menuFile.addAction(self.actionSave)
        self.menubar.addAction(self.menuFile.menuAction())

        self.re_translate_ui(main_window)
        QtCore.QMetaObject.connectSlotsByName(main_window)

    def re_translate_ui(self, main_window):
        _translate = QtCore.QCoreApplication.translate
        main_window.setWindowTitle(_translate("main_window", "main_window"))
        self.form_group_box.setTitle(str("MT-Measurement"))
        self.push_button.setText(
            _translate("main_window", "Automatic Selection"))
        self.mt_a_label.setText(_translate("main_window", "MT-A"))
        self.mt_d_label.setText(_translate("main_window", "MT-D"))
        self.mt_c_label.setText(_translate("main_window", "MT-C"))
        self.mt_b_label.setText(_translate("main_window", "MT-B"))
        self.menuFile.setTitle(_translate("main_window", "File"))
        self.actionOpen_2.setText(_translate("main_window", "Open"))
        self.actionSave.setText(_translate("main_window", "Save Image"))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    main_window = QtWidgets.QMainWindow()
    ui = UiMainWindow(main_window)
    main_window.show()
    sys.exit(app.exec_())
