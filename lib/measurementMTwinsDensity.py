# !/usr/bin/python3
# -*- coding: utf-8 -*-
"""Measure the MT density.

Created By Thah.

Edited by Ang at 26th May, 2016, ang.zhou@insa-rennes.fr
Adapted for both Python 2/3 and OSX/Winodws.

Edited by Ang since 12th June, 2016.
Use configparser instead of string reading.
Complete the comment strings.

Edited by Ang since 21th July, 2016.
Auto-trace function was added.
"""

from __future__ import print_function

import configparser
import csv
import logging.config
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as scio
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from matplotlib.colors import LogNorm

from . import read_input_intensity as rii

try:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename
except ImportError:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename

logging.config.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,  # this fixes the problem

    'formatters': {
        'standard': {
            'format': '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        },
    },
    'handlers': {
        'default': {
            'filename': os.path.join(
                os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
            'level': 'INFO',
            'class': 'logging.handlers.RotatingFileHandler',
            'formatter': 'standard',
        },
    },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'DEBUG',
            'propagate': True
        }
    }
})


def _read_config(configdic):
    """Read Config Files."""
    config = configparser.ConfigParser()
    config.read([
        os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
        os.path.join(configdic, 'config.ini')
    ])

    logging.info('reda')
    print(dict(config['db']))
    preference_dict = dict(config['db'])

    sample = config.get('db', 'sample')
    phi_max = int(config.getfloat('db', 'phi'))  # x
    khi_max = int(config.getfloat('db', 'chi'))  # y
    v_min = config.getint('db', 'v_min')
    v_max = config.getint('db', 'v_max')
    neighborhood_size = config.getint('db', 'neighborhood_size')
    threshold = config.getint('db', 'threshold')
    square_size = list(map(int, config.get('db', 'square_size').split(',')))
    square_size = np.asarray(square_size)

    para = (
        phi_max, khi_max, v_min, v_max, sample,
        neighborhood_size, threshold, square_size)

    return preference_dict


def beam_intensity(directory, sample):
    # Get the input beam intensity.
    config = configparser.ConfigParser()
    config.read([
        os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
        os.path.join(directory, 'config.ini')
    ])

    default_beam_intensity_file_list = [
        os.path.join(directory, 'beam1mm.raw'),
        os.path.join(directory, 'beam8mm.raw')]
    if all([os.path.isfile(i) for i in default_beam_intensity_file_list]):
        inp_intensity = np.mean(
            [rii.main(i) for i in default_beam_intensity_file_list])

    else:
        Tk().withdraw()
        beam_intensity_file_list = askopenfilename(
            title='Choose Intensity File...',
            initialdir=directory,
            filetypes=[("Raw files", "*.raw")]
        )
        if beam_intensity_file_list[0]:
            beam_intensity_file = beam_intensity_file_list[0]
            inp_intensity = np.mean(list(map(rii.main, beam_intensity_file)))
        else:
            try:
                inp_intensity = config.getint('db', 'beam_intensity')
            except:
                print("Intensity Error!")
                inp_intensity = 16000 * 8940

    inp_intensity = float(inp_intensity)

    return inp_intensity


def calculation(p, x_offset=0, y_offset=0, sqa=1):
    """Calculation integration.

    Args:
    x: X coordinate parameters.beamIntensity(directory)
    y: Y coordinate parameters.
    square_size: Square side length.
    x_offset: X coordinate error.
    y_offset: Y coordinate error.

    Returns:
    MT_seul: region path chosen to plot.
    """
    [data, x, y, square_size, khi_max, phi_max] = p
    # Initialization Beam Intensity(directory)
    MT_int = []
    MT_points = []
    bruits = []  # contient l'intensity des bruits autour des MT
    aire_bruit = []  # contient nombre de points dans des carrus bruits
    MT_seul = []
    # ----bruit de fond --
    for i in range(len(x)):
        bs = square_size + 2
        # choisir un endroit pour compter le bruit
        square(x[i], y[i], bs, x_offset, y_offset, sqa)
        b, p = _Sq(data, x[i], y[i], bs, x_offset, y_offset)
        bruits.append(b)
        aire_bruit.append(p)

    # --------------------
    for i in range(len(x)):
        square(x[i], y[i], square_size, x_offset, y_offset)
        intensity, points = _Sq(
            data, x[i], y[i], square_size, x_offset, y_offset)  # -bruit
        MT_int.append(intensity)
        MT_points.append(points)

    for i in range(len(x)):
        molecular = (aire_bruit[i] - MT_points[i] + 0.001)
        bf_moyen = (bruits[i] - MT_int[i]) / molecular
        mt_seul = MT_int[i] - bf_moyen * MT_points[i]
        MT_seul.append(mt_seul)

    return bruits, aire_bruit, MT_int, MT_points, MT_seul


def file2array(fichier):
    """Transform file to array."""
    data = scio.loadmat(fichier)
    matrice = np.asarray(data['int'])
    return matrice.transpose()


def square(x, y, size, x_offset, y_offset, sqa=1):
    """Calculation intergration.

    Args:
    x: X coordinate parameters.
    y: Y coordinate parameters.
    size: Square side length.
    x_offset: X coordinate error.
    y_offset: Y coordinate error.

    Returns:
    NaN
    """
    xA, yA = x - size[0] / 2.0 + x_offset, y - size[1] / 2.0 + y_offset
    xB, yB = x + size[0] / 2.0 + x_offset, y - size[1] / 2.0 + y_offset
    xC, yC = x + size[1] / 2.0 + x_offset, y + size[1] / 2.0 + y_offset
    xD, yD = x - size[1] / 2.0 + x_offset, y + size[1] / 2.0 + y_offset
    x = [xA, xB, xC, xD, xA]
    y = [yA, yB, yC, yD, yA]
    x = np.asarray(x)
    y = np.asarray(y)
    if sqa:
        plt.plot(x, y, linewidth=2)


def _Sq(data, x, y, size, x_offset=0, y_offset=0):
    """Calculation integration.

    Args:
    x: X coordinate parameters.
    y: Y coordinate parameters.
    size: Square side length.
    x_offset: X coordinate error.
    y_offset: Y coordinate error.

    Returns:
    intensity: sum of region selected.
    """
    x_min = int(np.floor(x - size[0] / 2.0 + x_offset))
    x_max = int(np.floor(x + size[0] / 2.0 + x_offset))
    y_min = int(np.floor(y - size[1] / 2.0 + y_offset))
    y_max = int(np.floor(y + size[1] / 2.0 + y_offset))

    data = data[y_min:y_max, x_min:x_max]
    intensity = np.sum(data)
    b = len(data)
    p = len(data[0])

    return intensity, b * p


def correction(sample, chi, thickness):
    theta = 14.22  # for GaP
    omega = 14.22  # for GaP
    if thickness:
        thickness = thickness
    else:
        thickness = 900. / 10000.
        print("Thickness was not provided")

    e_Angle = 90. - np.rad2deg(
        np.arccos(
            np.cos(np.deg2rad(chi)) * np.sin(np.deg2rad(theta))))
    i_Angle = 90. - np.rad2deg(
        np.arccos(
            np.cos(np.deg2rad(chi)) * np.sin(np.deg2rad(omega))))
    offset = e_Angle - i_Angle
    eta = 1. / 37.6152  # 1/um at 8.05keV (CXRO)
    p = np.sin(np.deg2rad(e_Angle + offset))
    q = np.sin(np.deg2rad(e_Angle - offset))
    coeff_B = p / (p + q)
    coeff_C = 1. - np.exp(-eta * thickness * (1. / p + 1. / q))
    coeff = coeff_B * (1. / eta) * coeff_C

    return coeff


# def count(data, square, khi_max, phi_max):
#     """
#     Args:
#     data - the image data.
#     square - the region contains peaks
#     Return:
#     intensity - the intensity sum of peak region.
#     points - the pixel number of peak region.
#     """
#     intensity = 0
#     points = 0
#
#     for x in range(khi_max):
#         for y in range(phi_max):
#             if square.contains_point([x, y]):
#                 intensity = intensity + data[x, y]
#                 points = points + 1
#     return intensity, points


def main(directory, int_file, para, showImage=1, ecriture_log=1):
    """Main Function."""
    logger = logging.getLogger(__name__)

    config = configparser.ConfigParser()
    config.read([
        os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
        os.path.join(directory, 'config.ini')
    ])
    (phi_max, khi_max, v_min, v_max, sample,
     neighborhood_size, threshold, square_size) = para
    phi_max = int(phi_max)
    khi_max = int(khi_max)
    log_file = os.path.join(directory, (sample + ".log"))

    # ----- Load data -------

    data = file2array(int_file)
    data = data.T

    # ---- Search local maxima (peaks) ---------

    data_max = filters.maximum_filter(data, neighborhood_size)
    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    x, y = [], []  # x,y contain positions of local maxima

    if not x:
        x = [90, 200, 270, 340]
        y = [10, 16, 19, 15.5]
        xp = [np.linspace(xL - 20, xL + 20, 41) for xL in x]
        yp = [np.linspace(yL - 5, yL + 5, 11) for yL in y]
        index = np.zeros((len(xp[1]), len(yp[1])))
        for i in range(len(xp)):
            for j in range(len(xp[1])):
                xD = [xp[0][j], xp[1][j], xp[2][j], xp[3][j]]
                for k in range(len(yp[1])):
                    yD = [yp[0][k], yp[1][k], yp[2][k], yp[3][k]]
                    intensity, _ = _Sq(data, xD[i], yD[i], [2, 2])
                    index[j, k] = intensity
            xIndex, yIndex = np.unravel_index(index.argmax(), index.shape)
            x[i] = xp[i][xIndex]
            y[i] = yp[i][yIndex]
        p = [data, x, y, square_size, khi_max, phi_max]
        bruits, aire_bruit, MT_int, MT_points, MT_seul = calculation(
            p, x_offset=0, y_offset=0, sqa=0)

        x_Z = [90, 200, 270, 10]
        y_Z = [10, 16, 19, 15.5]
        # If the last point beyond the right border, select from the left
        xp = [np.linspace(xL - 10, xL + 10, 21) for xL in x_Z]
        yp = [np.linspace(yL - 5, yL + 5, 11) for yL in y_Z]
        index_Z = np.zeros((len(xp[1]), len(yp[1])))
        for i in range(len(xp)):
            for j in range(len(xp[1])):
                xD = [xp[0][j], xp[1][j], xp[2][j], xp[3][j]]
                for k in range(len(yp[1])):
                    yD = [yp[0][k], yp[1][k], yp[2][k], yp[3][k]]
                    intensity_Z, _ = _Sq(data, xD[i], yD[i], [2, 2])
                    index_Z[j, k] = intensity_Z
            xIndex_Z, yIndex_Z = np.unravel_index(
                index_Z.argmax(), index_Z.shape)
            x_Z[i] = xp[i][xIndex_Z]
            y_Z[i] = yp[i][yIndex_Z]
        p = [data, x_Z, y_Z, square_size, khi_max, phi_max]
        _, _, _, _, MT_seul_Z = calculation(
            p, x_offset=0, y_offset=0, sqa=0)
        if MT_seul_Z[3] > MT_seul[3]:
            x = x_Z
            y = y_Z

    plt.figure(figsize=(25, 5), dpi=200)
    ax = plt.subplot()
    im = plt.imshow(data, origin="lower", norm=LogNorm(vmin=v_min, vmax=v_max))
    # ------- Measurement of MT and GaP (111) reflections intensity -----
    plt.axis([0, 360, 0, khi_max + 1])

    p = [data, x, y, square_size, khi_max, phi_max]
    bruits, aire_bruit, MT_int, MT_points, MT_seul = calculation(
        p, x_offset=0, y_offset=0)

    plt.plot(x, y, 'ro')
    for i in range(len(x)):
        plt.text(x[i] + square_size[0] / 2 + 2,
                 y[i], "MT = " + str(MT_seul[i]), color="white")
        plt.text(30, 16,
                 "Size= " + str(square_size[0]) + "x" + str(square_size[1]),
                 color="white")

    titre = sample + " Mesure = " + str(MT_seul[0])
    plt.title(titre)
    plt.colorbar(im, ax=ax, extend='max', fraction=0.046, pad=0.04)
    savename = os.path.normpath(
        os.path.join(directory, ("MT_density_" + sample + ".png")))
    plt.savefig(
        savename,
        dpi=500,
        bbox_inches='tight')
    if ecriture_log:
        headers = ["MT-A", "MT-C", "MT-B", "MT-D"]
        MT_seul1 = []
        MT_seul1.append(MT_seul[2])
        MT_seul1.append(MT_seul[3])
        MT_seul1.append(MT_seul[0])
        MT_seul1.append(MT_seul[1])

        row = MT_seul1
        with open(log_file.replace('.log', '_MT.csv'), 'w') as fp:
            f_csv = csv.writer(fp)
            f_csv.writerow(headers)
            f_csv.writerow(row)

        with open(log_file, "w") as f:
            f.write("x: %s" % str(x))
            f.write(os.linesep)
            f.write("y: %s " % str(y))
            f.write(os.linesep)
            f.write("Bruit de fond: %s" % str(bruits))
            f.write(os.linesep)
            f.write(
                "Nombre de points dans chaque carres bruits: %s" % str(
                    aire_bruit))
            f.write(os.linesep)
            f.write("MT intensity: %s" % str(MT_int))
            f.write(os.linesep)
            f.write("MT points: %s" % str(MT_points))
            f.write(os.linesep)
            f.write("MT sans bruits: %s" % str(MT_seul))
            f.write(os.linesep)

    if showImage:
        plt.show()
    inp_intensity = beam_intensity(directory, sample)
    # Fix the error of MT intensity.
    MT_seul1 = [i * 10000 / inp_intensity for i in MT_seul1]
    # Create the table contains the MT intensity.
    mt_table_file = os.path.join(directory, '{0}_result.csv'.format(sample))
    if config.has_option('db', 'thickness'):
        thickness = float(config.get('db', 'thickness'))
    else:
        thickness = None
    with open(mt_table_file, 'w') as tableTeX:
        eta = [correction(sample, x, thickness) for x in y]
        MT_seul1[0] = MT_seul1[0] * 0.939691064
        MT_seul1[1] = MT_seul1[1] * 0.426843274 / (eta[1] / eta[2])
        MT_seul1[2] = MT_seul1[2] * 0.72278158 / (eta[0] / eta[2])
        MT_seul1[3] = MT_seul1[3] * 0.666790711 / (eta[3] / eta[2])
        spamwriter = csv.writer(tableTeX, dialect='excel')
        spamwriter.writerow(MT_seul1)
    return savename


if __name__ == '__main__':
    import logging

    samplePath = os.path.abspath(
        os.path.join(os.path.dirname(sys.argv[0]), 'sample'))
    para = _read_config(samplePath)
    main(samplePath, os.path.join(samplePath, "pf_int"), para)
    print('This is a sample display.')
