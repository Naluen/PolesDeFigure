"""Plot the Figure_de_Poles.

Script written by Tra NGUYEN thanh-tra.nguyen@insa-rennes.fr
-----04/03/2013 -----
Edited by Ang ZHOU ang.zhou@insa-rennes.fr.
Inproved the adaption to PYTHON 3.
-----26/05/2016 -----
Edited by Ang ZHOU.
Added GUI and configparser reader.
-----22/07/2016 -----
Edited by Ang.
Fixed the bug of PF fig position.
"""
# -*- coding: latin-1 -*-


from __future__ import print_function, unicode_literals

import codecs
import logging
import os
import re
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as scio
from matplotlib.colors import LogNorm

from lib import measurementMTwinsDensity as measurementMTwinsDensity

PY2 = sys.version_info[0] == 2

if PY2:
    # import Tkinter as tk
    # import tkFileDialog as filedialog
    from lib.bruker import convert_raw_to_uxd
    import ConfigParser as configparser
else:
    # import tkinter as tk
    # from tkinter import filedialog
    from lib.bruker3 import convert_raw_to_uxd
    import configparser

logger = logging.getLogger(__name__)


def _generate_data_file(raw_file):
    """Generate Date Files."""
    uxd_file = os.path.join(raw_file.replace("raw", "uxd"))
    try:
        data = codecs.open(uxd_file, 'r', encoding='utf-8', errors='ignore')
    except IOError:
        convert_raw_to_uxd(raw_file, uxd_file)
        data = codecs.open(uxd_file, 'r', encoding='utf-8', errors='ignore')
    int_path = uxd_file.replace(".uxd", "_int.mat")

    step = -1  # flag
    khiData = []
    intData = []
    phiData = []
    flagOfMatrx = -1
    for lines in data.readlines():
        if not lines.strip():
            continue
        elif lines.startswith("_"):
            if lines.startswith("_KHI"):
                line = lines.split("=")
                khi = line[1]
                khi = khi.strip('\n')
                khiData.append(float(khi))
            elif lines.startswith("_STEPTIME"):
                line = lines.split("=")
                steptime = float(line[1])
            elif lines.startswith("_STEPS"):
                line = lines.split("=")
                steps = int(line[1])
                if step == -1:
                    step = steps
                    if step == 0:
                        break
                    else:
                        pass
                else:
                    step = steps
            else:
                continue
        elif lines.startswith(";"):
            if lines.startswith("; ( Data for Range number"):
                flagOfMatrx = flagOfMatrx + 1
                intData.append([])
                phiData.append([])
        else:
            line = list(map(float, lines.split()))
            intData[flagOfMatrx].append(line[1] / steptime)
            phiData[flagOfMatrx].append(line[0])
    scio.savemat(int_path, {'int': intData, 'phi': phiData, 'chi': khiData})

    data.close()
    return int_path


def _shift_angel(intensity, alpha):
    """Shift the intensity fig.

    shift the intensity fig  towards left for alpha degree and add the part
    which exceeded to the right.

    Args:
    intensity: input array.
    alpha: shift degree.

    Returns:
    intensity1: array shifted.
    """
    (xlength, ylength) = intensity.shape
    intensity1 = np.zeros(shape=intensity.shape)
    if alpha is not 0:
        intensity1[:(xlength - alpha) % xlength,
                   :] = intensity[(xlength + alpha) % xlength:, :]
        intensity1[(xlength - alpha) % xlength:,
                   :] = intensity[:(xlength + alpha) % xlength, :]
    else:
        intensity1 = intensity

    return intensity1


def _load_data(int_path, phi_offset):
    """Return the pathname of the KOS root directory."""

    data = scio.loadmat(int_path)
    intensity = np.asarray(data['int']).transpose()
    pp = data['phi'][1, :]
    pp = pp + phi_offset
    cc = data['chi']

    phi = np.asarray(pp)
    phi = phi.T
    phiCopy = phi
    phi = np.radians(phi)
    chi = np.asarray(cc)
    chi = chi.T

    return intensity, phi, chi, phiCopy


def _readCONFIG(directory):
    if os.path.isfile(directory):
        directory = os.path.dirname(directory)
    config = configparser.ConfigParser()
    config.read([
        os.path.join(os.path.dirname(sys.argv[0]), 'config.ini'),
        os.path.join(directory, 'config.ini')
    ])
    if not os.path.isfile(os.path.join(directory, 'config.ini')):
        shutil.copy(
            os.path.join(
                os.path.dirname(sys.argv[0]), 'config.ini'), directory)
    return config


def _readConfigDict(directory):
    """Read the config."""
    print(directory)
    config = _readCONFIG(os.path.dirname(directory))

    if os.path.isfile(directory):
        raw_file = directory
    else:
        raw_file = os.path.join(directory, config.get('db', 'raw_file'))

    directory = os.path.split(raw_file)[0]
    directory = os.path.normcase(directory)
    v_min = config.getfloat('db', 'v_min')
    v_max = config.getfloat('db', 'v_max')
    logscale = config.getint('db', 'logscale')
    phi_offset = config.getfloat('db', 'phi_offset')
    neighborhood_size = config.getint('db', 'neighborhood_size')
    threshold = config.getint('db', 'threshold')
    square_size = list(map(int, config.get('db', 'square_size').split(',')))
    square_size = np.asarray(square_size)

    # Recongnize the sample automatically
    try:
        sample = (re.findall(r'S\d\d\d\d', directory, flags=re.IGNORECASE))
        sample = sample[-1]
        if len(sample) == 0:
            sample = (re.findall(r'\d\d\d\d',
                                 directory, flags=re.IGNORECASE))[-1]
    except:
        # searchOBJ=re.findall(r'S\d\d\d\d', directory, flags=re.IGNORECASE)
        if config.has_option('db', 'sample'):
            sample = config.get('db', 'sample')
        else:
            print("Error, no sample name.")

    return (
        raw_file, directory, v_min, v_max, logscale, phi_offset,
        sample, neighborhood_size, threshold, square_size)


def _setConfig(config, phi, phi_offset, chi, directory, sample):
    """Set Config."""
    config.set("db", "phi", str(phi[-1] - phi_offset))
    config.set("db", "chi", str(chi[-1, -1]))
    config.set("db", "directory", directory)
    config.set("db", "tiki", sample)
    config.write(open(os.path.join(directory, 'config.ini'), "w"))


def plot2D(directory, measurement=True, removeCache=True, showImage=False):
    """Plot 2D Images."""
    # ----- plot 2D -------------------------------------------------------
    (raw_file, directory, v_min, v_max, logscale,
     phi_offset, sample, neighborhood_size,
     threshold, square_size) = _readConfigDict(directory)

    int_path = _generate_data_file(raw_file)

    # uxd_file = raw_file.replace("raw", "uxd")
    values, phi, chi, phiCopy = _load_data(int_path, phi_offset)
    config = _readCONFIG(directory)
    _setConfig(config, phiCopy, phi_offset, chi, directory, sample)

    plt.figure(figsize=(25, 5), dpi=200)
    ax2d = plt.subplot(111)
    (xlength, ylength) = values.shape
    phi_offset = int(phi_offset)
    values2d = values
    values2d = _shift_angel(values2d, 360 - int(phiCopy[-1]))
    im = plt.imshow(
        values2d.T, origin="lower", norm=LogNorm(vmin=v_min, vmax=v_max))
    plt.colorbar(im, ax=ax2d, extend='max', fraction=0.046, pad=0.04)

    ax2d.axis([0, xlength, 0, ylength])
    plt.title(sample + "\n")

    savename = os.path.join(directory, sample + "_2D.png")
    plt.savefig(
        savename,
        dpi=1000,
        bbox_inches='tight')

    # if the image will be shown, default no.
    if showImage:
        plt.show()
    plt.close()

    data = scio.loadmat(int_path)
    data['int'] = values2d.T
    scio.savemat(int_path, data)

    # Measure MT Intensity.
    para = (
        phiCopy[-1] - phi_offset, chi[-1], v_min, v_max, sample,
        neighborhood_size, threshold, square_size)
    if measurement:
        para = measurementMTwinsDensity.main(
            directory, int_path, para, showImage)
    # Remove the cache files,default yes.
    if removeCache:
        os.remove(int_path)

    return para


def plotPF(directory, showImage=0):
    # -- Polar plot... ------------------------------------------------
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

    config = _readCONFIG(directory)
    (raw_file, directory, v_min, v_max, logscale,
     phi_offset, sample, neighborhood_size,
     threshold, square_size) = _readConfigDict(directory)

    ctn_number = config.getint('db', 'ctn_number')
    step = (v_max - v_min) / ctn_number
    contour_levels = np.arange(v_min, v_max, step)

    int_path = _generate_data_file(raw_file)
    values, phi, chi, phiCopy = _load_data(int_path, phi_offset)
    r, theta = np.meshgrid(chi, phi)
    _setConfig(config, phiCopy, phi_offset, chi, directory, sample)

    ax.grid(color="white")
    ax.tick_params(axis="y", labelsize=15, labelcolor="white")
    ax.tick_params(axis="x", labelsize=15)

    if logscale:
        cax = ax.contourf(theta, r, values, contour_levels,
                          norm=LogNorm(vmin=v_min, vmax=v_max))
        cblabel = r'$ Counts\ per\ second $'
        fm = "%.f"
        savename = os.path.join(directory, sample + "_PF.png")
    else:
        cax = ax.contourf(theta, r, values, contour_levels)
        cblabel = r'$Counts\ per\ second $'
        fm = "%.f"
        savename = os.path.join(directory, sample + "_PF.png")

    ax.set_rmin(0.)

    cb = fig.colorbar(cax, pad=0.1, format=fm, extend='max')
    # make a color bar

    # ---- colorbar tick font size
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(15)

    cb.set_label(cblabel, size=15)  # put a label on CB
    plt.title(sample + "\n")
    plt.savefig(
        savename,
        dpi=1000,
        bbox_inches='tight')
    # if the image will be shown, default no.
    if showImage:
        plt.show()
    plt.close('all')

    os.remove(int_path)
    return savename


def main(directory):
    """Read data from config.ini.

    It will first scan the direcory of raw file.
    If none, it will use the default ini in program directory.

    Args:
        directory: The directory of raw file.
    """

    plot2D(directory)
    plotPF(directory)
    print('Finished')

    return directory


if __name__ == '__main__':
    main(os.path.abspath(
        os.path.join(os.path.dirname(sys.argv[0]), 'sample')), 1, 0, 0)
    print('This is a sample display.')
