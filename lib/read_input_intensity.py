#!/usr/bin/env python3
"""Read the beamfile to get the InputIntensity."""
import os
import sys


from lib.bruker3 import convert_raw_to_uxd


def main(beamfile):
    """Read the beamfile to get the InputIntensity.

    Arg:
    beamfile: The input beamfile like beam1mm or beam8mm.

    Return:
    maxInt : the beam intensity maxium.
    """
    try:
        data = open(beamfile.replace('.raw', '.uxd'), 'r')
    except IOError:
        convert_raw_to_uxd(beamfile, beamfile.replace('.raw', '.uxd'))
        data = open(beamfile.replace('.raw', '.uxd'), 'r')

    intData = []
    phiData = []
    for lines in data.readlines():
        if not lines.strip():
            continue
        elif lines.startswith("_"):
                continue
        elif lines.startswith(";"):
                continue
        else:
            line = list(map(float, lines.split()))
            intData.append(line[1])
            phiData.append(line[0])
        maxInt = max(intData)*8940
    return maxInt


if __name__ == '__main__':

    beamfile = os.path.join(os.path.abspath('sample'), 'beam1mm.raw')
    main(beamfile)
    print('This is a sample display.')
