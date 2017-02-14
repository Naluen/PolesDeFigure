import os
import sys
import codecs
import numpy as np
from lib.bruker3 import DatasetDiffractPlusV3
import re

ds = DatasetDiffractPlusV3(open(r'C:\Users\ang\Documents\GitHub\PolesDeFigure\lib\sample\PF.raw', 'rb'))
data = (ds.pretty_format(print_header=True))
data = data.split('\n')
step_time_matrix ={float(i.split("=")[1].strip('\n')) for i in data if i.startswith("_STEPTIME")}
scan_data_matrix = np.asarray(
            [list(map(float, i.split())) for i in data if re.match('\d', i)]
        )
print(scan_data_matrix)
print(step_time_matrix)
int_data_matrix = scan_data_matrix[:, 1]/next(iter(step_time_matrix))
print(int_data_matrix)