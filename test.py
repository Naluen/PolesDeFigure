import os
import sys
import codecs

with open(
    os.path.join(
        os.path.dirname(sys.argv[0]), 'lib', 'sample', 'PF.raw'),
    'rb'
) as fp:
    data = fp.read(480)
print(codecs.decode(data, encoding='utf-8', errors='ignore'))