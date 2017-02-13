from unittest import TestCase
from ..Poles import raw_file_reader

class TestPolesFIGURE(TestCase):
    def test_raw_file_reader(self):
        data = raw_file_reader('..sample/PF.raw')
        print(data)
