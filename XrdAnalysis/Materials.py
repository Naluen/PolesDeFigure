import abc

class Material(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass


class GaP(Material):
    def __init__(self):
        super(GaP, self).__init__()
        self.hkl = {
            "002": 32,
            "004": 68,
            "006": 115,
            "-2-24": 87
        }