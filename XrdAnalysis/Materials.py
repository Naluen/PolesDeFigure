import abc


class Material(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass


class GaP(Material):
    def __init__(self):
        super(GaP, self).__init__()

    @property
    def hkl(self):
        hkl = {
            "002": 32,
            "004": 68,
            "006": 115,
            "-2-24": 87
        }
        return hkl

    @property
    def form_name(self):
        return 'GaP'

    @property
    def color(self):
        return '#99E1D9'


class Si(Material):
    @property
    def form_name(self):
        return "Si(6\u00b0 off)"

    @property
    def color(self):
        return "#2F323A"


class AlGaP(Material):
    @property
    def form_name(self):
        return "$\mathrm{Al_{0.2}GaP}$"

    @property
    def color(self):
        return "#DBDBDB"


class GaPN(Material):
    @property
    def form_name(self):
        return "GaPN"

    @property
    def color(self):
        return "#0F7173"
