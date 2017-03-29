import logging

import h5py

from RawReader import RawReader as RawReader


class H5Reader(RawReader):
    def __init__(self, h5_file):
        super(H5Reader, self).__init__(h5_file)
        if not h5_file.endswith('.h5'):
            raise TypeError("Illegal h5 file name.")

    def create_file(self):
        with open(self.file, 'a'):
            file_handle = h5py.File(self.file)
        file_handle.attrs['HDF5_Version'] = h5py.version.hdf5_version
        file_handle.attrs['h5py_version'] = h5py.version.version

    def read_raw(self, raw_file):
        # Require data.
        instance = RawReader(raw_file)
        data_dict, plot_dict = instance.matrix_data()
        scan_dict = instance.get_head()
        # Require data set.
        file_handle = h5py.File(self.file, 'a')
        sample_name_str = plot_dict['db']['sample']
        sample_handle = file_handle.require_group(sample_name_str)
        lr_handle = sample_handle.require_group(scan_dict['_TYPE'])
        # Record data.
        for key, value in data_dict.items():
            if key in lr_handle:
                del lr_handle[key]

            lr_handle.create_dataset(
                key,
                data=value
            )
        for key, value in scan_dict.items():
            try:
                lr_handle.attrs.modify(key, value)
            except TypeError:
                pass
        if plot_dict:
            for key, value in plot_dict['db'].items():
                try:
                    lr_handle.attrs.modify(key, value)
                except TypeError:
                    pass
        # Close data set.
        file_handle.close()

    def matrix_data(self):
        pass


if __name__ == '__main__':
    # Setup Log
    logging.basicConfig(
        # filename=os.path.join(
        #     os.path.dirname(sys.argv[0]), 'log', __name__ + '.log'),
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
    )

    lib_file_name = r'mylib.h5'
    logging.info("File {0} was chosen.".format(lib_file_name))
    sample = H5Reader(lib_file_name)
    sample.read_rsm(
        r'C:\Users\ang\Dropbox\Experimental_Data\SSMBE\S1866\HR\002.raw')
    logging.info(
        "Finished!\n"
        "--------------------------------------------------------------------"
    )
