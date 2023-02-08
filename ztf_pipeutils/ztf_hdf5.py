import h5py
import ast
from astropy.table import Table, Column, vstack
import astropy.units as u
import astropy
import numpy as np
import simsurvey
import ztf_pipeutils.simsurvey_tools as sst
import matplotlib.pylab as plt
import pandas as pd
import os


class Write_LightCurve_cumul:
    "Definition of a class that records light curves and/or their meta data in an hdf5 file"

    def __init__(self, outputDir='dataLC', file_data='Data.hdf5', file_meta='Meta.hdf5', path_prefix='SN', List=['z', 't0',
                                                                                                                 'x0', 'x1', 'c', 'mwebv', 'ra', 'dec', 'mwebv_sfd98', 'idx_orig'], **kwargs):
        """
        Parameters
        --------------
        outputDir: str, opt
          output directory (default: dataLC)
        file_data : str, opt
            Name of the data file that you want to write (default='Data.hdf5').
        file_meta : str, opt
            Name of the meta data file that you want to write (default='Meta.hdf5').
        path_prefix: str, opt
          prefix for the path in hdf5 files (default: SN)
        """

        self.List = List

        self.outputDir = outputDir
        self.path_prefix = path_prefix

        self.data_name = file_data
        self.meta_file = os.path.join(self.outputDir, file_meta)
        self.data_file = os.path.join(self.outputDir, file_data)

        if not os.path.exists(self.outputDir):
            os.makedirs(self.outputDir)
        if os.path.isfile(self.data_file):
            os.remove(self.data_file)
        if os.path.isfile(self.meta_file):
            os.remove(self.meta_file)

        self.file_data = h5py.File('{}'.format(self.data_file), 'w')
        #self.file_meta = h5py.File('{}'.format(self.meta_file), 'w')

        self.metaData = Table()

    def write(self, lc, metaSN, metaRej, path='SN', pathRej='Rej', serialize_meta=True):
        """
        method to save light curves and metadata

        Parameters
        --------------
        lc: 
          SN light curves
        metaSN: 
          metadata corresponding to LC
        metaRej: 
          metadata corresponding to rejected simu

        """
        for j, llc in enumerate(lc):
            if llc is not None:
                if llc.meta is not None:
                    # update metadata with path
                    mmeta = metaSN[j]
                    thepath = '{}_{}'.format(path, j)
                    metadict = {name: mmeta[name]
                                for name in mmeta.dtype.names}
                    metadict['path'] = thepath
                    # update total metaData
                    t = Table([metadict])
                    self.metaData = vstack(
                        [self.metaData, t])
                    # write lc
                    llc.meta = metadict
                    astropy.io.misc.hdf5.write_table_hdf5(
                        llc, self.file_data, path=thepath, overwrite=True, serialize_meta=serialize_meta)

        if metaRej is not None:
            meta_rej = Table(metaRej)
            pp = ['{}_{}'.format(pathRej, i)
                  for i in range(len(meta_rej))]
            col = Column(pp, name='path')  # shape=(2,)
            meta_rej.add_column(col)
            self.metaData = vstack([self.metaData, meta_rej])

        self.metaData.meta["directory"] = self.outputDir
        self.metaData.meta["file_name"] = self.data_name
        #print('meta data', self.metaData.meta, self.metaData)

        if os.path.isfile(self.meta_file):
            os.remove(self.meta_file)
        file_meta = h5py.File('{}'.format(self.meta_file), 'w')

        astropy.io.misc.hdf5.write_table_hdf5(
            self.metaData, file_meta, path='meta', overwrite=True, serialize_meta=serialize_meta)


class Write_LightCurve:
    "Definition of a class that records light curves and/or their meta data in an hdf5 file"

    def __init__(self, outputDir='dataLC', file_data='Data.hdf5', file_meta='Meta.hdf5', path_prefix='SN', List=['z', 't0',
                                                                                                                 'x0', 'x1', 'c', 'mwebv', 'ra', 'dec', 'mwebv_sfd98', 'idx_orig'], **kwargs):
        """
        Parameters
        --------------
        outputDir: str, opt
          output directory (default: dataLC)
        file_data : str, opt
            Name of the data file that you want to write (default='Data.hdf5').
        file_meta : str, opt
            Name of the meta data file that you want to write (default='Meta.hdf5').
        path_prefix: str, opt
          prefix for the path in hdf5 files (default: SN)
        """

        self.Summary = Table()
        self.List = List

        self.outputDir = outputDir
        self.path_prefix = path_prefix

        self.data_name = file_data
        self.meta_file = os.path.join(self.outputDir, file_meta)
        self.data_file = os.path.join(self.outputDir, file_data)

        if not os.path.exists(self.outputDir):
            os.makedirs(self.outputDir)
        if os.path.isfile(self.data_file):
            os.remove(self.data_file)
        if os.path.isfile(self.meta_file):
            os.remove(self.meta_file)

        self.file_data = h5py.File('{}'.format(self.data_file), 'w')
        self.file_meta = h5py.File('{}'.format(self.meta_file), 'w')

    def write_data(self, lc, meta, meta_rejected, path='', path_rej='', serialize_meta=True):
        """
        Parameters
        ----------
         lc : LightcurveCollection
            List of AstropyTable with simulated lightcurve.
        metadata: numpy array
           metadata for lc
        path : str-
            Name of the key of the data file that you want to write.

        """
        if isinstance(lc, astropy.table.table.Table):

            # update metadata with path
            metadict = {name: meta[name] for name in meta.dtype.names}
            metadict['path'] = path
            self.Summary = vstack([self.Summary, Table([metadict])])
            # write lc
            lc.meta = metadict
            astropy.io.misc.hdf5.write_table_hdf5(
                lc, self.file_data, path=path, overwrite=True, serialize_meta=serialize_meta)
        else:
            #print('lc is not an astropy.table.table.Table type', type(lc), lc)
            inum = -1
            meta_rej = None
            path_rej = []
            for j, llc in enumerate(lc):
                if llc is not None:
                    if llc.meta is not None:
                        path_j = '{}_{}'.format(path, j)
                        self.write_lc(llc, path_j, inum, meta, j)
                    meta_rej = self.concat_metarej(meta_rejected, j, meta_rej)
                    path_rej.append('{}'.format(j))

            self.write_meta(meta_rej, path_rej)

    def write_lc(self, llc, path, inum, meta, j):

        for i, lc_b in enumerate(llc):
            if len(lc_b) > 0:
                inum += 1
                simpath = '{}{}_{}'.format(
                    self.path_prefix, path, inum)
                self.write_data(
                    lc_b, meta[j][i], None, simpath)

    def concat_metarej(self, meta_rejected, j, meta_rej):

        pp = None
        if meta_rejected is not None:
            pp = meta_rejected[j]
        if pp is not None:
            if meta_rej is None:
                meta_rej = pp
            else:
                meta_rej = np.concatenate((meta_rej, pp))

        return meta_rej

    def write_meta(self, meta_rej, path_rej=''):
        """
        write meta data in hdf5 file

        Parameters
        ---------------
        meta_rej: numpy array
          metadata of rejected LC (i.e. not simulated)

        """

        if meta_rej is not None and len(meta_rej) > 0:

            meta_rej = Table(meta_rej)
            pp = ['Rej{}_{}'.format(path_rej[i], i)
                  for i in range(len(meta_rej))]
            col = Column(pp, name='path')  # shape=(2,)
            meta_rej.add_column(col)
            self.Summary = vstack([self.Summary, meta_rej])

        self.Summary.meta["directory"] = self.outputDir
        self.Summary.meta["file_name"] = self.data_name
        print('meta data', self.Summary.meta, self.Summary)

        astropy.io.misc.hdf5.write_table_hdf5(self.Summary, self.file_meta, path='meta', append=True,
                                              overwrite=True, serialize_meta=True)


class Read_LightCurve:

    def __init__(self, file_name='Data.hdf5', inputDir='dataLC'):
        """
        Parameters
        ----------
        file_name : str
            Name of the hdf5 file that you want to read.
        """
        self.file_name = file_name
        self.file = h5py.File('{}/{}'.format(inputDir, file_name), 'r')

    def get_path(self):
        """
        Method to return the list of keys of the hdf5 file

        Returns
        ----------
        list(str): list of keys (aka paths)

        """

        return list(self.file.keys())

    def get_table(self, path):
        """
        Parameters
        ----------
        path : str
            hdf5 path for light curve.

        Returns
        -------
        AstropyTable
            Returns the reading of an .hdf5 file as an AstropyTable.
        """

        tab = Table()
        try:
            tab = astropy.io.misc.hdf5.read_table_hdf5(
                self.file, path=path, character_as_bytes=False)
        except (OSError, KeyError):
            pass

        return tab


# In[257]:


class Plot_LightCurve:

    def __init__(self, lc_sn=None):
        """
        Parameters
        ----------
        lc_sn : AstropyTable, ops
            AstropyTable containing light curves and their metadata (default=None).
        """
        self.lc_sn = lc_sn

    def Plot(self, file_name, path, figsize):  # axs
        """
        Parameters
        ----------
        file_name : str
            name of the file to process.
        path : str
            hdf5 path for light curve.
        figsize : tuple (float, float)
            Size of my figure.
        """
        Self = Read_LightCurve()
        Data = Self.Read_hdf5_file(file_name=file_name, path=path)

        band = ['ztfr', 'ztfg', 'ztfi']
        dico = dict(zip(band, [(0, 0), (0, 1), (1, 0)]))

        fig, axs = plt.subplots(2, 2, figsize=figsize)
        for color in band:

            mask_band = Data['band'] == color
            new_lc_sn = Data[mask_band]
            time = new_lc_sn['time']
            flux = new_lc_sn['flux']
            err = new_lc_sn['fluxerr']

            ipos = dico[color][0]
            jpos = dico[color][1]

            axs[ipos, jpos].errorbar(time, flux, yerr=err, fmt='.',
                                     ecolor='red', capsize=5, alpha=0.5)
            axs[ipos, jpos].plot(time, flux)
            axs[ipos, jpos].set_title('Light Curve for {} band'.format(color))
            axs[ipos, jpos].grid(True)
            plt.subplots_adjust(hspace=0.3)

    def Histo(self, file_name, path, figsize):
        """
        Parameters
        ----------
        file_name : str
            name of the file to process
        path : str
            hdf5 path for light curve
        figsize : tuple (float, float)
            Size of my figure.
        """

        Self = Read_LightCurve()
        Data = Self.Read_hdf5_file(file_name=file_name, path=path)
        var = ['z', 'x0', 'x1', 'c']
        dico = dict(zip(var, [(0, 0), (0, 1), (1, 0), (1, 1)]))

        fig, axs = plt.subplots(2, 2, figsize=figsize)
        for variables in var:
            ipos = dico[variables][0]
            jpos = dico[variables][1]

            axs[ipos, jpos].hist(Data[variables], alpha=0.3)
            axs[ipos, jpos].set_title('Histo var {}'.format(variables))
