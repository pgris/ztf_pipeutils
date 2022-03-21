import multiprocessing
import numpy as np
import pandas as pd
from astropy.table import Table, vstack
import operator
import yaml
from collections import MutableMapping
import os


def multiproc(data, params, func, nproc, gather=True):
    """
    Function to perform multiprocessing
    Parameters
    ---------------
    data: array
      data to process
    params: dict
      fixed parameters of func
    func: function
      function to apply for multiprocessing
    nproc: int
      number of processes
    """
    nproc = min([len(data), nproc])
    # multiprocessing parameters
    nz = len(data)
    t = np.linspace(0, nz, nproc+1, dtype='int')
    # print('multi', nz, t)
    result_queue = multiprocessing.Queue()

    procs = [multiprocessing.Process(name='Subprocess-'+str(j), target=func,
                                     args=(data[t[j]:t[j+1]], params, j, result_queue))
             for j in range(nproc)]

    for p in procs:
        p.start()

    resultdict = {}
    # get the results in a dict

    for i in range(nproc):
        resultdict.update(result_queue.get())

    for p in multiprocessing.active_children():
        p.join()

    if gather:
        restot = gather_results(resultdict)
    else:
        restot = resultdict
    return restot


def gather_results(resultdict):
    """
    Function to gather results of a directory
    Parameters
    ----------------
    resultdict: dict
      dictory of data
    Returns
    ----------
    gathered results. The type is determined from resultdict.
    Supported types: pd.core.frame.DataFrame, Table, np.ndarray,
    np.recarray, int
    """
    supported_types = ['pd.core.frame.DataFrame', 'Table', 'np.ndarray',
                       'np.recarray', 'int', 'dict']

    # get outputtype here
    first_value = None
    for key, vals in resultdict.items():
        if vals is not None:
            first_value = vals
            break

    restot = None
    if first_value is None:
        return restot

    if isinstance(first_value, pd.core.frame.DataFrame):
        restot = pd.DataFrame()

        def concat(a, b):
            return pd.concat((a, b), sort=False)

    if isinstance(first_value, Table):
        restot = Table()

        def concat(a, b):
            return vstack([a, b])

    if isinstance(first_value, np.ndarray) or isinstance(first_value, np.recarray):
        restot = []

        def concat(a, b):
            if isinstance(a, list):
                return b
            else:
                return np.concatenate((a, b))

    if isinstance(first_value, int):
        restot = 0

        def concat(a, b):
            return operator.add(a, b)

    if isinstance(first_value, dict):
        restot = {}

        def concat(a, b):
            return dict(a, **b)

    if restot is None:
        print('Sorry to bother you but: unknown data type', type(first_value))
        print('Supported types', supported_types)
        return restot

    # gather the results
    for key, vals in resultdict.items():
        restot = concat(restot, vals)

    return restot


def make_dict_from_config(path, config_file):
    """
    Function to make a dict from a configuration file
    Parameters
    ---------------
    path: str
      path dir to the config file
    config_file: str
       config file name
    Returns
    ----------
    dict with the config file infos
    """

    # open and load the file here
    ffile = open('{}/{}'.format(path, config_file), 'r')
    line = ffile.read().splitlines()
    ffile.close()

    # process the result
    params = {}
    for i, ll in enumerate(line):
        if ll != '' and ll[0] != '#':
            spla = ll.split('#')
            lspl = spla[0].split(' ')
            lspl = ' '.join(lspl).split()
            n = len(lspl)
            keym = ''
            lim = n-2
            for io, keya in enumerate([lspl[i] for i in range(lim)]):
                keym += keya
                if io != lim-1:
                    keym += '_'
            params[keym] = (lspl[n-1], lspl[n-2], spla[1])
    return params


def make_dict_from_optparse(thedict):
    """
    Function to make a nested dict from a dict
    The idea is to split the original dict key(delimiter: _)  to as many keys
    Parameters
    --------------
    thedict: dict
    Returns
    ----------
    final dict
    """

    params = {}
    for key, vals in thedict.items():
        lspl = key.split('_')
        n = len(lspl)
        mystr = ''
        myclose = ''
        for keya in [lspl[i] for i in range(n)]:
            mystr += '{\''+keya + '\':'
            myclose += ' }'

        if vals[0] != 'str':
            dd = '{} {} {}'.format(mystr, eval(
                '{}({})'.format(vals[0], vals[1])), myclose)
        else:
            dd = '{} \'{}\' {}'.format(mystr, vals[1], myclose)

        thedict = eval(dd)
        params = recursive_merge(params, thedict)

    return params


def recursive_merge(d1, d2):
    """
    Update two dicts of dicts recursively, 
    if either mapping has leaves that are non-dicts, 
    the second s leaf overwrites the first s.
    Parameters
    ---------------
    d1, d2: dicts to merge
    Returns
    -----------
    merged dict
    """
    for k, v in d1.items():  # in Python 2, use .iteritems()!
        if k in d2:
            # this next check is the only difference!
            if all(isinstance(e, MutableMapping) for e in (v, d2[k])):
                d2[k] = recursive_merge(v, d2[k])
                # we could further check types and merge as appropriate here.
    d3 = d1.copy()
    d3.update(d2)
    return d3


def dump_in_yaml(opts, confDict, dirOut, nameOut, prefix='simu'):

    # load the new values
    newDict = {}
    for key, vals in confDict.items():
        newval = eval('opts.{}'.format(key))
        newDict[key] = (vals[0], newval)

    # new dict with configuration params
    yaml_params = make_dict_from_optparse(newDict)
    nameOut = '{}_{}'.format(prefix, nameOut)
    yaml_name = '{}/{}'.format(dirOut, nameOut)
    with open(yaml_name, 'w') as f:
        data = yaml.dump(yaml_params, f)


def checkDir(outDir):
    """
    function to check whether a directory exist
    and create it if necessary
    """
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
