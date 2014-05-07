#! /usr/bin/env python
# *-* coding: utf-8 *-*
"""
The module provides functions to handle the testset provided by LingPy.  The
current testset consists of test files stored in the folder
``lingpy/test/tests/``. These files are formatted according to the requirements
for the different methods. Each specific format comes along with a specific
file extension. Currently, there are the following five extensions:

    1. ``lxs`` -- input files for the :py:class:`LexStat <lingpy.lexstat.LexStat>` class.

    2. ``psq``, ``psq`` -- input files for the :py:class:`Pairwise <lingpy.compare.Pairwise>` class.

    3. ``msq``, ``msa`` -- input files for the :py:class:`Multiple <lingpy.compare.Multiple>` class.

All these files can be easily accessed with help of some specific functions
defined in this module.

"""
from __future__ import division,print_function
import os
from glob import glob

def get_file(filename):
    """
    Return a path to the filename in the testset.

    Parameters
    ----------
    filename : str
        The name of the file (with extension) in the testset.

    Examples
    --------
    >>> from lingpy import *
    >>> get_file('SLAV.lxs')
    '/usr/local/lib/python2.6/dist-packages/lingpy/test/tests/lxs/SLAV.lxs'
    """
    path = os.path.split(os.path.abspath(__file__))[0] + '/tests/'

    name,ext = filename.split('.')

    new_path = path + ext + '/' + name + '.' + ext

    return new_path

def list_files(filetype,dataset='*'):
    """
    List all files in the testset that correspond to a certain filetype.

    Parameters
    ----------
    filetype : { 'lxs', 'msa', 'msq', 'psa', 'psq' }
        The extension of the files that shall be listed. 

    dataset : str (default='*')
        A string which can be used to specify the dataset closer. One can use
        the Unix wildcard syntax in order to narrow down which files to look
        for.
    
    Examples
    --------
    >>> from lingpy import *
    >>> list_files('msa','sindial*')
    sindial_3_1.msa
    sindial_3_2.msa
    sindial_3_3.msa
    sindial_6_1.msa
    """

    path = os.path.split(os.path.abspath(__file__))[0] + '/tests/'
    

    files = sorted(glob(path+filetype+'/'+dataset))

    for f in files: print(os.path.split(f)[1])

