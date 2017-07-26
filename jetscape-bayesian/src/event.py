# Author: Matthias Heinz
"""
High-energy nucleus-nucleus event data simulation module.

Abstracted interface to PYTHIA for jet simulation and TRENTO for background
simulation. Full event records are saved to timestamped files, and methods
return the absolute path to the file for later use.

class EventGenerator
--------------------

Abstracted interface to PYTHIA to quickly simulate jet event records for many
events with the same settings. The default settings are chosen to produce jets
at RHIC energies. It has the following methods::

    jet_generator = event.EventGenerator(save_directory, preset,
                                         pythia_settings)
    filename = jet_generator.generate()

class BackgroundGenerator
-------------------------

Abstracted interface to TRENTO to quickly simulate background event records for
many events. The default settings are chosen to match multiplicities at RHIC
energies. It has the following methods::

    background_generator = event.BackgroundGenerator(save_directory,
                                                     projectile1, projectile2)
    filename = background_generator.generate()

Changelog:

2017.07.18
    Initial creation
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
import os
import errno
import numpy as np
from time import gmtime, strftime
import pythia8
import subprocess
from future.moves.urllib.request import urlopen
from future.moves.urllib.error import HTTPError
from future.moves.urllib.error import URLError
from shutil import copyfileobj


class EventGenerator(object):

    def __init__(self, save_directory=None, preset='RHIC',
                 pythia_settings=None):
        pass

    def generate(self):
        return ''


class BackgroundGenerator(object):

    def __init__(self, save_directory=None, projectile1='Au',
                 projectile2='Au'):
        pass

    def generate(self):
        return ''


# ----------------------------- Internal Methods ---------------------------- #


def _load_pdg_data(make_dict=True):
    # Set up URL and filename format using <current year> as the default
    year = int(strftime("%Y", gmtime()))
    url = 'http://pdg.lbl.gov/{0}/mcdata/mass_width_{0}.mcd'
    filename = 'mass_width_{}.mcd'

    # Create data directory if it doesn't exist already
    save_dir = '~/.pdg_data/'
    save_dir = _create_dir(save_dir)

    # Try getting data from URL with <current year>
    try:
        response = urlopen(url.format(year))
    except HTTPError:
        # If <current year> does not exist, look for previous year
        # (previous year is guaranteed to exist)
        try:
            reponse = urlopen(url.format(year - 1))
            year -= 1
        except HTTPError:
            print('Warning:')
            print('\tUnable to get most recent particle data.')
            print('\tPDG site is unexpectedly causing an HTTPError.')
            print('\tAttempting to load a local copy of the data.\n')
            response = None
        except URLError:
            print('Warning:')
            print('\tUnable to get most recent particle data.')
            print('\tYou are probably not connected to the internet.')
            print('\tAttempting to load a local copy of the data.\n')
            response = None
    except URLError:
        print('Warning:')
        print('\tUnable to get most recent particle data.')
        print('\tYou are probably not connected to the internet.')
        print('\tAttempting to load a local copy of the data.\n')
        response = None

    # If response exists, save it to file
    if response:
        with open(save_dir + '/' + filename.format(year), 'wb') as f:
            copyfileobj(response, f)
        response.close()

    # Open most recent local copy of mass data
    file = None
    while year >= 2017 and not file:
        try:
            file = open(save_dir + '/' + filename.format(year), 'r')
        # Note: terrible catch-all required for Python 2 and 3 compatiblity
        #       should look into alternatives that don't surpress
        #       unexpected behavior
        except Exception:
            year -= 1

    # If file failed to open, return False
    if not file:
        print('Warning:')
        print('\tUnable to load data from local file.\n')
        return False

    # Load file contents and close file
    content = file.readlines()
    file.close()

    # Strip trailing newlines
    content = [x.rstrip() for x in content]

    # Remove header lines and reformat data
    data = []
    for x in content:
        if x[0] != '*':
            line = _format_line(x)
            for i in range(len(line[0])):
                temp_data = [line[0][i]] + line[1:6] + [line[6][i]]
                data.append(temp_data)

    if make_dict:
        # Make data a dict for faster retrieval
        data_dict = {}
        for x in data:
            data_dict[x[5] + x[6]] = x
        data = data_dict

    return data


def _format_line(line):
    # Extract particle ids as list
    pids = line[0:32].split()

    # Extract mass
    mass = float(line[33:51])

    # Extract errors on mass as list
    m_err = [float(x) for x in line[52:69].split()]

    # Extract width if it exists
    width = line[70:88]
    if width.isspace():
        width = None
    else:
        width = float(width)

    # Extract errors on width as list if they exist
    w_err = line[89:106]
    if w_err.isspace():
        w_err = None
    else:
        w_err = [float(x) for x in w_err.split()]

    # Extract name
    name = line[107:].split()[0]

    # Allow the rest to be list of charges
    charges = line[107:].split()[1].strip().split(',')

    return [pids, mass, m_err, width, w_err, name, charges]


def _generate_cmf():
    pass


def _create_dir(path):
    if path[0] == '~':
        path = path[2:]
        home = os.path.expanduser('~')
        path = os.path.join(home, path)
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return str(os.path.abspath(path))
