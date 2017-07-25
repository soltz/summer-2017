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

def _load_pdg_data():
    pass

def _generate_cmf():
    pass
