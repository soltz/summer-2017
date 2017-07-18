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
