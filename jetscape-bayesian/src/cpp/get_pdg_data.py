from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
import os
import errno
import json
from time import gmtime, strftime
import numpy as np
try:
    from urllib.request import urlopen
    from urllib.error import HTTPError
    from urllib.error import URLError
except ImportError:
    from urllib2 import urlopen
    from urllib2 import HTTPError
    from urllib2 import URLError
from shutil import copyfileobj


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


def _format_line(line):
    # Extract particle ids as list
    pids = line[0:32].split()
    pids = [int(i) for i in pids]

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


def has_anti_particle(pid):
    # Handle quarks and leptons
    if 1 <= pid <= 20:
        return True

    # Handle bosons (24 == W+, 34 == W2+, 37 == H+)
    if pid in [24, 34, 37]:
        return True
    elif 20 < pid <= 40:
        return False

    pid = str(pid)

    # If baryon
    if len(pid) >= 4 and pid[-4] != '0':
        return True

    # If quark numbers are different
    quark1 = pid[-2]
    quark2 = pid[-3]
    if quark1 != quark2:
        return True
    return False


# Set up URL and filename format using <current year> as the default
year = int(strftime("%Y", gmtime()))
url = 'http://pdg.lbl.gov/{0}/mcdata/mass_width_{0}.mcd'

# Create data directory if it doesn't exist already
save_dir = '~/.pdg_data/'
save_dir = _create_dir(save_dir)

whitelist = np.loadtxt('{}/whitelist.txt'.format(save_dir), dtype=int)
filter_part = True
if whitelist[0] == 0:
    filter_part = False

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
        response = None
    except URLError:
        response = None
except URLError:
    response = None

# If no response, something went wrong and nothing is left to be done
if not response:
    exit()

# Load file contents and close file
# reponse = response.decode('utf-8')
content = response.readlines()
response.close()

# Strip trailing newlines
content = [x.rstrip() for x in content]

# Set dict to get opposite charge
switch_charge = {'++': '--', '+': '-', '0': '0', '-': '+', '--': '++',
                 '-1/3': '+1/3', '+2/3': '-2/3'}

# Remove header lines and reformat data
# data = {}
data = []
for x in content:
    x = x.decode('utf-8')
    if x[0] != '*':
        line = _format_line(x)
        if not filter_part or line[0][0] in whitelist:
            for i in range(len(line[0])):
                temp_data = {}
                temp_data['pid'] = line[0][i]
                temp_data['mass'] = line[1]
                temp_data['mass_err'] = line[2]
                temp_data['width'] = line[3]
                temp_data['width_err'] = line[4]
                temp_data['name'] = line[5]
                temp_data['charge'] = line[6][i]
                # temp_data = [line[0][i]] + line[1:6] + [line[6][i]]
                # data[temp_data['pid']] = temp_data
                data.append(temp_data)

                if has_anti_particle(line[0][i]):
                    # print('Anti-particle found: {} - {} / charge: {}'.format(line[0][i], line[5], line[6][i]))
                    temp_data = {}
                    temp_data['pid'] = -1 * line[0][i]
                    temp_data['mass'] = line[1]
                    temp_data['mass_err'] = line[2]
                    temp_data['width'] = line[3]
                    temp_data['width_err'] = line[4]
                    temp_data['name'] = line[5]
                    temp_data['charge'] = switch_charge[line[6][i]]
                    # data[temp_data['pid']] = temp_data
                    data.append(temp_data)

# Save JSON file
with open('{}/mass_width.json'.format(save_dir), 'w') as f:
    json.dump(data, f)
