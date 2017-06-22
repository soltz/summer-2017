from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
import os
import errno
import math
import random
import pythia8
import numpy as np
import argparse
import subprocess
from time import gmtime, strftime


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


def _save_event(event, filename):
    with open(filename, 'w') as f:
        f.write('--------  PYTHIA Event Listing  (complete event)  ---------------------------------------------------------------------------------\n\n')
        f.write('    no        id   name            status     mothers   daughters     colours      p_x        p_y        p_z         e          m  \n')
        for i in range(event.size()):
            no = i
            id = event[i].id()
            name = event[i].name()
            status = event[i].status()
            mothers = (event[i].mother1(), event[i].mother2())
            daughters = (event[i].daughter1(), event[i].daughter2())
            colors = (event[i].col(), event[i].acol())
            px = event[i].px()
            py = event[i].py()
            pz = event[i].pz()
            e = event[i].e()
            m = event[i].m()
            f.write('{:>6}{:>10}   {:<12}{:>10}{:>6}{:>6}{:>6}{:>6}{:>6}{:>6}{:>11.3f}{:>11.3f}{:>11.3f}{:>10.3f}{:>11.3f}\n'.format(no, id, name, status, mothers[0], mothers[1], daughters[0], daughters[1], colors[0], colors[1], px, py, pz, e, m))
        f.write('--------  End PYTHIA Event Listing  -----------------------------------------------------------------------------------------------\n')
    print('Filtered event saved to {}'.format(filename))


# -------------------------- ARGUMENT PARSING ------------------------ #

prog_description = '''
Generates output below for PYTHIA jets in TRENTO background with
slowjet reconstruction.'''.replace('\n', '', 1)
prog_epilogue = 'Writes output to output.txt, see header for format.'

# Initialize argument parser
parser = argparse.ArgumentParser(
    description=prog_description,
    usage='%(prog)s [options]',
    epilog=prog_epilogue,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add arguments
parser.add_argument('-f', '--file', type=str, metavar='TRENTO_FILE',
                    default='generate',
                    help='path to TRENTO data file or turn \"off\"')
parser.add_argument('-e', '--eCM', type=float, default=200.0,
                    metavar='BEAM_ENERGY',
                    help='PYTHIA beam center-of-mass energy (in GeV)')
parser.add_argument('-n', '--nevt', type=int, default=10,
                    metavar='NUMBER_OF_EVENTS',
                    help='number of PYTHIA and TRENTO events to generate')
parser.add_argument('-c', '--QCDoff', action='store_true',
                    help='turn PYTHIA hard QCD processes off')
parser.add_argument('-q', '--QEDoff', action='store_true',
                    help='turn PYTHIA hard QED processes off')
parser.add_argument('-o', '--output', type=str, metavar='OUTPUT_FILE',
                    default='auto',
                    help='path to output file')
parser.add_argument('-s', '--seed', type=int, default=-1,
                    help='seed for random numbers in PYTHIA')
parser.add_argument('-p', '--SJpTmin', type=int, default=10,
                    help='slowjet minimum pT')
parser.add_argument('-r', '--SJradius', type=float, default=0.5,
                    help='slowjet radius')
parser.add_argument('-u', '--quench', type=float, metavar='QUENCH_FACTOR',
                    default=1.0,
                    help='scaling factor for momentum of non-photon jet')
parser.add_argument('-y', '--pTHatMin', type=float, default=20.0,
                    help='PYTHIA minimum jet pT')
parser.add_argument('-z', '--pTHatMax', type=float, default=25.0,
                    help='PYTHIA maximum jet pT')
parser.add_argument('-w', '--filter', action='store_true',
                    help='Enable filtering results')
parser.add_argument('-x', '--filter-file', type=str,
                    metavar='FILTER_FILE', default='auto',
                    help='Set path to file for filtered results')

# Parse arguments
args = parser.parse_args()

# Check validity of arguments
if args.file != 'off' and args.file != 'generate':
    if not os.path.isfile(args.file):
        raise ValueError('TRENTO file {} not found.'.format(args.file))

if args.eCM <= 0.0:
    raise ValueError('Beam energy must be positive')

if args.nevt <= 0:
    raise ValueError('Number of events must be positive')

if args.SJpTmin < 0:
    raise ValueError('Minimum slowjet pT must be non-negative')

if args.SJradius <= 0:
    raise ValueError('slowjet radius must be positive')

if args.quench <= 0 or args.quench > 1.0:
    raise ValueError('Quench factor must be between 0 and 1')

if args.pTHatMin < 0:
    raise ValueError('Minimum pT for PYTHIA must be non-negative')

if args.pTHatMax <= args.pTHatMin:
    raise ValueError('Maximum pT must be larger than minimum')

if args.output == args.file:
    raise ValueError('Cannot read from and output to same file')

if args.filter and args.filter_file == args.output and args.output != 'auto':
    raise ValueError('Filter file cannot be the same as output file')

if args.filter and args.filter_file == args.file:
    raise ValueError('Cannot read from and filter to the same file')

if not args.filter and args.filter_file != 'auto':
    print('Warning: Filter flag is not set so filter file will not be created')

if args.output == 'auto' or args.filter:
    # Set up data directory
    cwd = os.getcwd()
    path = os.path.join(cwd, 'data')
    date_time = strftime("%Y.%m.%d_%H.%M.%S", gmtime())
    path = os.path.join(path, date_time)
    _create_dir(path)
    filtered_events_path = os.path.join(path, 'filtered_events')
    _create_dir(filtered_events_path)

# Generate appropriate output and filter file paths
if args.output == 'auto':
    output_path = os.path.join(path, '{}_{}.txt'.format(args.quench,
                                                        args.nevt))
else:
    output_path = args.output
if args.filter_file == 'auto':
    filter_path = os.path.join(path, '{}_{}_filtered.txt'.format(args.quench,
                                                                 args.nevt))
else:
    filter_path = args.filter_file

# Check last to prevent changing files on the system when other
# arguments might still cause a premature exit
try:
    f = open(output_path, 'w')
except Exception:
    raise ValueError('Error while opening output file')

if args.filter:
    try:
        g = open(filter_path, 'w')
    except Exception:
        raise ValueError('Error while opening file for filtered results')

# Default settings for TRENTO data
trento_on = True
trento_seed = 0
pythia_on = True
trento_file = args.file

# PYTHIA settings (default values)
eCM = args.eCM
pTHatMin = args.pTHatMin
pTHatMax = args.pTHatMax
seed = args.seed
bool_on_dict = {False: 'on', True: 'off'}
QCD = bool_on_dict[args.QCDoff]
QED = bool_on_dict[args.QEDoff]
quench = args.quench
nevt = args.nevt
SJpTmin = args.SJpTmin
SJradius = args.SJradius
output = args.output
filter_results = args.filter
filter_file = filter_path

# PYTHIA without QCD and QED does nothing
if args.QCDoff and args.QEDoff:
    pythia_on = False

# ------------------------- DATA LOADING ----------------------------- #

# Load/generate TRENTO data for 100,000 Au Au events
if trento_file == 'generate':
    proc = subprocess.Popen('trento Au Au {}'.format(nevt).split(),
                            stdout=subprocess.PIPE)
    trento_data = np.array([l.split() for l in proc.stdout], dtype=float)
elif trento_file == 'off':
    trento_on = False
else:
    try:
        trento_data = np.loadtxt(trento_file)
    except Exception as e:
        # At this point the file should exist
        # Being unable to open the file means it is being used
        # Or the user does not have permission to open it
        print('Unable to open TRENTO file.')
        raise e

# trento_data = [[event_number, impact_param, Npart, mult, e2, e3, e4, e5]]
# Conversion from fit to PHENIX 200 GeV data
mult_scale = 4.66

# -------------------------- PHYSICAL CONSTANTS ---------------------- #

# temperature in GeV
T = 0.15

# ???
deta = 4

# set pion mass to 0.14 GeV
mpi = 0.14

# rho_0 scaling parameter for radial flow
# from Retiere and Lisa, PRC70.044907 (2004), table 2
rho_0 = 0.85

# e2 scaling parameter for elliptic flow
# from Alver and Roland, PRC81.054905 (2010), fig 4
rho_2 = 0.15

# --------------------------- PYTHIA RUNS ---------------------------- #

# Initialize PYTHIA
pythia = pythia8.Pythia()

# Set PYTHIA output mode to quiet
pythia.readString('Print:quiet = on')

# Set eCM
pythia.readString('Beams:eCM = {}'.format(eCM))

# Set QCD
pythia.readString('HardQCD:all = {}'.format(QCD))

# Set QED
pythia.readString('PromptPhoton:all = {}'.format(QED))

# Set minimum transverse momentum (?)
pythia.readString('PhaseSpace:pTHatMin = {}'.format(pTHatMin))

# Set maximum transverse momentum (?)
pythia.readString('PhaseSpace:pTHatMax = {}'.format(pTHatMax))

# Enable seeding and set seed
pythia.readString('Random:setseed = on')
pythia.readString('Random:seed = {}'.format(seed))

# Initialize PYTHIA with given settings
pythia.init()

# Prepare output file
header = '''# Input parameters: to be listed here
#       trento-stat          |      pythia-truth                 | slowjet
# iev mult xgj quench b Np e2   part1 pT eta phi part2 pT eta phi   Njets pT eta phi pT eta phi
'''

f.write(header)
if filter_results:
    g.write(header)

# Loop over events, start with PYTHIA then add TRENTO to PYTHIA event
i = 0
problem_count = 0
while i < nevt:
    if trento_on:
        mult = int(mult_scale * trento_data[i][3])
        b = trento_data[i][1]
        Npart = int(trento_data[i][2])
        e2 = trento_data[i][4]
    else:
        mult = 0
        b = 0
        Npart = 0
        e2 = 0
    problem = False

    # Reset event data ???
    pythia.event.reset()
    if pythia_on:
        # Generate new PYTHIA event ???
        pythia.next()

        # for j in range(pythia.event.size()):
        #     print(pythia.event[j].list(False, False, 3))

        # Initial hard partons are stored in pythia.event [5] and [6] UNUSED
        # daughters5 = []
        # daughters5.extend(pythia.event[5].daughterList())
        # for j in daughters5:
        #     if j != 0:
        #         daughters5.extend(pythia.event[j].daughterList())
        # daughters6 = []

        # Pick pure and quenched jets
        # If there is a gamma jet, make that pure
        # Otherwise, pick one of the parton jets to be pure
        # Initial hard partons are stored in pythia.event [5] and [6]
        if pythia.event[5].id() == 22:
            purejet = pythia.event[5]
            quenchedjet = pythia.event[6]
        elif pythia.event[6].id() == 22:
            purejet = pythia.event[6]
            quenchedjet = pythia.event[5]
        else:
            if random.randint(1, 2) == 1:
                purejet = pythia.event[5]
                quenchedjet = pythia.event[6]
            else:
                purejet = pythia.event[6]
                quenchedjet = pythia.event[5]

        # Get PYTHIA event multiplicity
        pythia_mult = 0
        for j in range(pythia.event.size()):
            if pythia.event[j].isFinal():
                pythia_mult += 1

        # TODO introduce checks on hard-scattered particle assignments
        # First check 0 < abs(jet.id()) < 100
        for jet in (purejet, quenchedjet):
            if abs(jet.id()) < 1 or abs(jet.id()) > 99:
                problem = True
                print('Problem with jet id found: {}'.format(jet.id()))

        # Rescale quenchedjet daughters by momentum, not same as energy
        daughters = []
        daughters.extend(quenchedjet.daughterList())
        for j in daughters:
            if j != 0:
                daughters.extend(pythia.event[j].daughterList())

        for j in daughters:
            prt = pythia.event[j]
            px = quench * prt.px()
            py = quench * prt.py()
            pz = quench * prt.pz()
            prt_mass = prt.m()
            prt_e = (prt_mass**2 + px**2 + py**2 + pz**2)**0.5
            prt.px(px)
            prt.py(py)
            prt.pz(pz)
            prt.e(prt_e)

    # The TRENTO particle production below is un-physical,
    # put together by last year's summer student.
    # It will need to be re-written at some point.
    if trento_on:
        for j in range(mult):

            p_sample = 1
            while p_sample > 0.99:
                p_sample = np.random.random()

            particle_dict = {'pi+': (0.140, 211),
                             'pi-': (0.140, -211),
                             'pi0': (0.135, 111),
                             'K+': (0.494, 321),
                             'K-': (0.494, -321),
                             'p': (0.938, 2212),
                             'pbar': (0.938, -2212),
                             'n': (0.940, 2112),
                             'nbar': (0.940, -2112)}

            # Particle selected randomly (still needs to be changed)
            if p_sample <= 0.11:
                particle_name = 'pi+'
            elif 0.11 < p_sample <= 0.22:
                particle_name = 'pi-'
            elif 0.22 < p_sample <= 0.33:
                particle_name = 'pi0'
            elif 0.33 < p_sample <= 0.44:
                particle_name = 'K+'
            elif 0.44 < p_sample <= 0.55:
                particle_name = 'K-'
            elif 0.55 < p_sample <= 0.66:
                particle_name = 'p'
            elif 0.66 < p_sample <= 0.77:
                particle_name = 'pbar'
            elif 0.77 < p_sample <= 0.88:
                particle_name = 'n'
            elif 0.88 < p_sample <= 0.99:
                particle_name = 'nbar'
            else:
                raise Exception('Invalid particle value: {}'.format(p_sample))
            mass, pid = particle_dict[particle_name]

            r1, r2, r3 = np.random.random(3)

            # Compute 3-momentum magnitude
            pmag = -1 * T * math.log(r1 * r2 * r3)

            # Compute relativistic energy
            E = math.sqrt(pmag * pmag + mass * mass)

            # Check results against boltzmann distribution
            b_sample = np.random.random()
            while b_sample >= math.exp((pmag - E) / T):
                # Get new values until check is passed
                r1, r2, r3 = np.random.random(3)
                pmag = -1 * T * math.log(r1 * r2 * r3)
                E = math.sqrt(pmag * pmag + mass * mass)
                b_sample = np.random.random()

            # Compute theta and phi based on "Scott's Tick for Boltzmann 3D"
            ctheta = math.log(r1 / r2) / math.log(r1 * r2)
            stheta = math.sqrt(1 - ctheta * ctheta)
            phi = T * T * (math.log(r1 * r2)**2) / (pmag * pmag)
            if phi > 1.0:
                raise ValueError('phi={}, out of range'.format(phi))
            phi = 2.0 * math.pi * phi

            # Compute 3-mom components
            pz = pmag * ctheta
            px = pmag * stheta * math.cos(phi)
            py = pmag * stheta * math.sin(phi)

            # Append TRENTO event to PYTHIA event record
            # particle id, status, mothers, daughters, px, py, pz, E, m,
            # scaleIn, polIn
            pythia.event.append(pid, 200, 0, 0, px, py, pz, E, mass, 0.0, 9.0)

    # Initialize and call SlowJet
    etaMax = 4.
    nSel = 2
    massSet = 2
    slowJet = pythia8.SlowJet(-1, SJradius, SJpTmin, etaMax, nSel, massSet)
    slowJet.analyze(pythia.event)
    Njets = slowJet.sizeJet()
    if (Njets == 1):
        xgj = 0
    else:
        # xgj = slowJet.pT(1)/slowJet.pT(0)
        # pythia.event.list()
        # print('Jet 0')
        # for i in slowJet.constituents(0):
        #     print(i)
        #     # print('no: {} id: {} status: {}'.format(0, i.id(), i.status()))
        # print('Jet 1')
        # for i in slowJet.constituents(1):
        #     print(i)
        # #     print('no: {} id: {} status: {}'.format(0, i.id(), i.status()))
        if slowJet.constituents(0)[0] in purejet.daughterListRecursive():
            purejet_index = 0
            xgj = slowJet.pT(1)/slowJet.pT(0)
        elif slowJet.constituents(1)[0] in purejet.daughterListRecursive():
            purejet_index = 1
            xgj = slowJet.pT(0)/slowJet.pT(1)
        else:
            print('Neither jet was identified as pure')
            xgj = 0
            problem = True

    # Prepare output, stats include TRENTO information
    # (values set to zero if not turned on)
    stats_output = '{: 5d}'.format(i) \
                 + '{: 6d}'.format(mult + pythia_mult) \
                 + '{: 5.3f}'.format(xgj) \
                 + '{: 4.2f}'.format(quench) \
                 + '{: 6.2f}'.format(b) \
                 + '{: 4d}'.format(Npart) \
                 + '{: 5.3f}'.format(e2)

    if pythia_on:
        gamma_output = '{: 3d}'.format(purejet.id()) \
                     + '{: 5.3f}'.format(purejet.pT()) \
                     + '{: 8.3f}'.format(purejet.eta()) \
                     + '{: 8.3f}'.format(purejet.phi())
        quark_output = '{: 3d}'.format(quenchedjet.id()) \
                     + '{: 5.3f}'.format(quenchedjet.pT()) \
                     + '{: 8.3f}'.format(quenchedjet.eta()) \
                     + '{: 8.3f}'.format(quenchedjet.phi())
        pythia_output = gamma_output + quark_output
    else:
        pythia_output = '0 0 0 0 0 0 0 0'
        problem = True
        print('Pythia is not on for some reason')

    slowjet_output = '{: 3d}'.format(Njets)
    if Njets > 0:
        slowjet_output += ' {: 5.3f}'.format(slowJet.pT(0)) \
                        + ' {: 5.3f} '.format(slowJet.p(0).eta()) \
                        + ' {: 5.3f}'.format(slowJet.phi(0))
        if Njets > 1:
            slowjet_output += ' {: 5.3f}'.format(slowJet.pT(1)) \
                            + ' {: 5.3f} '.format(slowJet.p(1).eta()) \
                            + ' {: 5.3f}'.format(slowJet.phi(1))
        else:
            slowjet_output += ' 0 0 0'
            problem = True
            print('Second jet not found')
    else:
        slowjet_output += ' 0 0 0 0 0 0'
        problem = True
        print('No jets found')

    output = '{} {} {}\n'.format(stats_output, pythia_output,
                                 slowjet_output)
    if filter_results and problem:
        g.write(output)
        filename = '{}_{}_event.txt'.format(i, problem_count)
        filename = os.path.join(filtered_events_path, filename)
        _save_event(pythia.event, filename)
        problem_count += 1
    else:
        f.write(output)
        i += 1
        problem_count = 0

f.close()
if args.filter:
    g.close()
