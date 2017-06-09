from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
import math
import pythia8
import numpy as np
import argparse


prog_description = '''
Generates output below for pythia jets in trento background with
slowjet reconstruction.'''.replace('\n', '', 1)
prog_epilogue = 'Writes output to output.txt, see header for format.'

# Initialize argument parser
parser = argparse.ArgumentParser(
    description=prog_description,
    epilog=prog_epilogue,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add arguments
parser.add_argument('-f', '--file', type=str,
                    default='../data/AuAu_200GeV_100k.txt',
                    help='path to trento data file or turn \"off\"')
parser.add_argument('-e', '--eCM', type=float, default=200.0,
                    help='PYTHIA beam center-of-mass energy (in GeV)')
parser.add_argument('-n', '--nevt', type=int, default=10,
                    help='number of PYTHIA and trento events to generate')
parser.add_argument('-c', '--QCDoff', action='store_true',
                    help='turn PYTHIA hard QCD processes off')
parser.add_argument('-q', '--QEDoff', action='store_true',
                    help='turn PYTHIA hard QED processes off')
parser.add_argument('-o', '--output', type=str, default='output.txt',
                    help='path to output file')
parser.add_argument('-s', '--seed', type=int, default=-1,
                    help='seed for random numbers in PYTHIA')
parser.add_argument('-p', '--SJpTmin', type=int, default=10,
                    help='slowjet minimum pT')
parser.add_argument('-r', '--SJradius', type=float, default=0.5,
                    help='slowjet radius')
parser.add_argument('-u', '--quench', type=float, default=1.0,
                    help='scaling factor for momentum of non-photon jet')
parser.add_argument('-y', '--pTHatMin', type=float, default=20.0,
                    help='PYTHIA minimum jet pT')
parser.add_argument('-z', '--pTHatMax', type=float, default=25.0,
                    help='PYTHIA maximum jet pT')

# Parse arguments
args = parser.parse_args()

# Check validity of arguments
if args.file != 'off':
    try:
        with open(args.file, 'r') as f:
            pass
    except FileNotFoundError:
        raise ValueError('trento file {} not found.'.format(args.file))

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

# Check last to prevent changing files on the system when other
# arguments might still cause a premature exit
try:
    f = open(args.output, 'w')
except Exception:
    raise ValueError('Error while opening output file.')

# Default settings for Trento data
trento_on = True
trento_seed = 0
pythia_on = True
trento_file = args.file

# PYTHIA settings (default values)
eCM = args.eCM
pTHatMin = args.pTHatMin
pTHatMax = args.pTHatMax
seed = args.seed
QCD = not args.QCDoff
QED = not args.QEDoff
quench = args.quench
nevt = args.nevt
SJpTmin = args.SJpTmin
SJradius = args.SJradius
output = args.output

# PYTHIA without QCD and QED does nothing ???
if args.QCDoff and args.QEDoff:
    pythia_on = False

# Load trento data for 100,000 Au Au events
if trento_file == 'off':
    trento_on = False
    b = np.zeros(1e5)
    Npart = np.zeros(1e5)
    mult = np.zeros(1e5)
    e2 = np.zeros(1e5)
else:
    try:
        data = np.loadtxt(trento_file)
    except FileNotFoundError:
        print('\nFile ' + trento_file + ' not found, exiting...' + '\n')
        sys.exit()

    # Conversion from fit to PHENIX 200 GeV data
    mult_scale = 4.66

    # data = [[event_number, impact_param, Npart, mult, e2, e3, e4, e5]]
    b = data[:, 1]
    Npart = data[:, 2]
    mult = mult_scale * data[:, 3]
    e2 = data[:, 4]

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

# Initialize Pythia
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

# Prepare output file ???
header = '''# Input parameters: to be listed here
#       trento-stat          |      pythia-truth                 | slowjet
# iev mult xgj quench b Np e2   part1 pT eta phi part2 pT eta phi   Njets pT eta phi pT eta phi
'''

f = open(output, 'w')
f.write(header)

# Loop over events, start with pythia then add trento to pythia event
for i in range(trento_seed, trento_seed + nevt):

    # Reset event data ???
    pythia.event.reset()
    if pythia_on:
        # Generate new PYTHIA event ???
        pythia.next()

        # Initial hard partons are stored in pythia.event [5] and [6]
        daughters5 = []
        daughters5.extend(pythia.event[5].daughterList())
        for j in daughters5:
            if j != 0:
                daughters5.extend(pythia.event[j].daughterList())
        daughters6 = []

        # Check for gamma in [5] position
        # For QCD events quarkjet is always [5]
        if pythia.event[5].name() == 'gamma':
            gammajet = pythia.event[5]
            quarkjet = pythia.event[6]
        else:
            quarkjet = pythia.event[5]
            gammajet = pythia.event[6]

        if gammajet.id() != 22:
            print('Danger, Will Robinson!')

        # TODO introduce checks on hard-scattered particle assignments

        # Rescale daughters by momentum, not same as energy
        daughters = []
        daughters.extend(quarkjet.daughterList())
        for j in daughters:
            if j != 0:
                daughters.extend(pythia.event[j].daughterList())

        for j in daughters:
            prt = pythia.event[j]
            px = quench*prt.px()
            py = quench*prt.py()
            pz = quench*prt.pz()
            prt_mass = prt.m()
            prt_e = (prt_mass**2 + px**2 + py**2 + pz**2)**0.5
            prt.px(px)
            prt.py(py)
            prt.pz(pz)
            prt.e(prt_e)

    # The Trento particle production below is un-physical,
    # put together by last year's summer student.
    # It will need to be re-written at some point.
    if trento_on:
        for j in range(int(mult[i])):
            r1, r2, r3, r4, r5 = np.random.random(5)
            while r5 > 0.99:
                r5 = np.random.random(1)

            # pT = transverse momentum
            pT_r1 = T*(math.sqrt(-2*math.log(r1)))

            # phi = azimuthal angle
            phi_r2 = 2*(math.pi)*(r2 - 0.5)

            # eta = pseudo-rapidity
            eta_r3 = deta*(r3 - 0.5)

            # rho = normalized radial distance
            rho_r4 = r4**0.5

            particle_dict = {'pi+': (0.140, 211),
                             'pi-': (0.140, -211),
                             'pi0': (0.135, 111),
                             'K+': (0.494, 321),
                             'K-': (0.494, -321),
                             'p': (0.938, 2212),
                             'pbar': (0.938, -2212),
                             'n': (0.940, 2112),
                             'nbar': (0.940, -2112)}
            # particle selected randomly
            if r5 <= 0.11:
                particle_name = 'pi+'
            elif 0.11 < r5 <= 0.22:
                particle_name = 'pi-'
            elif 0.22 < r5 <= 0.33:
                particle_name = 'pi0'
            elif 0.33 < r5 <= 0.44:
                particle_name = 'K+'
            elif 0.44 < r5 <= 0.55:
                particle_name = 'K-'
            elif 0.55 < r5 <= 0.66:
                particle_name = 'p'
            elif 0.66 < r5 <= 0.77:
                particle_name = 'pbar'
            elif 0.77 < r5 <= 0.88:
                particle_name = 'n'
            elif 0.88 < r5 <= 0.99:
                particle_name = 'nbar'
            else:
                raise Exception('Invalid particle value: {}'.format(r5))
            mass, pid = particle_dict[particle_name]

            # calculate initial transverse rapidity (yT)
            eT = (mass*mass+pT_r1*pT_r1)**0.5
            yT = 0.5 * np.log((eT+pT_r1)/(eT-pT_r1))
            pT_initial = pT_r1
            yT_initial = yT

            # apply flow as additive boost to transverse rapidity
            yBoost = rho_r4*rho_0 + rho_2*e2[i]*np.cos(2*phi_r2)
            yT = yT_initial + yBoost

            # convert back to pT
            pT_wflow = mass*np.cosh(yT)

            # add particles to the pythia event list
            px = pT_wflow * math.cos(phi_r2)
            py = pT_wflow * math.sin(phi_r2)
            pz = pT_wflow * math.sinh(eta_r3)
            E = (pT_wflow**2 + pz**2 + mass**2)**0.5
            pythia.event.append(pid, 200, 0, 0, px, py, pz, E, mass, 0.0,
                                9.0)

    # Initialize and call SlowJet
    etaMax = 4.
    nSel = 2
    massSet = 2
    slowJet = pythia8.SlowJet(-1, SJradius, SJpTmin, etaMax, nSel, massSet)
    slowJet.analyze(pythia.event)
    Njets = slowJet.sizeJet()
    if (Njets == 1):
        xgj = 0.
    else:
        xgj = slowJet.pT(1)/slowJet.pT(0)

    # Prepare output, stats include trento information
    # (values set to zero if not turned on)
    stats_output = '{: 5d}'.format(i) \
                 + '{: 6d}'.format(int(mult[i])) \
                 + '{: 5.3f}'.format(xgj) \
                 + '{: 4.2f}'.format(quench) \
                 + '{: 6.2f}'.format(b[i]) \
                 + '{: 4d}'.format(int(Npart[i])) \
                 + '{: 5.3f}'.format(e2[i])

    if pythia_on:
        gamma_output = '{: 3d}'.format(gammajet.id()) \
                     + '{: 5.3f}'.format(gammajet.pT()) \
                     + '{: 8.3f}'.format(gammajet.eta()) \
                     + '{: 8.3f}'.format(gammajet.phi())
        quark_output = '{: 3d}'.format(quarkjet.id()) \
                     + '{: 5.3f}'.format(quarkjet.pT()) \
                     + '{: 8.3f}'.format(quarkjet.eta()) \
                     + '{: 8.3f}'.format(quarkjet.phi())
        pythia_output = gamma_output + quark_output
    else:
        pythia_output = '0 0 0 0 0 0 0 0'

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
    else:
        slowjet_output += ' 0 0 0 0 0 0'

    output = '{} {} {}\n'.format(stats_output, pythia_output,
                                 slowjet_output)
    f.write(output)
