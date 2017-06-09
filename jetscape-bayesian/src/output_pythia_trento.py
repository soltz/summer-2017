from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
import math
import pythia8
import numpy as np
import getopt
import argparse


def usage():
    output = '''
Generates output below for pythia jets in trento background with
slowjet reconstruction:

Usage: python output_pythia_trento_slowjet.py [options]

Writes output to output.txt, see header for format.

Options:
    -h, --help       : this message
    -f, --file       = set trento data file [../data/AuAu_200GeV_100k.txt]
                       or turn [off]
    -e, --eCM        = pythia beam center-of-mass energy (GeV) [200.0]
    -n, --nevt       = number of pythia + trento events to generate [10]
    -c, --QCDoff     : turn pythia hard QCD processes off
    -q, --QEDoff     : turn pythia hard QED processes off
    -o, --output     = output file [output.txt]
    -s, --seed       = pythia initial random number seed [-1]
    -p, --SJpTmin    = slowjet minimum pT [10]
    -r, --SJradius   = slowJet radius [0.5]
    -u, --quench     = scaling factor for momentum of non-photon jet,
                       QED only [1.0]
    -y, --pTHatMin   = pythia minimum jet pT [20.0]
    -z, --pTHatMax   = pythia maximum jet pT [25.0]
    '''
    print(output)


# Parse command line and set defaults
#   (see http://docs.python.org/library/getopt.html)
# If unrecognized option is passed, a GetoptError will be raised and
# caught. The error will be shown as 'option -a not recognized',
# followed by the proper usage. Then the program exits.
try:
    opts, args = getopt.getopt(sys.argv[1:], 'hf:e:n:cqo:s:p:r:u:y:z:l',
                               ['help', 'file=', 'eCM=', 'nevt=',
                                'pTHatMax=', 'pTHatMin=', 'output=',
                                'seed=', 'QCD', 'QED', 'quench=',
                                'SJpTmin=', 'SJradius='])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

# Default settings for Trento data
trento_on = True
trento_seed = 0
pythia_on = True
trento_file = '../data/AuAu_200GeV_100k.txt'

# PYTHIA settings (default values)
eCM = 200.0
pTHatMin = 20.0
pTHatMax = 25.0
seed = -1
QCD = 'on'
QED = 'on'
quench = 1.0
nevt = 10
SJpTmin = 10
SJradius = 0.5
output = 'output.txt'

# Iterate over passed options and override defaults where applicable
for o, a in opts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-f', '--file'):
        trento_file = str(a)
    elif o in ('-e', '--eCM'):
        eCM = float(a)
    elif o in ('-n', '--nevt'):
        nevt = int(a)
    elif o in ('-y', '--pTHatMin'):
        pTHatMin = float(a)
    elif o in ('-z', '--pTHatMax'):
        pTHatMax = float(a)
    elif o in ('-o', '--output'):
        output = a
    elif o in ('-s', '--seed'):
        seed = int(a)
    elif o in ('-c', '--QCDoff'):
        QCD = 'off'
    elif o in ('-q', '--QEDoff'):
        QED = 'off'
    elif o in ('-u', '--quench'):
        quench = float(a)
    elif o in ('-p', '--SJpTmin'):
        SJpTmin = float(a)
    elif o in ('-r', '--SJradius'):
        SJradius = float(a)
    else:
        raise Exception('unhandled option')

# PYTHIA without QCD and QED does nothing ???
if (QCD == 'off') and (QED == 'off'):
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
