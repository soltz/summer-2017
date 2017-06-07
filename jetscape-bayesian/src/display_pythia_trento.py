import math
import pythia8
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import getopt, sys

def usage():
    print 'Shows the frequency of slowJet jet pT when the pT of pythia events is restricted to 20-25 GeV/c.'
    print 'Usage: python restricted_jetpT.py [options]'
    print '   -h, --help      : this message'
    print '   -t, --trento     : turn trento background off'
    print '   -d, --trentoseed     = change first trento event that is viewed[0]'
    print '   -o, --pythia     : turn pythia events off'
    print '   -f, --file     = set trento data file [AuAu_200GeV_100k.txt]'
    print '   -e, --eCM     = pythia beam center-of-mass energy (GeV) [200.0]'
    print '   -n, --pTHatMin     = pythia minimum jet pT [20.0]'
    print '   -x, --pTHatMax     = pythia maximum jet pT [25.0]'
    print '   -s, --seed     = pythia initial random number seed [-1]'
    print '   -c, --QCD     : turn pythia hard QCD processes off'
    print '   -q, --QED     : turn pythia hard QED processes on'
    print '   -u, --quench     = scaling factor for momentum of non-photon jet, QED only [1.0]'
    print '   -p, --pTjetMin     = minimum slowJet pT [15]'
    print '   -r, --radius     = slowJet radius [0.7]'
    print '   -R, --radopt     : require slowjet to find only 2 jets by varying radius'
    print '   -P, --pTopt     : require slowjet to find only 2 jets by varying pTjetMin'
    print '   -b, --bins     = number of histogram bins on each axis [20]'
    print '   -l, --labels     : turn plot color labels off'

def main():

#   Parse command line and set defaults (see http://docs.python.org/library/getopt.html)
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'htd:of:e:n:x:s:cqu:p:r:RPb:l', \
              ['help','trento','trentoseed=','pythia','file=','eCM=','pTHatMin=','pTHatMax=','seed=','QCD','QED','quench=','pTjetMin=','radius=','radopt','pTopt','bins=','labels'])
    except getopt.GetoptError, err:
        print str(err) # will print something like 'option -a not recognized'
        usage()
        sys.exit(2)

    trento = True
    trento_seed = 0
    pythia_on = True
    trento_file = '../../../data/AuAu_200GeV_100k.txt'
    
    # pythia settings
    eCM  = 200.0
    pTHatMin  = 20.0
    pTHatMax  = 25.0
    seed  = -1
    QCD = 'on'
    QED = 'off'
    quench = 1.0

    # slowJet settings
    pTjetMin_in = 15
    radius_in = 0.7
    pT_opt = False
    rad_opt = False

    # plot settings
    bins = 20
    labels = True

    for o, a in opts:
        if o in ('-h', '--help'):
            usage()
            sys.exit()
        elif o in ('-t', '--trento'):
            trento = False
        elif o in ('-d', '--trentoseed'):
            trento_seed = int(a)
        elif o in ('-o', '--pythia'):
            pythia_on = False
        elif o in ('-f', '--file'):
            trento_file = str(a)
        elif o in ('-e', '--eCM'):
            eCM = float(a)
        elif o in ('-n', '--pTHatMin'):
            pTHatMin = float(a)
        elif o in ('-x', '--pTHatMax'):
            pTHatMax = float(a)
        elif o in ('-s', '--seed'):
            seed = int(a)
        elif o in ('-c', '--QCD'):
            QCD = 'off'
        elif o in ('-q', '--QED'):
            QED = 'on'
        elif o in ('-u', '--quench'):
            quench = float(a)
        elif o in ('-p', '--pTjetMin'):
            pTjetMin_in = float(a)
        elif o in ('-r', '--radius'):
            radius_in = float(a)
        elif o in ('-R', '--radopt'):
            rad_opt = True
        elif o in ('-P', '--pTopt'):
            pT_opt = True
        elif o in ('-b', '--bins'):
            bins = int(a)
        elif o in ('-l', '--labels'):
            labels = False
        else:
            assert False, 'unhandled option'
            
    #   load trento data for 100,000 Au Au events
    data = np.loadtxt(trento_file)

    if trento:
        #   data = [[event_number, impact_param, Npart, mult, e2, e3, e4, e5],...]
        #   create a list for the initial entropy of the events
        mult = data[:,3]
        e2 = data[:,4]
        
        #   apply proportionality constanst to convert initial entropy to charged multiplicity
        #   constant was calculated by fitting the trento data to phenix data
        fit_par = 4.65905256
        if trento_file == 'AuAu_200GeV_100k.txt':
            fit_par = 4.65905256
        if trento_file == 'AuAu_130GeV_100k.txt':
            fit_par = 3.96639009
        if trento_file == 'AuAu_62p4GeV_100k.txt':
            fit_par = 3.13381932
        if trento_file == 'AuAu_39GeV_100k.txt':
            fit_par = 2.54804101
        if trento_file == 'AuAu_27GeV_100k.txt':
            fit_par = 2.16041632
        if trento_file == 'AuAu_19p6GeV_100k.txt':
            fit_par = 1.95015905
        if trento_file == 'AuAu_15p0GeV_100k.txt':
            fit_par = 1.6460897
        if trento_file == 'AuAu_7p7GeV_100k.txt':
            fit_par = 1.27382896
        mult = fit_par * mult
        
        #   change the constituents of mult to int values
        for i in range(len(mult)):
            mult[i] = round(mult[i])
        mult = mult.astype(np.int64)
        
    T = 0.15 # temperature in GeV
    deta = 4
    # set pion mass to 0.14 GeV
    mpi = 0.14
    
    # rho_0 scaling parameter for radial flow
    # from Retiere and Lisa, PRC70.044907 (2004), table 2
    rho_0 = 0.85
    
    # e2 scaling parameter for elliptic flow
    # from Alver and Roland, PRC81.054905 (2010), fig 4
    rho_2 = 0.15
    
    #   Initialize Pythia
    pythia = pythia8.Pythia()
        
    eCM = str(eCM)
    set_eCM = "Beams:eCM = " + eCM
    pythia.readString(set_eCM)
    
    set_QCD = "HardQCD:all = " + QCD
    pythia.readString(set_QCD)
    
    set_QED = "PromptPhoton:all = " + QED
    pythia.readString(set_QED)
    
    pTHatMin = str(pTHatMin)
    set_pTHatMin = "PhaseSpace:pTHatMin = " + pTHatMin
    pythia.readString(set_pTHatMin)

    pTHatMax = str(pTHatMax)
    set_pTHatMax = "PhaseSpace:pTHatMax = " + pTHatMax
    pythia.readString(set_pTHatMax)
    
    pythia.readString("Random:setseed = on")
        
    seed = str(seed)
    set_seed = "Random:seed = " + seed
    pythia.readString(set_seed)
    
    pythia.init()
    
    for i in range(trento_seed,100000): 
    
        pT = []
        phi = []
        eta = []
    
    #   empty lists for pythia particle info
        pjet_eta = []
        pjet_phi = []
        pjet_pT = []
        pjet_eT = []
            
        pbg_eta = []
        pbg_phi = []
        pbg_pT = []
        pbg_eT = []
    
    #   empty lists for trento particle info
        tjet_eta = []
        tjet_phi = []
        tjet_pT = []
        tjet_eT = []
            
        tbg_eta = []
        tbg_phi = []
        tbg_pT = []
        tbg_eT = []
        
        pythia.event.reset()
        if pythia_on:
            pythia.next()

        #   The daughters of the initial hard process are recorded below
            if QCD == 'on':
                daughters5 = []
                daughters5.extend(pythia.event[5].daughterList())
                for j in daughters5:
                   if j != 0:
                        daughters5.extend(pythia.event[j].daughterList())
                daughters6 = []
                daughters6.extend(pythia.event[6].daughterList())
                for j in daughters6:
                   if j != 0:
                        daughters6.extend(pythia.event[j].daughterList())
                daughters5.extend(daughters6)
                daughters = daughters5
                
            if QED == 'on':
                daughters5 = []
                daughters5.extend(pythia.event[5].daughterList())
                for j in daughters5:
                   if j != 0:
                        daughters5.extend(pythia.event[j].daughterList())
                daughters6 = []
                daughters6.extend(pythia.event[6].daughterList())
                for j in daughters6:
                   if j != 0:
                        daughters6.extend(pythia.event[j].daughterList())
                        
                # Non-photon QED jets are rescaled
                if len(daughters5) > len(daughters6):
                    for j in daughters5:
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
                if len(daughters6) > len(daughters5):
                    for j in daughters6:
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

                daughters5.extend(daughters6)
                daughters = daughters5

        if trento:
            for j in range(mult[i]):
                r1, r2, r3, r4, r5= np.random.random(5)
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

                # particle selected randomly
                if r5 <= 0.11:
                    mass = 0.140
                    pid = 211 # pi+
                if r5 > 0.11 and r5 <= 0.22:
                    mass = 0.140
                    pid = -211 # pi-
                if r5 > 0.22 and r5 <= 0.33:
                    mass = 0.135
                    pid = 111 # pi0
                if r5 > 0.33 and r5 <= 0.44:
                    mass = 0.494
                    pid = 321 # K+
                if r5 > 0.44 and r5 <= 0.55:
                    mass = 0.494
                    pid = -321 # K-
                if r5 > 0.55 and r5 <= 0.66:
                    mass = 0.938
                    pid = 2212 # p
                if r5 > 0.66 and r5 <= 0.77:
                    mass = 0.938
                    pid = -2212 # pbar
                if r5 > 0.77 and r5 <= 0.88:
                    mass = 0.940
                    pid = 2112 # n
                if r5 > 0.88 and r5 <= 0.99:
                    mass = 0.940
                    pid = -2112 # nbar
        
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
        
                pT.append(pT_wflow)
                phi.append(phi_r2)
                eta.append(eta_r3)
        
                # add particles to the event list
                px = pT_wflow * math.cos(phi_r2)
                py = pT_wflow * math.sin(phi_r2)
                pz = pT_wflow * math.sinh(eta_r3)
                E = (pT_wflow**2 + pz**2 + mass**2)**0.5
                pythia.event.append(pid, 200, 0, 0, px, py, pz, E, mass, 0., 9.)
        #   Initialize SlowJet
        etaMax = 4.
        nSel = 2
        massSet = 2
        radius = radius_in
        pTjetMin = pTjetMin_in

        slowJet = pythia8.SlowJet( -1, radius, pTjetMin, etaMax, nSel, massSet);
        slowJet.analyze(pythia.event)
        jets_found = slowJet.sizeJet()

        #   Option to raise or lower radius until only 2 jets are found
        if rad_opt and not pT_opt:
            while jets_found != 2 and radius > 0 and radius < 1.0:
                slowJet = pythia8.SlowJet( -1, radius, pTjetMin, etaMax, nSel, massSet);
                slowJet.analyze(pythia.event)
                jets_found = slowJet.sizeJet()
            
                if jets_found > 10:
                    radius = radius - 0.1
                elif jets_found > 3:
                    radius = radius - 0.01
                elif jets_found > 2:
                    radius = radius - 0.001
                elif jets_found < 2:
                    radius = radius + 0.01

        #   Option to raise or lower pTjetMin until only 2 jets are found
        if pT_opt and not rad_opt:
            while jets_found != 2 and pTjetMin > 0:
                slowJet = pythia8.SlowJet( -1, radius, pTjetMin, etaMax, nSel, massSet);
                slowJet.analyze(pythia.event)
                jets_found = slowJet.sizeJet()
           
                if jets_found > 10:
                    pTjetMin = pTjetMin + 1.0
                elif jets_found > 6:
                    pTjetMin = pTjetMin + 0.1
                elif jets_found > 2:
                    pTjetMin = pTjetMin + 0.01
                elif jets_found < 2:
                    pTjetMin = pTjetMin - 0.01

        #   Option to raise or lower pTjetMin until only 2 jets are found
        #   If pTjetMin exceeds 40, option will move to radius
        if pT_opt and rad_opt:
            while jets_found != 2 and pTjetMin < 40.0 and pTjetMin > 0:
                slowJet = pythia8.SlowJet( -1, radius, pTjetMin, etaMax, nSel, massSet);
                slowJet.analyze(pythia.event)
                jets_found = slowJet.sizeJet()
           
                if jets_found > 10:
                    pTjetMin = pTjetMin + 1.0
                elif jets_found > 6:
                    pTjetMin = pTjetMin + 0.1
                elif jets_found > 2:
                    pTjetMin = pTjetMin + 0.01
                elif jets_found < 2:
                    pTjetMin = pTjetMin - 0.01
                    
            while jets_found != 2 and radius > 0 and radius < 1:
                slowJet = pythia8.SlowJet( -1, radius, pTjetMin, etaMax, nSel, massSet);
                slowJet.analyze(pythia.event)
                jets_found = slowJet.sizeJet()
            
                if jets_found > 10:
                    radius = radius - 0.1
                elif jets_found > 3:
                    radius = radius - 0.01
                elif jets_found > 2:
                    radius = radius - 0.001
                elif jets_found < 2:
                    radius = radius + 0.01

    
    #   Extract constituents and convert nested list for ease of manipulation
        slowJetPrtList = [[] for j in range(slowJet.sizeJet())]
        for j in range(slowJet.sizeJet()):
            slowJetPrtList[j] = list(slowJet.constituents(j))
    
        if slowJet.sizeJet() > 0:
            totJetPrtList = np.concatenate(slowJetPrtList)
        else: totJetPrtList = [0]

        if pythia_on:                            
            for j in range(pythia.event.size()):
                prt = pythia.event[j]
                if prt.isFinal():
                    prt_eta = prt.eta()
                    prt_phi = prt.phi()
                    prt_pT = prt.pT()
                    prt_eT = prt.eT()
                    if j in totJetPrtList and j in daughters:
                        pjet_eta.append(prt_eta)
                        pjet_phi.append(prt_phi)
                        pjet_pT.append(prt_pT)
                        pjet_eT.append(prt_eT)
                    elif j in totJetPrtList:
                        tjet_eta.append(prt_eta)
                        tjet_phi.append(prt_phi)
                        tjet_pT.append(prt_pT)
                        tjet_eT.append(prt_eT)
                    elif j in daughters:
                        pbg_eta.append(prt_eta)
                        pbg_phi.append(prt_phi)
                        pbg_pT.append(prt_pT)
                        pbg_eT.append(prt_eT)
                    else:
                        tbg_eta.append(prt_eta)
                        tbg_phi.append(prt_phi)
                        tbg_pT.append(prt_pT)
                        tbg_eT.append(prt_eT)

    #   Create bins using np.histogram2d(), the number of bins produced is the square of the value "bins"
        tot_bins = bins**2
    
    #   Data from trento background is binned
        tbg_h, xedges, yedges = np.histogram2d(tbg_eta,tbg_phi,bins=bins,range=[[-2,2], [-(np.pi), np.pi]],weights=tbg_eT)
        tbg_h = np.concatenate(tbg_h)
    
    #   Data from trento particles identified as jets is binned    
        tjet_h, tjet_xedges, tjet_yedges = np.histogram2d(tjet_eta,tjet_phi,bins=[xedges,yedges],weights=tjet_eT)
        tjet_h = np.concatenate(tjet_h)
    
    #   Data from pythia background is binned
        pbg_h, pbg_xedges, pbg_yedges = np.histogram2d(pbg_eta,pbg_phi,bins=[xedges,yedges],weights=pbg_eT)
        pbg_h = np.concatenate(pbg_h)
    
    #   Data from pythia particles identified as jets is binned
        pjet_h, pjet_xedges, pjet_yedges = np.histogram2d(pjet_eta,pjet_phi,bins=[xedges,yedges],weights=pjet_eT)
        pjet_h = np.concatenate(pjet_h)
    
    #   x and y coordinates for bincenters are created using bin edges
        xcenters = 0.5*(xedges[1:]+xedges[:-1])
        xcenters = np.repeat(xcenters,bins)
        ycenters = 0.5*(yedges[1:]+yedges[:-1])
        ycenters = np.tile(ycenters,bins)
    
    #   arrays of bin centers are extended to match the length of the data sets
        xcenters = list(xcenters)
        ycenters = list(ycenters)
        xcenters.extend(xcenters)
        ycenters.extend(ycenters)
        xcenters.extend(xcenters)
        ycenters.extend(ycenters)
            
    #   Bars associated with jets will sit on top of bars associated with background
    #   the order bars will be stacked is trento bg, pythia bg, trento jet, pythia jet     
    
    #   the base for trento bg will be all zeroes
        tbg_base = np.zeros(len(tjet_h))
        tbg_base = list(tbg_base)
    
    #   the base for pythia bg will be the top of the trento bg
        pbg_base = tbg_h
        pbg_base = list(pbg_base)
        tbg_base.extend(pbg_base)
        pbg_base = np.array(pbg_base)
    
    #   the base for trento jets will be the top of trento bg and pythia bg
        tjet_base = pbg_base + pbg_h
        tjet_base = list(tjet_base)
        tbg_base.extend(tjet_base)
        tjet_base = np.array(tjet_base)
    
    #   the base for pythia jets will be on top of trento bg, pythia bg, and trento jets
        pjet_base = tjet_base + tjet_h
        pjet_base = list(pjet_base)
        tbg_base.extend(pjet_base)
    
        base = tbg_base
            
        tbg_h = list(tbg_h)
        pbg_h = list(pbg_h)
        tjet_h = list(tjet_h)
        pjet_h = list(pjet_h)
        tbg_h.extend(pbg_h)
        tbg_h.extend(tjet_h)
        tbg_h.extend(pjet_h)
    
        h = tbg_h
            
    #   Set size of bars in plot
        d_eta = (xedges[1]-xedges[0])
        d_phi = (yedges[1]-yedges[0])
    
    #   Create an array of colors where bars associated with jet particles are red
    #   Background bars are blue and bars with no value are white
        jet_colors=[]
        
        tbg_col = 'blue'
        pbg_col = 'red'
        tjet_col = 'yellow'
        pjet_col = 'green'
        
        for j in h[0:tot_bins]:
            if j > 0:
                jet_colors.append(tbg_col)
            else:
                jet_colors.append('w')
        for j in range(len(pbg_h)):
            if pbg_h[j] > 0:
                jet_colors.append(pbg_col)
            elif h[j] > 0:
                jet_colors.append(tbg_col)
            else:
                jet_colors.append('w')
        for j in range(len(tjet_h)):
            if tjet_h[j] > 0:
                jet_colors.append(tjet_col)
            elif pbg_h[j] > 0:
                jet_colors.append(pbg_col)
            elif h[j] > 0:
                jet_colors.append(tbg_col)
            else:
                jet_colors.append('w')
        for j in range(len(pjet_h)):
            if pjet_h[j] > 0:
                jet_colors.append(pjet_col)
            elif tjet_h[j] > 0:
                jet_colors.append(tjet_col)
            elif pbg_h[j] > 0:
                jet_colors.append(pbg_col)
            elif h[j] > 0:
                jet_colors.append(tbg_col)
            else:
                jet_colors.append('w')
    
    #   End of event loop. Statistics. Histogram. Done.
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.bar3d(xcenters,ycenters,base,d_eta,d_phi,h,color=jet_colors,alpha=1)
    
        ax.set_xlabel('\n$\eta$',fontsize=20)
        ax.set_ylabel('\n$\phi$',fontsize=20)
        ax.set_zlabel('eT (GeV)')

        if trento:
            title = '\nslowJet in pythia with trento background' + '\ntrento multiplicity = ' + str(mult[i])
        if not pythia_on:
            title = '\nslowJet in trento background without pythia' + '\ntrento multiplicity = ' + str(mult[i])
        if not trento:
            title = 'slowJet in pythia without trento background'
        plt.title(title, loc = 'left')
        
        pi = np.pi
        plt.ylim(-pi,pi)

        m = max(h)
        index = [j for j, k in enumerate(h) if k == m]

        if pythia_on:
            if slowJet.sizeJet() >= 2:
                for j in range(2):
                    pT_jet = slowJet.pT(j)
                    pT_label = 'slowJet pT = %.2f' % pT_jet
                    ax.text(slowJet.y(j), slowJet.phi(j), max(h) + base[index[0]] +1, pT_label, horizontalalignment='center')
            else:
                for j in range(slowJet.sizeJet()):
                    pT_jet = slowJet.pT(j)
                    pT_label = 'slowJet pT = %.2f' % pT_jet
                    ax.text(slowJet.y(j), slowJet.phi(j), max(h) + base[index[0]] + 1, pT_label, horizontalalignment='center')

        if labels:
            blue_proxy = plt.Rectangle((0, 0), 1, 1, fc="b")
            red_proxy = plt.Rectangle((0, 0), 1, 1, fc="r")
            yellow_proxy = plt.Rectangle((0, 0), 1, 1, fc="y")
            green_proxy = plt.Rectangle((0, 0), 1, 1, fc="g")
            ax.legend([green_proxy,yellow_proxy,red_proxy,blue_proxy],['true jets','false jets','missed jets','background'])
        plt.show(block = False)
    
        pythia.event.list()
        print "radius: ", radius
        print "pTjetMin: ", pTjetMin
        print "jets found: ", slowJet.sizeJet()
        if trento:
            print "trento multiplicity: ", mult[i]
        
        query = raw_input("q to quit, p to save to png,  <CR> to continue: ")
        if (query=='q'):
            break
        if (query=='p'):
            query = raw_input("name this png file (do not include extension .png): ")
            filename = query + '.png'
            plt.savefig(filename)
            query = raw_input("q to quit, <CR> to continue: ")
            if (query=='q'):
                break
        else: continue
    
if __name__ == '__main__':main()
