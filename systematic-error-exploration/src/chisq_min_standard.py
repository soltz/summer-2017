from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from hepdata_manager import HEPData_Manager
import numpy as np
import iminuit
from scipy import stats


# Load data in
alice_event = 'ins1394678'
data_manager = HEPData_Manager()
if not data_manager.has_data(alice_event):
    data_manager.get_data(alice_event)
    print('Downloaded and saved {} data!'.format(alice_event))

table1_data = data_manager.load_table(alice_event, 'Table1')

# Global variables for use in minimization
global p, t, d, err_s, err_c, err_a

# Functions to be potentially minimized
def chi_squ_one_err(a):
    chisq = 0

    for i in range(len(p)):
        chisq += (d[i] - t[i] - a * err_c[i])**2 / err_s[i]**2

    chisq += a**2

    return chisq

def chi_squ_two_err(a, b):
    chisq = 0

    for i in range(len(p)):
        chisq += (d[i] - t[i] - a * err_c[i] - b * err_a[i])**2 / err_s[i]**2

    chisq += a**2
    chisq += b**2

    return chisq

def redmer_chi_sq(a, b):
    chisq = 0

    for i in range(len(p)):
        chisq += (d[i] + a * err_c[i] + b)**2 / err_s[i]**2

    chisq += a**2

    subsum = 0

    for i in range(len(p)):
        subsum += b**2 / err_a[i]**2
    subsum = subsum / len(p)

    chisq += subsum

    return chisq

def redmer_chi_sq_alt(a, b):
    err_n = 0
    for i in range(len(err_a)):
        err_n += 1 / err_a[i]**2
    err_n = err_n / len(err_a)
    err_n = np.sqrt(1 / err_n)

    chisq = 0

    for i in range(len(indep)):
        chisq += (dep[i] + a * corr_err[i] + b * err_n)**2 / stat_err[i]**2

    chisq += a**2
    chisq += b**2

    return chisq

indep = [(x['high'] + x['low'])/2
         for x in table1_data['independent_variables'][0]['values']]

theory = [0 for x in indep]

dep = [x['value'] for x in table1_data['dependent_variables'][0]['values']]

stat_err = [x['errors'][0]['symerror']
            for x in table1_data['dependent_variables'][0]['values']]

corr_err = [x['errors'][2]['symerror']
            for x in table1_data['dependent_variables'][0]['values']]

acorr_err = [x['errors'][1]['symerror']
             for x in table1_data['dependent_variables'][0]['values']]

# --- Full range ---

# Set globals before minimization
p = indep
t = theory
d = dep
err_s = stat_err
err_c = corr_err
err_a = acorr_err

# Minimize and print p
m = iminuit.Minuit(redmer_chi_sq)
m.migrad()
chisq = m.fval
dof = len(p) - len(iminuit.describe(redmer_chi_sq))
print('30 GeV - 100 GeV')
print('p value = {}'.format(1 - stats.chi2.cdf(chisq, dof)))

# --- 30-60 GeV ---

# Set globals before minimization
p = indep[:3]
t = theory[:3]
d = dep[:3]
err_s = stat_err[:3]
err_c = corr_err[:3]
err_a = acorr_err[:3]

# Minimize and print p
m = iminuit.Minuit(redmer_chi_sq)
m.migrad()
chisq = m.fval
dof = len(p) - len(iminuit.describe(redmer_chi_sq))
print('30 GeV - 60 GeV')
print('p value = {}'.format(1 - stats.chi2.cdf(chisq, dof)))

# --- 60-100 GeV ---

# Set globals before minimization
p = indep[3:]
t = theory[3:]
d = dep[3:]
err_s = stat_err[3:]
err_c = corr_err[3:]
err_a = acorr_err[3:]

# Minimize and print p
m = iminuit.Minuit(redmer_chi_sq)
m.migrad()
chisq = m.fval
dof = len(p) - len(iminuit.describe(redmer_chi_sq))
print('60 GeV - 100 GeV')
print('p value = {}'.format(1 - stats.chi2.cdf(chisq, dof)))
