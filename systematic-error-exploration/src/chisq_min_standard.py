from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from hepdata_manager import HEPData_Manager
import numpy as np
import iminuit
from scipy import stats


alice_event = 'ins1394678'
data_manager = HEPData_Manager()
if not data_manager.has_data(alice_event):
    data_manager.get_data(alice_event)
    print('Downloaded and saved {} data!'.format(alice_event))

table1_data = data_manager.load_table(alice_event, 'Table1')

global indep, dep, stat_err, corr_err, acorr_err, theory

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

print(indep)
print(corr_err)
print(acorr_err)

def chi_squ_one_err(a):
    chisq = 0

    for i in range(len(indep)):
        chisq += (dep[i] - theory[i] - np.sqrt(0.75) * a * corr_err[i])**2 \
               / stat_err[i]**2

    chisq += a**2

    return chisq

def chi_squ_two_err(a, b):
    chisq = 0

    for i in range(len(indep)):
        chisq += (dep[i] - theory[i] - a * corr_err[i] - b * acorr_err[i])**2 \
               / stat_err[i]**2

    chisq += a**2
    chisq += b**2

    return chisq

def redmer_chi_sq(a, b):
    chisq = 0

    for i in range(len(indep)):
        chisq += (dep[i] + a * corr_err[i] + b)**2 / stat_err[i]**2

    chisq += a**2

    subsum = 0

    for i in range(len(indep)):
        subsum += b**2 / acorr_err[i]**2
    subsum = subsum / len(indep)

    chisq += subsum

    return chisq

global sigma_new
sigma_new = 0
for i in range(len(acorr_err)):
    sigma_new += 1 / acorr_err[i]**2
sigma_new = sigma_new / len(acorr_err)
sigma_new = np.sqrt(1 / sigma_new)

def redmer_chi_sq_alt(a, b):
    chisq = 0

    for i in range(len(indep)):
        chisq += (dep[i] + a * corr_err[i] + b * sigma_new)**2 / stat_err[i]**2

    chisq += a**2
    chisq += b**2

    return chisq

print(iminuit.describe(redmer_chi_sq))

m = iminuit.Minuit(redmer_chi_sq)
m.migrad()
# print('args', m.args)
chisq = m.fval
print('chisq = {}'.format(chisq))
dof = len(indep) - len(iminuit.describe(redmer_chi_sq))
print('chisq/dof (dof = # of data - free params = {}) = {}'.format(dof, 
                                                                   chisq / dof))
print('p value = {}'.format(1 - stats.chi2.cdf(chisq, dof)))
