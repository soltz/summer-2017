from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from hepdata_manager import HEPData_Manager
import numpy as np
import iminuit


alice_event = 'ins1394678'
data_manager = HEPData_Manager()
if not data_manager.has_data(alice_event):
    data_manager.get_data(alice_event)
    print('Downloaded and saved {} data!'.format(alice_event))

table1_data = data_manager.load_table(alice_event, 'Table1')

global indep, dep, stat_err, corr_err, theory

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

def chi_squ_one_err(a):
    chisq = 0

    for i in range(len(indep)):
        chisq += (dep[i] - theory[i] - np.sqrt(0.75) * a * corr_err[i])**2 / stat_err[i]**2

    chisq += a**2

    return chisq

print(iminuit.describe(chi_squ_one_err))

m = iminuit.Minuit(chi_squ_one_err)
m.migrad()
# print('args', m.args)
print('chisq = {}'.format(m.fval))
dof = len(indep) - len(iminuit.describe(chi_squ_one_err))
print('chisq/dof (dof = # of data - free params = {}) = {}'.format(dof, m.fval/dof))
