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

stat_cov_mat = np.diag([x**2 for x in stat_err])

corr_cov_mat = np.array([[0.75 * x * y for x in corr_err] for y in corr_err])

acorr_cov_mat = np.array([[0.5 * x * y for x in acorr_err] for y in acorr_err])

cov_mat = np.add(stat_cov_mat, corr_cov_mat)

cov_mat = np.add(cov_mat, acorr_cov_mat)

del_vec = np.array([dep[i] - theory[i] for i in range(len(dep))])

chisq = np.dot(del_vec, np.dot(np.linalg.inv(cov_mat), del_vec))

print('chisq = {}'.format(chisq))

dof = len(indep) - 2

print('chisq/dof (dof = # of data - free params = {}) = {}'.format(dof, 
                                                                   chisq/dof))
print('p value = {}'.format(1 - stats.chi2.cdf(chisq, dof)))
