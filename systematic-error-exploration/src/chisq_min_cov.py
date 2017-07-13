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

def calc_p_val(dep, theory, stat_err, corr_err, acorr_err):
    stat_cov_mat = np.diag([x**2 for x in stat_err])
    corr_cov_mat = np.array([[x * y for x in corr_err] for y in corr_err])
    err_n = 0
    for i in range(len(acorr_err)):
        err_n += 1 / acorr_err[i]**2
    err_n = err_n / len(acorr_err)
    err_n = np.sqrt(1 / err_n)
    acorr_cov_mat = np.array([[err_n * err_n for x in acorr_err] for y in acorr_err])
    cov_mat = np.add(stat_cov_mat, corr_cov_mat)
    cov_mat = np.add(cov_mat, acorr_cov_mat)
    del_vec = np.array([dep[i] - theory[i] for i in range(len(dep))])

    chisq = np.dot(del_vec, np.dot(np.linalg.inv(cov_mat), del_vec))

    dof = len(dep) - 2

    return 1 - stats.chi2.cdf(chisq, dof)

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

# Set variables before minimization
p = indep
t = theory
d = dep
err_s = stat_err
err_c = corr_err
err_a = acorr_err

# Minimize and print p
p_val = calc_p_val(d, t, err_s, err_c, err_a)
print('30 GeV - 100 GeV')
print('p value = {}'.format(p_val))

# --- 30-60 GeV ---

# Set variables before minimization
p = indep[:3]
t = theory[:3]
d = dep[:3]
err_s = stat_err[:3]
err_c = corr_err[:3]
err_a = acorr_err[:3]

# Minimize and print p
p_val = calc_p_val(d, t, err_s, err_c, err_a)
print('30 GeV - 60 GeV')
print('p value = {}'.format(p_val))

# --- 60-100 GeV ---

# Set variables before minimization
p = indep[3:]
t = theory[3:]
d = dep[3:]
err_s = stat_err[3:]
err_c = corr_err[3:]
err_a = acorr_err[3:]

# Minimize and print p
p_val = calc_p_val(d, t, err_s, err_c, err_a)
print('60 GeV - 100 GeV')
print('p value = {}'.format(p_val))
