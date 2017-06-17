from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
from fudge.gnd.covariances import base
from xData import array, gridded
from xData import axes as axesModule, values as valuesModule, link as linkModule
from hepdata_manager import HEPData_Manager
import numpy as np


def get_GND_covariance(arr, energyBounds):
    GNDarray = array.full(shape=arr.shape, data=arr[ np.tril_indices(len(arr)) ], symmetry='lower')

    axes = axesModule.axes(labelsUnits = { 0 : ('matrix_elements', ''), 
                                                1 : ('column_energy_bounds', 'GeV/c'),
                                                2 : ('row_energy_bounds', 'GeV/c') })
    axes[2] = axesModule.grid(axes[2].label, axes[2].index, axes[2].unit, 
                style = axesModule.boundariesGridToken, values = valuesModule.values(energyBounds))
    axes[1] = axesModule.grid(axes[1].label, axes[1].index, axes[1].unit, 
            style = axesModule.linkGridToken, values = linkModule.link(link = axes[2].values, relative = True))

    GNDgridded = gridded.gridded2d(axes, GNDarray)

    covariance = base.covarianceMatrix("example", type="absolute", matrix=GNDgridded)
    return covariance

alice_event = 'ins1394678'
data_manager = HEPData_Manager()
if not data_manager.has_data(alice_event):
    data_manager.get_data(alice_event) # INVESTIGATE WHY THIS FAILS IN PYTHON2
    print('Downloaded and saved {} data!'.format(alice_event))

table1_data = data_manager.load_table(alice_event, 'Table1')

indep = [(x['high'] + x['low'])/2
         for x in table1_data['independent_variables'][0]['values']]

dep = [x['value'] for x in table1_data['dependent_variables'][0]['values']]

stat_err = [x['errors'][0]['symerror']
            for x in table1_data['dependent_variables'][0]['values']]

corr_err = [x['errors'][2]['symerror']
            for x in table1_data['dependent_variables'][0]['values']]

acorr_err = [x['errors'][1]['symerror']
             for x in table1_data['dependent_variables'][0]['values']]

stat_cov_mat = np.diag([x**2 for x in stat_err])

corr_cov_mat = np.array([[0.75 * x * y for x in corr_err] for y in corr_err])

stat_cov = get_GND_covariance(stat_cov_mat, indep)
corr_cov = get_GND_covariance(corr_cov_mat, indep)
print('\n'.join(stat_cov.toXMLList()))
print('\n'.join(corr_cov.toXMLList()))
