from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import errno
from fudge.gnd.covariances import base
from fudge.gnd import styles
from fudge.gnd import physicalQuantity
from fudge.gnd.covariances import covarianceSuite
from fudge.gnd.covariances import section
from fudge.gnd.covariances import mixed
from xData import gridded
from xData import array
from xData import gridded
from xData import axes as axesModule
from xData import values as valuesModule
from xData import link as linkModule
from hepdata_manager import HEPData_Manager
import numpy as np


def get_GND_covariance(arr, energyBounds, name='default'):
    GNDarray = array.full(shape=arr.shape, data=arr[np.tril_indices(len(arr))],
                          symmetry='lower')

    axes = axesModule.axes(labelsUnits={ 0 : ('matrix_elements', ''), 
                                         1 : ('column_energy_bounds', 'GeV/c'),
                                         2 : ('row_energy_bounds', 'GeV/c') })
    axes[2] = axesModule.grid(axes[2].label, axes[2].index, axes[2].unit, 
                              style=axesModule.boundariesGridToken,
                              values=valuesModule.values(energyBounds))
    axes[1] = axesModule.grid(axes[1].label, axes[1].index, axes[1].unit, 
                              style=axesModule.linkGridToken,
                              values=linkModule.link(link=axes[2].values,
                                                     relative=True))

    GNDgridded = gridded.gridded2d(axes, GNDarray)

    covariance = base.covarianceMatrix(name, type="absolute",
                                       matrix=GNDgridded)
    return covariance

def generate_full_covariance(covariances):
    CS = covarianceSuite.covarianceSuite(projectile="Au197", target="Au197", evaluation="RHIC data v0.1")
    # 'evaluation' on previous line is like a version number for the covarianceSuite

    CS.styles.add(
        styles.evaluated(
            label="eval", derivedFrom="",
            temperature=physicalQuantity.temperature(1e9, 'K'),
            library="ALICE HEP data", version="0.1"))

    mixedData = mixed.mixedForm(label="eval")

    # assuming you have two covarianceMatrix instances that you want to combine:
    for i in covariances:
        mixedData.addComponent(i)
    # mixedData.addComponent( covariance2 )

    rowData = section.rowData( path="/xpath/to/data" )
    columnData = section.columnData( path="/xpath/to/other/data" )  # only necessary if not equal to rowData
    thisSection = section.section(label="Pb-Pb", rowData=rowData, columnData=columnData)
    thisSection.add(mixedData)

    CS.addSection(thisSection)

    return CS

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

alice_event = 'ins1394678'
data_manager = HEPData_Manager()
if not data_manager.has_data(alice_event):
    data_manager.get_data(alice_event)
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

stat_cov = get_GND_covariance(stat_cov_mat, indep, 'statistical')
corr_cov = get_GND_covariance(corr_cov_mat, indep, 'correlated')
# print('\n'.join(stat_cov.toXMLList()))
# print('\n'.join(corr_cov.toXMLList()))

CS = generate_full_covariance([stat_cov, corr_cov])
_create_dir('./data')
CS.saveToFile("./data/covarianceExample.xml")

