import re
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from scipy.integrate import simps
import xml.etree.ElementTree as ET

eV2Ha = 1 / 27.21138397
n_electrons = 383
HOMO_index = int(math.floor(n_electrons/2))
LUMO_index = HOMO_index + 1

_VALS = {'max_values': [],
         'min_values': [],
         'HOMO_values': [],
         'LUMO_values': [],
         'weight_list': [],
         'eigenvalue_master': [],
         'E_min': None,
         'E_max': None}

def get_fermi_energy(filepath):
    # write function that takes in pw.out file and returns Fermi energy
    # of most converged slab
    return 1.6017

def read_surface_file_xml(filepath):

    tree = ET.parse(filepath)
    root = tree.getroot()

    ks_energy_roots = []
    for descendant in root.iter():
        if descendant.tag == 'ks_energies':
            ks_energy_roots.append(descendant)

    for energy_root in ks_energy_roots:
        for child in energy_root:
            # find k-point weight
            if child.tag == 'k_point':
                _VALS['weight_list'].append(float(child.attrib['weight']))
            # extract all eigenvalues for each k-point
            if child.tag == 'eigenvalues':
                eigen_values = child.text.split()
                eigen_values = [float(i) for i in eigen_values]
                _VALS['min_values'].append(min(eigen_values))
                _VALS['max_values'].append(max(eigen_values))
                _VALS['HOMO_values'].append(eigen_values[HOMO_index])
                _VALS['LUMO_values'].append(eigen_values[LUMO_index])
                _VALS['eigenvalue_master'].append(eigen_values)

    _VALS['E_min'] = min(_VALS['min_values'])/eV2Ha
    _VALS['E_max'] = max(_VALS['max_values'])/eV2Ha
        
    return

def calculate_DOS(eigenvalue_list, E_max, E_min):

    qe_lambda_eV = [i*eV2Ha for i in eigenvalue_list]
    # qe_lambda_eV = eigenvalue_list*eV2Ha
    mu = 0.05 #in eV
    N = 20*len(qe_lambda_eV)

    DOS_eV, E_eV = eig2DOS_bulk(qe_lambda_eV, E_min, E_max, N, mu) # swapped emin and emax

    return DOS_eV, E_eV

def eig2DOS_bulk(lambda_val, E_min, E_max, num, sigma): #send the first eigen value 
    
    #we are sampling between minimum and maximum values of lambda
    E = np.linspace(E_min,E_max,num)
    x = colminusrow(E,lambda_val)

    DOS = np.sum(gauss_distribution(x,sigma), axis=1)
    
    return DOS, E

#function to generate gaussian distribution using output from colminusrow and sigma value
def gauss_distribution(x, s):
    p1 = -0.5*((x)/s)**2
    p2 = (s * np.sqrt(2*np.pi))
    f = np.exp(p1)/p2
    return f

def colminusrow(x, y):
    """
    1. get a 2D array- x is E (sampling points), y is list of eigen_values 
    2. xx has no. of columns as no. of eigen values, no. of rows for each element in sampling space
    3. in yy, each eigen value has a column of its own, no. of rows are no. of elements in sampling space
    4. for example there are 120 eigen values and 600 elements in sampling space, xx will have dimension of (600X120)
       and yy will have dimension of (600X120)
    """
    xx, yy = np.meshgrid(x, y, indexing = 'ij')
    xmy = xx - yy
    return xmy

if __name__ == '__main__':

    base = '/Users/jamesmccord/Dropbox (GaTech)/pBlock/fall2021/final_data/N/Base'
    file = 'data-file-schema.xml'
    file_path = os.path.join(base,file)

    E_fermi = get_fermi_energy(file_path)
    read_surface_file_xml(file_path)

    DOS_list = []
    for i in range(len(_VALS['eigenvalue_master'])):
        DOS_eV, E_eV = calculate_DOS(_VALS['eigenvalue_master'][i], _VALS['E_max']+5, _VALS['E_min']-5)
        weight = _VALS['weight_list'][i]
        DOS_list.append(DOS_eV*weight)

    DOS_array = np.array(DOS_list)
    DOS = np.sum(DOS_array, axis=0)

    plt.plot(E_eV, DOS)
    plt.show()

#
# DOS_array_100 = np.array(DOS_weighted_list_100)
# DOS_100 = np.sum(DOS_array_100, axis = 0)
# HOMO_value_100 = np.max(HOMO_values_100)*27.2114
# LUMO_value_100 = np.min(LUMO_values_100)*27.2114
#
# plt.figure(figsize=(20,10))
# #plt.xlim([0,20])
# plt.ylim([0,150])
# plt.plot(E_eV_100, DOS_100)
# plt.axvline(x=E_fermi_100, color='cyan', label='Fermi level')
# plt.axvline(x=HOMO_value_100, color='black', label='HOMO')
# plt.axvline(x=LUMO_value_100, color='green', label='LUMO')
# plt.xlabel('Energy (eV)')
# plt.ylabel('Density of States')
# plt.legend()

