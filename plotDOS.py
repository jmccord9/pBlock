import numpy as np
import os
import math
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from parsimonious.grammar import Grammar

eV2Ha = 1/27.21138397

_VALS = {'max_values': [],
         'min_values': [],
         'HOMO_values': [],
         'LUMO_values': [],
         'weight_list': [],
         'eigenvalue_master': [],
         'E_min': None,
         'E_max': None,
         'E_fermi': None,
         'HOMO_index': None,
         'LUMO_index': None,
         'n_electrons': None}

def parse_pw_out(filepath):
    # Parse pw.out for Fermi energy and number of electrons
    with open(filepath) as f:
        for line in f:
            li = line.strip()
            if li.startswith('number of electrons'):
                _VALS['n_electrons'] = int(float(li.split()[-1]))
            if li.startswith('the Fermi energy is'):
                fermi_string = li

    _VALS['E_fermi'] = float(fermi_string.split()[-2])
    _VALS['HOMO_index'] = int(math.ceil(_VALS['n_electrons']/2)) -1
    _VALS['LUMO_index'] = _VALS['HOMO_index'] + 1
    print("Hello")
    return

def read_surface_file_xml(filepath):
    # Parse QuantumEspresso output xml file to get k-points information.
    # Assumes standard QuantumEspresso xml tree structure and names.
    tree = ET.parse(filepath)
    root = tree.getroot()
    # Search through all of child nodes to find k-point information.
    ks_energy_roots = []
    for descendant in root.iter():
        if descendant.tag == 'ks_energies':
            ks_energy_roots.append(descendant)
    # Iterate through every k-point state to extract relevant information.
    for energy_root in ks_energy_roots:
        for child in energy_root:
            # find k-point weights
            if child.tag == 'k_point':
                _VALS['weight_list'].append(float(child.attrib['weight']))
            # extract all eigenvalues for each k-point
            if child.tag == 'eigenvalues':
                eigen_values = np.array(child.text.split()).astype(float)
                # eigen_values = [float(i) for i in child.text.split()]
                _VALS['min_values'].append(np.amin(eigen_values,axis=0))
                _VALS['max_values'].append(np.amax(eigen_values,axis=0))
                _VALS['HOMO_values'].append(eigen_values[_VALS['HOMO_index']])
                _VALS['LUMO_values'].append(eigen_values[_VALS['LUMO_index']])
                _VALS['eigenvalue_master'].append(eigen_values)

    _VALS['E_min'] = min(_VALS['min_values'])/eV2Ha
    _VALS['E_max'] = max(_VALS['max_values'])/eV2Ha
    return

def calculate_DOS(eigenvalue_list, E_max, E_min):
    # Calculate the DOS using the list of eigenvalues (eigenvalue_master).
    qe_lambda_eV = eigenvalue_list/eV2Ha # I don't think this needs to be multiplied...
    # qe_lambda_eV = [i*eV2Ha for i in eigenvalue_list]
    mu = 0.05 #in eV
    # N = 20*len(qe_lambda_eV)
    N = 20*np.size(qe_lambda_eV,0)
    DOS_eV, E_eV = eig2DOS_bulk(qe_lambda_eV, E_min, E_max, N, mu)
    return DOS_eV, E_eV

def eig2DOS_bulk(lambda_val, E_min, E_max, num, sigma):
    # Sample between minimum and maximum values of lambda.
    E = np.linspace(E_min,E_max,num)
    x = colminusrow(E,lambda_val)
    DOS = np.sum(gauss_distribution(x,sigma), axis=1)
    return DOS, E

def gauss_distribution(x, s):
    # Generate Gaussian distribution from variance and energy distribution values.
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
    xml_file, out_file = 'data-file-schema.xml', 'pw.out'
    xml_path, out_path = os.path.join(base,xml_file), os.path.join(base,out_file)
    fig_path = os.path.join('/Users/jamesmccord/Dropbox (GaTech)/pBlock/fall2021/final_data/N/Base','DOS_more_zoom.png')

    parse_pw_out(out_path)
    read_surface_file_xml(xml_path)

    DOS_list = []
    for i in range(len(_VALS['eigenvalue_master'])):
        DOS_eV, E_eV = calculate_DOS(_VALS['eigenvalue_master'][i], _VALS['E_max']+5, _VALS['E_min']-5)
        weight = _VALS['weight_list'][i]
        DOS_list.append(DOS_eV*weight)

    DOS_array = np.array(DOS_list)
    DOS = np.sum(DOS_array, axis=0)/_VALS['n_electrons']

    plt.plot(E_eV, DOS)
    plt.axvline(x=_VALS['E_fermi'], color='cyan', label='Fermi level')
    plt.axvline(x=np.max(_VALS['HOMO_values'])/eV2Ha, color='black', label='HOMO')
    plt.axvline(x=np.min(_VALS['LUMO_values'])/eV2Ha, color='green', label='LUMO')
    plt.legend()
    plt.xlabel('Energy (eV)')
    plt.ylabel('Density of States')
    plt.xlim([-5, 5])
    # plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.show()

