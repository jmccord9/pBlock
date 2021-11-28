#plotting DOS (python script)
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import xml.etree.ElementTree as ET

E_fermi = 0.8016 #in eV (look for this value in pw.out file)

filepath = "br_base+N2.xml" #specify the filepath to xml file

eV2Ha = 1 / 27.21138397

def read_surface_file_xml(filepath):
    """function to read xml file with eigen values and kpoint weights"""
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    max_values = []
    min_values = []
    HOMO_values = []
    LUMO_values = []
    weight_list = []
    
    for i in range(12,22): #this range will vary 
        values = list(root[3][9][i][0].attrib.values())
        weight_list.append(float(values[0]))

    eigen_values_master = []
    for i in range(12,22): #this for loop is reading the xml tree (there are 10 eigen values in my case, therefore, it goes from 12 to 22)
        eigen_values = []
        for j in range(106):#this number refers to number of eigen values (<eigenvalues size="106">) <= you will find this in the .xml file
            eigen_values.append(float(root[3][9][i][2].text.split()[j]))
        eigen_values_master.append(eigen_values)
    
    for i in range(len(eigen_values_master)):
        min_values.append(min(eigen_values_master[i]))#for sampling 
        max_values.append(max(eigen_values_master[i]))#for sampling
        HOMO_values.append(eigen_values_master[i][95])#homo for each k-point (eigen value for n_electrons/2-th state)
        LUMO_values.append(eigen_values_master[i][96])#lumo for each k-point (eigen value for n_electrons/2 + 1 -th state)

    E_min = min(min_values)/eV2Ha
    E_max = max(max_values)/eV2Ha
        
    return weight_list, eigen_values_master, HOMO_values, LUMO_values, E_min, E_max

#calculating the DOS of TiO2 bulk for each list of eigen values

def calculate_DOS(eigen_values_list, E_max, E_min):

    qe_lambda_eV = []
    
    #loop to convert eigen values to eV
    for i in range(len(eigen_values_list)):
        qe_lambda_eV.append(eigen_values_list[i]/eV2Ha)
        
    mu = 0.05 #in eV
    N = 20*len(qe_lambda_eV)

    #fixed E_max and E_min for sampling
    DOS_eV, E_eV = eig2DOS_bulk(qe_lambda_eV, E_max, E_min, N, mu)
    return DOS_eV, E_eV

#converting from eigen values to DOS for TiO2 bulk

def eig2DOS_bulk(lambda_val, E_min, E_max, num, sigma): #send the first eigen value 
    
    #we are sampling between minimum and maximum values of lambda
    
    E = np.linspace(E_min,E_max,num)
    DOS = np.sum(gauss_distribution(colminusrow(E, lambda_val), sigma), axis=1)
    
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

DOS_list_100 = []
for i in range(10): #10 kpoints
    DOS_eV_100, E_eV_100 = calculate_DOS(eigen_values_master_100[i], E_max_100+5, E_min_100-5)
    #print(DOS_eV)
    #print(np.count_nonzero(DOS_eV))
    DOS_list_100.append(DOS_eV_100)

DOS_weighted_list_100 = []
#print(len(DOS_list))   
for i in range(len(DOS_list_100)):
    DOS_weighted_list_100.append((DOS_list_100[i]*weight_list_100[i]))# sum of weights is equal to 2

DOS_array_100 = np.array(DOS_weighted_list_100)
DOS_100 = np.sum(DOS_array_100, axis = 0)
HOMO_value_100 = np.max(HOMO_values_100)*27.2114
LUMO_value_100 = np.min(LUMO_values_100)*27.2114

plt.figure(figsize=(20,10))
#plt.xlim([0,20])
plt.ylim([0,150])
plt.plot(E_eV_100, DOS_100)
plt.axvline(x=E_fermi_100, color='cyan', label='Fermi level')
plt.axvline(x=HOMO_value_100, color='black', label='HOMO')
plt.axvline(x=LUMO_value_100, color='green', label='LUMO')
plt.xlabel('Energy (eV)')
plt.ylabel('Density of States')
plt.legend()

