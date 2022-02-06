from ase.io import read
from ase.optimize import BFGSLineSearch as Opt
from espresso import iEspresso as Espresso
import numpy as np
import socket
#import fcntl
#import struct
#import psutil

sim = read('atoms.traj', index=-1) #Read in the structure built by the other script
calc = Espresso(atoms=sim,
                pw = 500,
                convergence = {"energy": 1e-6, \
                               "mixing": 0.05, \
                               "maxsteps": 750, \
                               "diag": "david"},
                smearing = "fd",
                sigma = 0.01,
                nbands = -10,
                xc = 'BEEF-vdW',
                kpts = (4,4,1),
                spinpol = True,
                dipole = {"status" : True},
                nosym=True)
#calc =  SocketIOCalculator(calc_, unixsocket=socket_)
sim.calc = calc
dyn = Opt(sim,trajectory='atoms.traj',logfile='opt.log',master=None)
dyn.run(fmax=0.03)
e = sim.get_potential_energy()
with open('energy.txt', 'w') as f:
    f.write(str(e))
sim.calc.close()
