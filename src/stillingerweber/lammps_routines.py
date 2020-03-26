from lammps import PyLammps
import numpy as np
from scipy.optimize import leastsq
from tqdm import tqdm
import io
import sys
import os

def initialize():
    text_trap = io.StringIO()
    sys.stdout = text_trap
    L = PyLammps()
    sys.stdout = sys.__stdout__
    L.command("echo log")
    L.command("units metal")
    L.command("atom_style atomic")
    L.command("boundary p p p")
    return L

def set_lattice(L, structure, alat, reps=(3,3,3)):
    L.lattice("%s %f orient x 1 0 0 orient y 0 1 0 orient z 0 0 1"%(structure, alat))
    L.region("box block 0 %d 0 %d 0 %d"%(reps[0], reps[1], reps[2]))
    L.create_box("1 box")
    L.create_atoms("1 box")
    return L

def set_potential(L, potfile, mass=28.08):
    L.pair_style("sw")
    L.pair_coeff("* * %s Si"%potfile)
    L.neighbor("1.0 bin")
    L.neigh_modify("every 1 delay 1 check yes")
    L.mass("* %f"%mass)
    return L

def routine_energy(sw, structure, alat, reps=(3,3,3), mass=28.08, return_volume = False):
    potfile = os.path.join(os.getcwd(), "si.sw")
    sw.write(potfile)
    L = initialize()
    L = set_lattice(L, structure, alat, reps=reps)
    L = set_potential(L, potfile, mass=mass)
    L.run(1)
    pe = L.eval('pe')/L.system.natoms
    if return_volume:
        vol = L.eval('vol')/L.system.natoms
        return pe, vol
    else:
        return pe

def routine_lattice_constant(sw, structure, alat, reps=(3,3,3), mass=28.08):
    potfile = os.path.join(os.getcwd(), "si.sw")
    sw.write(potfile)
    L = initialize()
    L = set_lattice(L, structure, alat, reps=reps)
    L = set_potential(L, potfile, mass=mass)
    L.fix("1 all box/relax iso 0. vmax 0.0001 nreset 1")
    L.minimize(" 1.0e-8 1.0e-8 100000000 100000000")
    vol = L.eval('vol')/(reps[0]*reps[1]*reps[2])
    lat = vol**(1/3)
    return lat

def birch_murnaghan(vol, ene):
    a, b, c = np.polyfit(vol, ene, 2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    Bp = 4.0
    eta = (vol/V0)**(1.0/3.0)
    E = E0 + 9.0*B0*V0/16.0 * (eta**2-1.0)**2 * (6.0 + Bp*(eta**2-1.0) - 4.0*eta**2)
    return E

def routine_ev_curve(sw, structure, alat, reps=(3,3,3), mass=28.08, variation = 0.01, points=100):
    alats = np.linspace((1-variation)*alat, (1+variation)*alat, points)
    energies = []
    volumes = []
    for c, val in tqdm(enumerate(alats)):
        e, v = routine_energy(sw, structure, alat, reps=reps, mass=mass, return_volume= True)
        energies.append(e)
        volumes.append(v)
    #now fit
    efit = birch_murnaghan(volumes, energies)
    return volumes, energies, efit
