from lammps import PyLammps
import numpy as np
from scipy.optimize import leastsq
from tqdm import tqdm
import io
import sys

class LammpsSw(PyLammps):
    def __init__(self):
        text_trap = io.StringIO()
        sys.stdout = text_trap
        super().__init__()
        sys.stdout = sys.__stdout__
        self.command("echo log")
        self.command("units metal")
        self.command("atom_style atomic")
        self.command("boundary p p p")

    def set_lattice(self, structure, alat, reps=(3,3,3)):
        self.lattice("%s %f orient x 1 0 0 orient y 0 1 0 orient z 0 0 1"%(structure, alat))
        self.region("box block 0 %d 0 %d 0 %d"%(reps[0], reps[1], reps[2]))
        self.create_box("1 box")
        self.create_atoms("1 box")

    def set_potential(self, potfile, mass=28.08):
        self.pair_style("sw")
        self.pair_coeff("* * %s Si"%potfile)
        self.neighbor("1.0 bin")
        self.neigh_modify("every 1 delay 1 check yes")
        self.mass("* %f"%mass)

    def routine_energy(self, structure, alat, potfile, mass=28.08, return_volume = False):
        self.set_lattice(structure, alat, reps=reps)
        self.set_potential(potfile, mass=mass)
        self.run(1)
        pe = self.eval('pe')/self.system.natoms
        if return_volume:
            vol = self.eval('vol')/self.system.natoms
            return pe, vol
        else:
            return pe

    def routine_lattice_constant(self, structure, alat, potfile, mass=28.08):
        self.set_lattice(structure, alat, reps=reps)
        self.set_potential(potfile, mass=mass)
        self.command("fix 1 all box/relax iso 0. vmax 0.0001 nreset 1")
        self.minimize(" 1.0e-8 1.0e-8 100000000 100000000")
        vol = self.eval('vol')/(reps[0]*reps[1]*reps[2])
        lat = vol**(1/3)
        return lat

    def birch_murnaghan(self, vol, ene):
        a, b, c = np.polyfit(vol, ene, 2)
        V0 = -b/(2*a)
        E0 = a*V0**2 + b*V0 + c
        B0 = 2*a*V0
        Bp = 4.0
        eta = (vol/V0)**(1.0/3.0)
        E = E0 + 9.0*B0*V0/16.0 * (eta**2-1.0)**2 * (6.0 + Bp*(eta**2-1.0) - 4.0*eta**2)
        return E

"""
    def routine_ev_curve(self, structure, alat, potfile, mass=28.08, variation = 0.01, points=100):
        alats = np.linspace((1-variation)*alat, (1+variation)*alat, points)
        energies = []
        volumes = []
        for c, val in tqdm(enumerate(alats)):
            e, v = routine_energy(structure, val, potfile, mass=mass, return_volume= True)
            energies.append(e)
            volumes.append(v)
        #now fit
        efit = birch_murnaghan(volumes, energies)
        return volumes, energies, efit
"""
