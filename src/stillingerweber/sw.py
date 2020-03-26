"""
SW class
"""

import numpy as np

class Sw:
    """
    A class to depict SW potential
    """
    def __init__(self, tag, **kwargs):
        self.tag = tag
        self.elei = 'Si'
        self.elej = 'Si'
        self.elek = 'Si'
        self.epsilon = 2.1683
        self.sigma = 2.0951
        self.a = 1.80
        self.lmbda = 21.0
        self.gamma = 1.20
        self.costheta0 = -0.333333333333
        self.A =  7.049556277
        self.B = 0.6022245584
        self.p = 4
        self.q = 0
        self.tol = 0

        allargs = ['epsilon', 'sigma', 'a', 'lambda', 'gamma', 'A', 'B', 'p', 'q', 'costheta0']
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in allargs)

    def write(self, outfile):
        """
        Write potential to a file
        """
        with open(outfile, 'w') as fout:
            fout.write("# %s\n"%str(self.tag))
            fout.write("# format of a single entry (one or more lines):\n")
            fout.write("#   element 1, element 2, element 3,\n")
            fout.write("#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol\n")
            fout.write("%s %s %s %f %f %f %f %f %f %f %f %f %f %f\n"%(self.elei, self.elej, self.elek, self.epsilon,
                                                                    self.sigma, self.a, self.lmbda, self.gamma, self.costheta0, self.A,
                                                                    self.B, self.p, self.q, self.tol))

    def phi2(self, rij, separate=False):
        """
        Calculate the two body term when rij is given
        If separate, return two dfferent functions
        rij should be numpy array
        """
        prefactor = self.A*self.epsilon
        term1 = self.B*(self.sigma**self.p)*(rij**(-1*self.p)) - (self.sigma**self.q)*(rij**(-1*self.q))
        term2 = np.exp(self.sigma/(rij-self.a*self.sigma))
        if separate:
            return prefactor*term1, term2
        else:
            return prefactor*term1*term2

    def phi3(self, rij, rik, costheta):
        """
        Calculate the three-body term for a fixed costheta
        """
        prefactor = self.lmbda*self.epsilon
        term1 = (costheta-self.costheta0)**2
        term2 = np.exp(self.gamma*self.sigma/(rij-self.a*self.sigma))
        term3 = np.exp(self.gamma*self.sigma/(rik-self.a*self.sigma))
        return prefactor*term1*term2*term3

    def energy(self, sys):
        """
        Calculate two-three body energy contribution
        """
        atoms = sys.atoms
        phi2sum = 0
        phi3sum = 0
        for i in range(1):
            atomi = atoms[i]
            for j in range(i+1, len(atoms)):
                atomj = atoms[j]
                rij, vec1 = sys.get_distance(atomi, atomj)
                phi2sum += 0.5*self.phi2(rij)
                for k in range(j+1, len(atoms)):
                    atomk = atoms[k]
                    rik, vec2 = sys.get_distance(atomi, atomk)
                    #find costheta
                    costheta = np.dot(np.array(vec1), np.array(vec2))/(rij*rik)
                    phi3sum += self.phi3(rij, rik, costheta)
        return phi2sum, phi3sum, phi2sum+phi3sum
