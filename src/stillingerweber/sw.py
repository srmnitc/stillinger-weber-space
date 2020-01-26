"""
SW class
"""
class Sw:
    """
    A class to depict SW potential
    """
    def __init__(self, tag, params = None):
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
        
        if params is not None:
            #params should be a dict with keywords as the variable params - epsilon, sigma, a , lambda
            self.epsilon = params['epsilon']
            self.sigma = params['sigma']
            self.a = params['a']
            self.lmbda = params['lambda']
    
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