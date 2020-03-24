import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

def birch_murn(params, x):
    """
    function that calculates energy using birch murnaghan
    equation of state.
    """
    E0 = params[0]
    V0 = params[1]
    B0 = params[2]
    BC = params[3]
    ratio = (V0/(x**3))**(2.0/3.0)
    term1 = BC*(ratio-1)**3
    term2 = (ratio-1)**2
    term3 = 6-4*ratio
    prefactor = (9*V0*B0/16.00)
    fullterm = E0 + prefactor*(term1 + term2*term3)
    return fullterm

def guess_params(x, e):
    """
    Guess initial params for fitting of birch murnaghan
    equation of state
    """
    a, b, c = np.polyfit(x**3, e, 2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    BC = 4.0
    params = [E0, V0, B0, BC]
    return params

def difference_fn(params, x, e):
    """
    Diff between fit and actual values
    to minimise - for birch murnaghan
    """
    return (e - birch_murn(params, x))

def fit_ev(x, e):
    """
    Fit a e-v file of the format -
    xlat , energy
    """
    #load file
    params = guess_params(x, e)
    birch_murn_params, ier = leastsq(difference_fn, params, args=(x,e))
    xnew = np.linspace(x.min(),x.max(),1000)
    enew = birch_murn(birch_murn_params, xnew)
    return xnew, enew
