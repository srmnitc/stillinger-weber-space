import numpy as np
from stillingerweber.sw import Sw
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib import colors as mcolors
from matplotlib.patches import Polygon
import ipywidgets as widgets
import random

def create_slider(mn, mx, step, val, name):
    """
    Create a slider
    """
    w1 = widgets.FloatSlider(
        value=val,
        min=mn,
        max=mx,
        step=step,
        description=name,
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
    )
    return w1

def make_plot(rij, rik, phi2, phi3, threebody=False):
    """
    Make sw plots
    """
    try:
        fig.close()
    except:
        pass
    
    if not threebody:
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(1,1,1)
        #ax1 = fig.add_subplot(1,2,2)
        ax.plot(rij, phi2, color="#AD1457", linewidth=3)
        ax.axhline(0, c="#FB8C00")
        #ax.set_ylim(-10, 100)
        #ax.set_xlim(0, 7)
        ax.set_xlabel("Distance")
        ax.set_ylabel("Energy")
        ax.set_title("$\phi_2$")
    else:
        fig = plt.figure(figsize=(6, 3))
        ax1 = fig.add_subplot(1,1,1)
        X, Y = np.meshgrid(rij, rik)
        pic = ax1.contour(X, Y, phi3, 20, cmap='magma')
        ax1.set_xlabel("$r_{ij}$")
        ax1.set_ylabel("$r_{ik}$")
        #ax1.set_ylim(0, 8)
        #ax1.set_xlim(0, 8)
        ax1.set_title("$\phi_3$")
        fig.colorbar(pic)

def run_sw_2d(epsilon=None, sigma=None, a=None, lmbda=None, 
          gamma=None, costheta0=-0.333333, A=None, B=None, 
          p=None, q=None, costheta=None):
    
    sw = Sw("tag", epsilon=epsilon, sigma=sigma, a=a, lmbda=lmbda,
           gamma=gamma, costheta0=costheta0, A=A, B=B,
           p=p, q=q)
    
    rij = np.linspace(1.0, sigma*a, 1000)
    #rik = np.linspace(1.0, sigma*a, 1000)
    phi2 = sw.phi2(rij)
    
    #X, Y = np.meshgrid(rij, rik)
    #phi3 = sw.phi3(X, Y, costheta)
    
    make_plot(rij, rij, phi2, phi2)

def run_sw_3d(epsilon=None, sigma=None, a=None, lmbda=None, 
          gamma=None, costheta0=-0.333333, A=None, B=None, 
          p=None, q=None, costheta=None):
    
    sw = Sw("tag", epsilon=epsilon, sigma=sigma, a=a, lmbda=lmbda,
           gamma=gamma, costheta0=costheta0, A=A, B=B,
           p=p, q=q)
    
    rij = np.linspace(1.0, sigma*a, 1000)
    rik = np.linspace(1.0, sigma*a, 1000)
    phi2 = sw.phi2(rij)
    
    X, Y = np.meshgrid(rij, rik)
    phi3 = sw.phi3(X, Y, costheta)
    
    make_plot(rij, rik, phi2, phi3, threebody=True)

def create_sw_widgets(threebody=False):
    params = ['epsilon', 'sigma', 'a', 'lambda', 'gamma', 'A', 'B', 'p', 'q', 'costheta']
    valmin = [0.1,   0.1,    0.1, 0,    0.1,  0.1,  0.1,  0,  0, -1.0]
    valmax = [5.0,   5.0,    3.0, 50,   3.0,  15,   7,    12, 12, 1.0]
    valitr = [0.1,   0.1,    0.1, 1,    0.1,  0.1,  0.1,  1,  1,  0.01]
    valdef = [2.183, 2.0951, 1.8, 21.0, 1.20, 7.05, 0.60, 4,  0,  -0.3333]
    ws = []
    for c, param in enumerate(params):
        w = create_slider(valmin[c], valmax[c], valitr[c], valdef[c], param)
        ws.append(w)
    
    if not threebody:
        out = widgets.interactive_output(run_sw_2d, {'epsilon': ws[0], 'sigma': ws[1], 'a': ws[2], 
                        'lmbda': ws[3], 'gamma': ws[4], 'A': ws[5],
                        'B': ws[6], 'p': ws[7], 'q': ws[8], 'costheta':ws[9]})
    else:
        out = widgets.interactive_output(run_sw_3d, {'epsilon': ws[0], 'sigma': ws[1], 'a': ws[2], 
                        'lmbda': ws[3], 'gamma': ws[4], 'A': ws[5],
                        'B': ws[6], 'p': ws[7], 'q': ws[8], 'costheta':ws[9]})

    ui = widgets.GridBox(ws, layout=widgets.Layout(grid_template_columns="repeat(2, 500px)"))
    display(ui, out)