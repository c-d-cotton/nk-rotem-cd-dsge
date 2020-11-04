#!/usr/bin/env python3
# PYTHON_PREAMBLE_START_STANDARD:{{{

# Christopher David Cotton (c)
# http://www.cdcotton.com

# modules needed for preamble
import importlib
import os
from pathlib import Path
import sys

# Get full real filename
__fullrealfile__ = os.path.abspath(__file__)

# Function to get git directory containing this file
def getprojectdir(filename):
    curlevel = filename
    while curlevel is not '/':
        curlevel = os.path.dirname(curlevel)
        if os.path.exists(curlevel + '/.git/'):
            return(curlevel + '/')
    return(None)

# Directory of project
__projectdir__ = Path(getprojectdir(__fullrealfile__))

# Function to call functions from files by their absolute path.
# Imports modules if they've not already been imported
# First argument is filename, second is function name, third is dictionary containing loaded modules.
modulesdict = {}
def importattr(modulefilename, func, modulesdict = modulesdict):
    # get modulefilename as string to prevent problems in <= python3.5 with pathlib -> os
    modulefilename = str(modulefilename)
    # if function in this file
    if modulefilename == __fullrealfile__:
        return(eval(func))
    else:
        # add file to moduledict if not there already
        if modulefilename not in modulesdict:
            # check filename exists
            if not os.path.isfile(modulefilename):
                raise Exception('Module not exists: ' + modulefilename + '. Function: ' + func + '. Filename called from: ' + __fullrealfile__ + '.')
            # add directory to path
            sys.path.append(os.path.dirname(modulefilename))
            # actually add module to moduledict
            modulesdict[modulefilename] = importlib.import_module(''.join(os.path.basename(modulefilename).split('.')[: -1]))

        # get the actual function from the file and return it
        return(getattr(modulesdict[modulefilename], func))

# PYTHON_PREAMBLE_END:}}}

import numpy as np

def getinputdict(loglineareqs = True):
    inputdict = {}

    inputdict['paramssdict'] = {'GAMMA': 1, 'BETA': 0.96 ** (1/4), 'ETA': 2, 'ALPHA': 0.3, 'RHO_A': 0.95, 'Abar': 1, 'Pistar': 1.02 ** (1/4), 'PHIpi': 1.5, 'MU': 100, 'SIGMA': 8, 'DELTA': 0.1}
    inputdict['states'] = ['A', 'K']
    inputdict['controls'] = ['C', 'R', 'W', 'L', 'Y', 'MC', 'Omega', 'I', 'Pi']
    inputdict['shocks'] = ['epsilon_I']

    # equations:{{{
    inputdict['equations'] = []

    # household
    if loglineareqs is True:
        inputdict['equations'].append('-GAMMA * C = R_p - GAMMA * C_p')
    else:
        inputdict['equations'].append('C^(-GAMMA) = BETA*R_p*C_p^(-GAMMA)')
    if loglineareqs is True:
        inputdict['equations'].append('-GAMMA * C = I - Pi_p - GAMMA * C_p')
    else:
        inputdict['equations'].append('C^(-GAMMA) = BETA*I/Pi_p*C_p^(-GAMMA)')
    if loglineareqs is True:
        inputdict['equations'].append('W - GAMMA * C = ETA * L')
    else:
        inputdict['equations'].append('W*C^(-GAMMA) = L^(ETA)')

    # firm production
    if loglineareqs is True:
        inputdict['equations'].append('MC = R_ss / (R_ss - 1 + DELTA) * R - A + (1 - ALPHA) * K - (1 - ALPHA) * L')
    else:
        inputdict['equations'].append('MC = (R - 1 + DELTA) / (ALPHA * A * K**(ALPHA-1) * L**(1-ALPHA) )')
    if loglineareqs is True:
        inputdict['equations'].append('MC = W - A - ALPHA * K + ALPHA * L')
    else:
        inputdict['equations'].append('MC = W / ((1 - ALPHA) * A * K**ALPHA * L**(-ALPHA) )')
    if loglineareqs is True:
        inputdict['equations'].append('Y = A + ALPHA * K + (1 - ALPHA) * L')
    else:
        inputdict['equations'].append('Y = A * K**ALPHA * L**(1-ALPHA)')
    if loglineareqs is True:
        inputdict['equations'].append('Omega - Y = - MC_ss / (1 - MC_ss) * MC')
    else:
        inputdict['equations'].append('Omega / Y = 1 - MC')

    # firm pricing
    if loglineareqs is True:
        inputdict['equations'].append('Pi = SIGMA * MC_ss / MU * MC + BETA * Pi_p')
    else:
        inputdict['equations'].append('MU * log(Pi) = SIGMA * MC - (SIGMA - 1) + BETA * MU * log(Pi_p)')
    
    # exogenous process
    if loglineareqs is True:
        inputdict['equations'].append('A_p = RHO_A * A')
    else:
        inputdict['equations'].append('log(A_p) = RHO_A*log(A) + (1 - RHO_A) * log(Abar)')

    # monetary policy
    if loglineareqs is True:
        inputdict['equations'].append('I = R_p + PHIpi * Pi + epsilon_I')
    else:
        inputdict['equations'].append('I = R_p * Pistar * (Pi / Pistar) ** PHIpi * exp(epsilon_I)')


    # resource
    if loglineareqs is True:
        inputdict['equations'].append('C_ss * C + K_ss * K_p = Y_ss * Y + (1 - DELTA) * K_ss * K')
    else:
        inputdict['equations'].append('C + K_p = Y + (1 - DELTA) * K')
        
    # equations:}}}

    p = inputdict['paramssdict']
    p['Pi'] = p['Pistar']
    p['MC'] = 1 / p['SIGMA'] * (p['SIGMA'] - 1 + (1 - p['BETA']) * p['MU'] * np.log(p['Pi']))

    p['A'] = p['Abar']
    p['R'] = 1/p['BETA']
    p['I'] = p['R'] * p['Pi']

    k = (p['MC'] * p['ALPHA'] * p['A'] / (p['R'] - 1 + p['DELTA'])) ** (1 / (1 - p['ALPHA']))
    p['W'] = p['A'] * p['MC'] * (1 - p['ALPHA']) * k ** p['ALPHA']
    y = p['A'] * k**p['ALPHA']
    c = y - p['DELTA'] * k
    p['L'] = (p['W'] * c**(-p['GAMMA'])) ** (1 / (p['GAMMA'] + p['ETA']))
    p['C'] = c * p['L']
    p['K'] = k * p['L']
    p['Y'] = y * p['L']

    p['Omega'] = p['Y'] * (1 - p['MC'])

    if loglineareqs is True:
        inputdict['loglineareqs'] = True
    else:
        inputdict['logvars'] = inputdict['states'] + inputdict['controls']
    inputdict['irfshocks'] = ['A', 'epsilon_I']

    # save stuff
    inputdict['savefolder'] = __projectdir__ / Path('temp/')

    # main vars
    inputdict['mainvars'] = ['C', 'R', 'Pi', 'I', 'K', 'L']

    return(inputdict)


def check():
    inputdict_loglin = getinputdict(loglineareqs = True)
    inputdict_log = getinputdict(loglineareqs = False)
    importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsgediff_func.py'), 'checksame_inputdict')(inputdict_loglin, inputdict_log)
    

def dsgefull():
    inputdict = getinputdict()
    importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsge_bkdiscrete_func.py'), 'discretelineardsgefull')(inputdict)


# Run:{{{1
check()
dsgefull()
