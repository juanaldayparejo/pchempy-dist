#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Dictionary of physics constants
"""

import numpy as np

phys_const = {
    'k_B': 1.38065e-23,              # J K-1 Boltzmann constant
    'sig_B': 5.67037e-8,             # W m-2 K-4 Stephan Boltzmann constant
    'R': 8.31446,                    # J mol-1 K-1 universal gas constant
    'G': 6.67199976E-11,             # m3 kg-1 s-2 universal gravitational constant
    'eps_LJ': 59.7*5.67037e-8,       # J depth of the Lennard-Jones potential well for H2
    'c_p': 14300,                    # J K-1 hydrogen specific heat
    'h': 6.62607,                    # Js Planck's constant
    'hbar': 1.05457,                 # Js
    'N_A': 6.02214e23,               # Avagadro's number (mol-1)
}
