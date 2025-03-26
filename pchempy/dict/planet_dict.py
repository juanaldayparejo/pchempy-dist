#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Dictionary with parameters specific for each planet
"""

import numpy as np


#Planetary constants
#################################################################################################

planet_const = {

    'Mars':{
        'Radius': 3397.2e3,              # Mars radius (m)
        'Mass': 6.4169e23,               #Mass (kg)
        'daytosec': 88775.,              # Length of a day (s)
        'rotrate': 2.0*np.pi/88775.,     # Rotation rate (rad s-1)
        'g0': 3.72,                      # Gravity at the surface (m s-2)
        'mugaz': 43.49,                  # Atmosphere molar mars (g mol-1)
        'yeartoday': 668.6,              # Martian year (days)
        'd_perihelion': 1.381436,        # Minimum Sun-Mars distance (AU)
        'd_aphelion': 1.66593,           # Maximum Sun-Mars distance (AU)
        'day_perihelion': 485.,          # Perihelion date (days since Ls=0 - N. Spring)
        'obliquity': 25.2,               # Obliquity (deg)
    }

}


#Diffusion coefficients
###################################################################################################

diffusion_coeff = {

    'Mars':{

        "39": {
            "name": "H2",
            "A": 2.23e17,
            "s": 0.75,
            "Btherm": -0.25,
                },
        "48": {
            "name": "H",
            "A": 8.4e17,
            "s": 0.597,
            "Btherm": -0.25,
                },

    },
}

#Upper boundary conditions
#######################################################################################################

#Upper boundary conditions (Type 1 - Fixed density (m-3) ; Type 2 - Fixed flux (m-2 s-1) ; Type 3 - Fixed velocity (m s-1) ; Type 4 - Effusion velocity)
#If species is not present it is assumed to be Type 2 with flux = 0.0

upper_bc = {

    'Mars':{

        "39": {
            "name": "H2",
            "type": 3,
            "value": 3.4e1/100., #m s-1
            #"value": 0.0/100., #m s-1
                },
        "45": {
            "name": "O",
            "type": 2,
            "value": 1.2e8*1.0e4, #m-2 s-1
            #"value": 0.0*1.0e4, #m-2 s-1
                },
        "48": {
            "name": "H",
            "type": 3,
            "value": 3.1e3/100., #m s-1
            #"value": 0.0*1.0e4, #m-2 s-1
                },

    }
}

#Lower boundary conditions
#######################################################################################################

#Lower boundary conditions (Type 1 - Fixed density (m-3) ; Type 2 - Fixed flux (m-2 s-1) ; Type 3 - Fixed velocity (m s-1))
#If species is not present it is assumed to be Type 2 with flux = 0.0

lower_bc = {

    'Mars':{
        "1": {
            "name": "H2O",
            "type": 2,
            "value": 0.0,   #m-2 s-1
                },
        "76": {
            "name": "Ar",
            "type": 2,
            "value": 0.0,   #density given by input profile
                },
        "22": {
            "name": "N2",
            "type": 2,
            "value": 0.0,   #density given by input profile
                },
        "2": {
            "name": "CO2",
            "type": 2,
            "value": 0.0,   #density given by input profile
                },
    }
}