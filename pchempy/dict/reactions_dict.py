#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Dictionary of gas IDs and their parameters
for pchempy
"""

import numpy as np

reaction_info = {
    
    '1': {
            "name": "O + O2 + CO2 -> O3 + CO2",
            "thirdbody": "dens"
        },
    
    '2': {
            "name": "O + O + CO2 -> O2 + CO2",
            "thirdbody": "dens"
        },
    
    '3': {
            "name": "O + O3 -> O2 + O2",
            "thirdbody": ""
        },
    
    '4': {
            "name": "O(1D) + CO2 -> O + CO2",
            "thirdbody": "co2"
        },
    
    '5': {
            "name": "O(1D) + H2O -> OH + OH",
            "thirdbody": ""
        },
    
    '6': {
            "name": "O(1D) + H2 -> OH + H",
            "thirdbody": ""
        },
    
    '7': {
            "name": "O(1D) + O2 -> O + O2",
            "thirdbody": "o2"
        },
    
    '8': {
            "name": "O(1D) + O3  -> O2 + O2",
            "thirdbody": ""
        },
    
    '9': {
            "name": "O(1D) + O3  -> O2 + O + O",
            "thirdbody": ""
        },
    
    '10': {
            "name": "O + HO2 -> OH + O2",
            "thirdbody": ""
        },
    
    '11': {
            "name": "O + OH -> O2 + H",
            "thirdbody": ""
        },
    
    '12': {
            "name": "H + O3 -> OH + O2",
            "thirdbody": ""
        },
    
    '13': {
            "name": "H + HO2 -> OH + OH",
            "thirdbody": ""
        },
    
    '14': {
            "name": "H + HO2 -> H2 + O2",
            "thirdbody": ""
        },
    
    '15': {
            "name": "H + HO2 -> H2O + O",
            "thirdbody": ""
        },
    
    '16': {
            "name": "OH + HO2 -> H2O + O2",
            "thirdbody": ""
        },
    
    '17': {
            "name": "HO2 + HO2 -> H2O2 + O2",
            "thirdbody": ""
        },
    
    '18': {
            "name": "OH + H2O2 -> H2O + HO2",
            "thirdbody": ""
        },
    
    '19': {
            "name": "OH + H2 -> H2O + H",
            "thirdbody": ""
        },
    
    '20': {
            "name": "H + O2 + CO2 -> HO2 + CO2",
            "thirdbody": "dens"
        },
    
    '21': {
            "name": "O + H2O2 -> OH + HO2",
            "thirdbody": ""
        },
    
    '22': {
            "name": "OH + OH -> H2O + O",
            "thirdbody": ""
        },
    
    '23': {
            "name": "OH + O3 -> HO2 + O2",
            "thirdbody": ""
        },
    
    '24': {
            "name": "HO2 + O3 -> OH + O2 + O2",
            "thirdbody": ""
        },
    
    '25': {
            "name": "HO2 + HO2 + CO2 -> H2O2 + O2 + CO2",
            "thirdbody": "dens"
        },
    
    '26': {
            "name": "OH + OH + CO2 -> H2O2 + CO2",
            "thirdbody": "dens"
        },
    
    '27': {
            "name": "H + H + CO2 -> H2 + CO2",
            "thirdbody": "dens"
        },
    
    '28': {
            "name": "NO2 + O -> NO + O2",
            "thirdbody": ""
        },
    
    '29': {
            "name": "NO + O3 -> NO2 + O2",
            "thirdbody": ""
        },
    
    '30': {
            "name": "NO + HO2 -> NO2 + OH",
            "thirdbody": ""
        },
    
    '31': {
            "name": "N + NO -> N2 + O",
            "thirdbody": ""
        },
    
    '32': {
            "name": "N + O2 -> NO + O",
            "thirdbody": ""
        },
    
    '33': {
            "name": "NO2 + H -> NO + OH",
            "thirdbody": ""
        },
    
    '34': {
            "name": "N + O -> NO",
            "thirdbody": ""
        },
    
    '35': {
            "name": "N + HO2 -> NO + OH",
            "thirdbody": ""
        },
    
    '36': {
            "name": "N + OH -> NO + H",
            "thirdbody": ""
        },
    
    '37': {
            "name": "N(2D) + O  -> N + O",
            "thirdbody": "o"
        },
    
    '38': {
            "name": "N(2D) + N2  -> N + N2",
            "thirdbody": "n2"
        },
    
    '39': {
            "name": "N(2D) + CO2  -> NO + CO",
            "thirdbody": ""
        },
    
    '40': {
            "name": "OH + CO -> CO2 + H",
            "thirdbody": ""
        },
    
    '41': {
            "name": "O + CO + M -> CO2 + M",
            "thirdbody": "dens"
        },
    
    '42': {
            "name": "",
            "thirdbody": ""
        },
    
    '43': {
            "name": "",
            "thirdbody": ""
        },
    
    '44': {
            "name": "",
            "thirdbody": ""
        },
    
    '45': {
            "name": "",
            "thirdbody": ""
        },
    
    '46': {
            "name": "",
            "thirdbody": ""
        },
    
    '47': {
            "name": "",
            "thirdbody": ""
        },
    
    '48': {
            "name": "",
            "thirdbody": ""
        },
    
    '49': {
            "name": "",
            "thirdbody": ""
        },
    
    '50': {
            "name": "",
            "thirdbody": ""
        },
    
    '51': {
            "name": "",
            "thirdbody": ""
        },
    
    '52': {
            "name": "",
            "thirdbody": ""
        },
    
    '53': {
            "name": "",
            "thirdbody": ""
        },
    
    '54': {
            "name": "",
            "thirdbody": ""
        },
    
    '55': {
            "name": "",
            "thirdbody": ""
        },  
}

reaction_c13_info = {
    '39': {
            "name": "N(2D) + (13C)O2  -> NO + CO",
            "thirdbody": ""
        },
    
    '40': {
            "name": "OH + (13C)O -> (13C)O2 + H",
            "thirdbody": ""
        },
    
    '41': {
            "name": "O + (13C)O + M -> (13C)O2 + M",
            "thirdbody": "dens"
        },
}