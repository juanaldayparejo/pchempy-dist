#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Dictionary of gas IDs and their parameters
for pchempy
"""

import numpy as np


#############################################################################
# DICTIONARY
#############################################################################

gas_id = {
    'H2O': 1,
    'CO2': 2,
    'O3': 3,
    'N2O': 4,
    'CO': 5,
    'CH4': 6,
    'O2': 7,
    'NO': 8,
    'SO2': 9,
    'NO2': 10,
    'NH3': 11,
    'HNO3': 12,
    'OH': 13,
    'HF': 14,
    'HCL': 15,
    'HBr': 16,
    'HI': 17,
    'ClO': 18,
    'OCS': 19,
    'H2CO': 20,
    'HOCl': 21,
    'N2': 22,
    'HCN': 23,
    'CH3Cl': 24,
    'H2O2': 25,
    'C2H2': 26,
    'C2H6': 27,
    'PH3': 28,
    'C2N2': 29,
    'C4H2': 30,
    'HC3N': 31,
    'C2H4': 32,
    'GeH4': 33,
    'C3H8': 34,
    'HCOOH': 35,
    'H2S': 36,
    'COF2': 37,
    'SF6': 38,
    'H2': 39,
    'He': 40,
    'AsH3': 41,
    'C3H4': 42,
    'ClONO2': 43,
    'HO2': 44,
    'O': 45,
    'NO+': 46,
    'CH3OH': 47,
    'H': 48,
    'C6H6': 49,
    'CH3CN': 50,
    'CH2NH': 51,
    'C2H3CN': 52,
    'HCP': 53,
    'CS': 54,
    'HC5N': 55,
    'HC7N': 56,
    'C2H5CN': 57,
    'CH3NH2': 58,
    'HNC': 59,
    'Na': 60,
    'K': 61,
    'TiO': 62,
    'VO': 63,
    'CH2CCH2': 64,
    'C4N2': 65,
    'C5H5N': 66,
    'C5H4N2': 67,
    'C7H8': 68,
    'C8H6': 69,
    'C5H5CN': 70,
    'HOBr': 71,
    'CH3Br': 72,
    'CF4': 73,
    'SO3': 74,
    'Ne': 75,
    'Ar': 76,
    'COCl2': 77,
    'SO': 78,
    'H2SO4': 79,
    'e–': 80,
    'H3+': 81,
    'FeH': 82,
    'AlO': 83,
    'AlCl': 84,
    'AlF': 85,
    'AlH': 86,
    'BeH': 87,
    'C2': 88,
    'CaF': 89,
    'CaH': 90,
    'H–': 91,
    'CaO': 92,
    'CH': 93,
    'CH3': 94,
    'CH3F': 95,
    'CN': 96,
    'CP': 97,
    'CrH': 98,
    'HD+': 99,
    'HeH+': 100,
    'KCl': 101,
    'KF': 102,
    'LiCl': 103,
    'LiF': 104,
    'LiH': 105,
    'LiH+': 106,
    'MgF': 107,
    'MgH': 108,
    'MgO': 109,
    'NaCl': 110,
    'NaF': 111,
    'NaH': 112,
    'NH': 113,
    'NS': 114,
    'OH+': 115,
    'cis-P2H2': 116,
    'trans-P2H2': 117,
    'PH': 118,
    'PN': 119,
    'PO': 120,
    'PS': 121,
    'ScH': 122,
    'SH': 123,
    'SiH': 124,
    'SiH2': 125,
    'SiH4': 126,
    'SiO': 127,
    'SiS': 129,
    'TiH': 130,
    'Cl2': 131,
    'ClO2': 132,
    'O(1D)': 133,
    'N': 134,
    'N(2D)': 135
}


gas_info = {
    "1": {
        "name": "H2O",
        "label": "H$_2$O",
        "isotope": {
            "1": {
                "name": 'H2(16O)',
                "mass": 18.0,
                "label": 'H$_2^{16}$O'
            },
            "2": {
                "name": 'H2(18O)',
                "mass": 20.0,
                "label": 'H$_2^{18}$O'
            },
            "3": {
                "name": 'H2(17O)',
                "mass": 19.0,
                "label": 'H$_2^{17}$O'
            },
            "4": {
                "name": 'HDO',
                "mass": 19.0,
                "label": 'HDO'
            },
        },
        "mmw": 18.
    },
    "2": {
        "name": "CO2",
        "label": "CO$_2$",
        "isotope": {
            "1": {
                "name": 'CO2',
                "mass": 44.0,
                "label": '$^{12}$C$^{16}$O$_2$'
            },
            "2": {
                "name": '(13C)O2',
                "mass": 45.0,
                "label": '$^{13}$C$^{16}$O$_2$'
            },
            "3": {
                "name": '(18O)CO',
                "mass": 46.0,
                "label": '$^{18}$O$^{12}$C$^{16}$O'
            },
            "4": {
                "name": '(17O)CO',
                "mass": 45.0,
                "label": '$^{17}$O$^{12}$C$^{16}$O'
            },
        },
        "mmw": 44.
    },
    "3": {
        "name": "O3",
        "label": "O$_3$",
        "isotope": {
            "1": {
                "name": 'O3',
                "mass": 48.0,
                "label": '$^{16}$O$_3$'
            },
            "2": {
                "name": '(18O)O2',
                "mass": 50.0,
                "label": '$^{18}$O$^{16}$O$_2$'
            },
            "3": {
                "name": '(18O)2O',
                "mass": 52.0,
                "label": '$^{18}$O$_2^{16}$O'
            },
            "4": {
                "name": '(18O)3',
                "mass": 54.0,
                "label": '$^{18}$O$_3$'
            },
        },
        "mmw": 48.
    },
    "5": {
        "name": "CO",
        "label": "CO",
        "isotope": {
            "1": {
                "name": 'CO',
                "mass": 28.0,
                "label": '$^{12}$C$^{16}$O'
            },
            "2": {
                "name": '(13C)O',
                "mass": 29.0,
                "label": '$^{13}$C$^{16}$O'
            },
            "3": {
                "name": 'C(18O)',
                "mass": 30.0,
                "label": '$^{12}$C$^{18}$O'
            },
            "4": {
                "name": 'C(17O)',
                "mass": 29.0,
                "label": '$^{12}$C$^{17}$O'
            },
        },
        "mmw": 28.
    },
    "7": {
        "name": "O2",
        "label": "O$_2$",
        "isotope": {
            "1": {
                "name": 'O2',
                "mass": 32.0,
                "label": '$^{16}$O$_2$'
            },
            "2": {
                "name": '(18O)O',
                "mass": 34.0,
                "label": '$^{18}$OO'
            },
            "3": {
                "name": '(18O)2',
                "mass": 36.0,
                "label": '$^{18}$O$_2$'
            },
        },
        "mmw": 32.
    },
    "8": {
        "name": "NO",
        "label": "NO",
        "isotope": {
            "1": {
                "name": 'NO',
                "mass": 30.0,
                "label": '$^{14}$N$^{16}$O'
            },
            "2": {
                "name": '(15)NO',
                "mass": 31.0,
                "label": '$^{15}$N$^{16}$O'
            },
            "3": {
                "name": 'N(18O)',
                "mass": 32.0,
                "label": '$^{14}$N$^{18}$O'
            },
            "4": {
                "name": 'N(17O)',
                "mass": 31.0,
                "label": '$^{14}$N$^{17}$O'
            },
        },
        "mmw": 30.
    },
    "10": {
        "name": "NO2",
        "label": "NO$_2$",
        "isotope": {
            "1": {
                "name": 'NO2',
                "mass": 46.0,
                "label": '$^{14}$N$^{16}$O$_2$'
            },
            "2": {
                "name": '(15)NO2',
                "mass": 47.0,
                "label": '$^{15}$N$^{16}$O$_2$'
            },
            "3": {
                "name": 'N(18O)O',
                "mass": 48.0,
                "label": '$^{14}$N$^{18}$O$^{16}$O'
            },
            "4": {
                "name": 'N(18O)2',
                "mass": 50.0,
                "label": '$^{14}$N$^{18}$O$_2$'
            },
        },
        "mmw": 46.
    },
    "13": {
        "name": "OH",
        "label": "OH",
        "isotope": {
            "1": {
                "name": 'OH',
                "mass": 17.0,
                "label": '$^{16}$OH'
            },
            "2": {
                "name": '(18O)H',
                "mass": 19.0,
                "label": '$^{18}$OH'
            },
            "3": {
                "name": 'OD',
                "mass": 18.0,
                "label": '$^{16}$OD'
            },
        },
        "mmw": 17.
    },
    "22": {
        "name": "N2",
        "label": "N$_2$",
        "isotope": {
            "1": {
                "name": 'N2',
                "mass": 28.0,
                "label": '$^{14}$N$_2$'
            },
            "2": {
                "name": '(15)NN',
                "mass": 29.0,
                "label": '$^{15}$N$^{14}$N'
            },
            "3": {
                "name": '(15N)2',
                "mass": 30.0,
                "label": '$^{15}$N$_2$'
            },
        },
        "mmw": 28.
    },
    "25": {
        "name": "H2O2",
        "label": "H$_2$O$_2$",
        "isotope": {
            "1": {
                "name": 'H2(16O)2',
                "mass": 34.0,
                "label": 'H$_2^{16}$O$_2$'
            },
            "2": {
                "name": 'H2(18O)O',
                "mass": 36.0,
                "label": 'H$_2^{18}$O$^{16}$O'
            },
            "3": {
                "name": 'H2(18O)2',
                "mass": 38.0,
                "label": 'H$_2^{18}$O$_2$'
            },
            "4": {
                "name": 'HDO2',
                "mass": 35.0,
                "label": 'HD$^{16}$O$_2$'
            },
            "5": {
                "name": 'D2O2',
                "mass": 36.0,
                "label": 'D$_2^{16}$O$_2$'
            },
        },
        "mmw": 34.
    },
    "39": {
        "name": "H2",
        "label": "H$_2$",
        "isotope": {
            "1": {
                "name": 'H2',
                "mass": 2.0,
                "label": 'H$_2$'
            },
            "2": {
                "name": 'HD',
                "mass": 3.0,
                "label": 'HD'
            },
            "3": {
                "name": 'D2',
                "mass": 4.0,
                "label": 'D$_2$'
            },
        },
        "mmw": 2.
    },
    "44": {
        "name": "HO2",
        "label": "HO$_2$",
        "isotope": {
            "1": {
                "name": 'HO2',
                "mass": 33.0,
                "label": 'H$^{16}$O$_2$'
            },
            "2": {
                "name": 'H(18O)O',
                "mass": 35.0,
                "label": 'H$^{18}$O$^{16}$O'
            },
            "3": {
                "name": 'H(18O)2',
                "mass": 37.0,
                "label": 'H$^{18}$O$_2$'
            },
            "4": {
                "name": 'DO2',
                "mass": 34.0,
                "label": 'D$^{16}$O$_2$'
            },
        },
        "mmw": 33.
    },
    "45": {
        "name": "O",
        "label": "O",
        "isotope": {
            "1": {
                "name": 'O',
                "mass": 16.0,
                "label": '$^{16}$O'
            },
            "2": {
                "name": '(18O)',
                "mass": 18.0,
                "label": '$^{18}$O'
            },
            "3": {
                "name": '(17O)',
                "mass": 17.0,
                "label": '$^{17}$O'
            },
        },
        "mmw": 16.
    },
    "48": {
        "name": "H",
        "label": "H",
        "isotope": {
            "1": {
                "name": 'H',
                "mass": 1.0,
                "label": 'H'
            },
            "2": {
                "name": 'D',
                "mass": 2.0,
                "label": 'D'
            },
        },
        "mmw": 1.
    },
    "76": {
        "name": "Ar",
        "label": "Ar",
        "isotope": {
            "1": {
                "name": 'Ar',
                "mass": 40.0,
                "label": '$^{40}$Ar'
            },
            "2": {
                "name": '(36Ar)',
                "mass": 36.0,
                "label": '$^{36}$Ar'
            },
            "3": {
                "name": '(38Ar)',
                "mass": 38.0,
                "label": '$^{38}$Ar'
            },
        },
        "mmw": 40.
    },
    "133": {
        "name": "O(1D)",
        "label": "O(1D)",
        "isotope": {
            "1": {
                "name": 'O(1D)',
                "mass": 16.0,
                "label": '$^{16}$O($^1D$)'
            },
            "2": {
                "name": '(18O)(1D)',
                "mass": 18.0,
                "label": '$^{18}$O($^1D$)'
            },
            "3": {
                "name": '(17O)(1D)',
                "mass": 17.0,
                "label": '$^{17}$O($^1D$)'
            },
        },
        "mmw": 16.
    },
    "134": {
        "name": "N",
        "label": "N",
        "isotope": {
            "1": {
                "name": 'N',
                "mass": 14.0,
                "label": '$^{14}$N'
            },
            "2": {
                "name": '(15N)',
                "mass": 15.0,
                "label": '$^{15}$N'
            },
        },
        "mmw": 14.
    },
    "135": {
        "name": "N(2D)",
        "label": "N(2D)",
        "isotope": {
            "1": {
                "name": 'N(2D)',
                "mass": 14.0,
                "label": '$^{14}$N($^2D$)'
            },
            "2": {
                "name": '(15N)(2D)',
                "mass": 15.0,
                "label": '$^{15}$N($^2D$)'
            },
        },
        "mmw": 14.
    },
}