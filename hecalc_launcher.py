# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 10:22:06 2021

@author: Peter
"""

import hecalc
from hecalc.GUI.main_GUI import launch_GUI

if __name__ == '__main__':
    GUI_or_CLI = input('Run [1] in command line or [2] in GUI? ')
    
    if GUI_or_CLI == '2':
        launch_GUI()
        
    elif GUI_or_CLI == '1':
        percent_precision = float(input('Enter desired % precision of Monte Carlo mean: ')) # user-defined relative precision, in percent
        decimals =  int(input('Enter number of decimals to output: ')) # number of decimals to report
        measured_U235 = input('Was 235U measured directly? [Y/N]')
        if measured_U235 == 'Y' or measured_U235 == 'y':
            measured_U235 = True
        else:
            measured_U235 = False
        mc = input('Run Monte Carlo? [Y/N]')
        if mc == 'Y' or mc == 'y':
            monteCarlo = True
        else:
            monteCarlo = False
        lin = input('Run linear error propagation? [Y/N]')
        if lin == 'Y' or lin == 'y':
            linear = True
        else:
            linear = False
        if monteCarlo:
            histograms = input('Output histograms? [Y/N] ') # Y -or- N to determine histogram output
            if histograms == 'Y' or histograms == 'y':
                histograms = True
                parameterize = input('Parameterize histogram distributions? [Y/N] ') # user decides whether to parameterize for skewed distributions
                if parameterize == 'y' or parameterize == 'Y':
                    parameterize = True
                else:
                    parameterize = False
            else:
                histograms = False
                parameterize = False
        else:
            histograms = False
            parameterize = False
        hecalc.main.hecalc_main(file=None, saveAs=None, percent_precision=percent_precision,
                                decimals=decimals, measured_U235=measured_U235, monteCarlo=monteCarlo,
                                linear=linear, histograms=histograms, parameterize=parameterize)