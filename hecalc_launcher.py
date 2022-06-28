# -*- coding: utf-8 -*-
# This file is part of HeCalc, which calculates (U-Th)/He dates and uncertainties
# Copyright (C) 2021 Peter E. Martin <pemartin92@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
This script is provided to allow command-line interface and GUI
support for HeCalc such that the package and GUI can be used without
the need for every user to incorporate individual functions in their
own code.

When this script is run directly, the user will be prompted with a
series of options:
-Whether to run in the command line or with a graphical user interface
    -If the GUI is selected, no additional command line steps are necessary
-What percent precision the user wishes to target for the mean of the Monte
Carlo results (assuming Monte Carlo modeling is chosen later)
-How many decimals to output. This option affects only the saved data file
and does not affect the statistical handling of the data
-Whether 235U was measured directly (as opposed to assumed based on 238U
measurement)-- this is very rare and should be chosen carefully
-What kinds of uncertainty modeling to perform (Monte Carlo and/or linear
uncertainty propagation)
-if Monte Carlo modeling is chosen, the user will also decide whether
to save out the histograms produced and whether to parameterize the histograms
"""

import hecalc
from hecalc.GUI.main_GUI import launch_GUI

if __name__ == '__main__':
    # Prompt user to decide whether to run GUI
    GUI_or_CLI = input('Run [1] in command line or [2] in GUI? ')
    
    # Run GUI is chosen
    if GUI_or_CLI == '2':
        launch_GUI()

    # Prompt user to decide on details of modeling
    elif GUI_or_CLI == '1':
        percent_precision = float(input('Enter desired % precision of Monte Carlo uncertainty: ')) # user-defined relative precision, in percent
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
