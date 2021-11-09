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

'''
Functions to return a (U-Th)/He date or multiple dates.
All functions will read in radionuclides, He, and Ft values
as integers, floats, or numpy arrays (if multiple dates are desired).

All functions return date(s) in years.
'''

import numpy as np

l_U238 = 1.55125e-10
l_U235 = 9.8485e-10
l_Th232 = 4.9475e-11
l_Sm147 = 6.539e-12
l_tot = [l_U238, l_U235, l_Th232, l_Sm147]

def iterated_date(He, t_guess,
                  U238=0, U235=None, Th232=0, Sm147=0,
                  Ft238=1, Ft235=1, Ft232=1, Ft147=1):
    '''
    Return (U-Th)/He date accurately calculated using the Newton-Raphson
    method. The iteration is run until successive iterations differ
    by less than one year. If an array of dates is being calculated,
    the iteration will be run until the averaged difference between iterations
    is less than one year, while eliminating any array elements that
    do not converge.

    Takes int, float, or numpy arrays as input.
    
    Returned date(s) are in years.
    
    There is a restriction on this function of 100 iterations.
    When calculations take longer than this, it is generally because
    the He age equation is not well-behaved (e.g., a radionuclide value
    is negative). In this case, a NaN value is returned.
    
    Required inputs: He, estimated t, and at least one radioniclide value.
    Optional inputs: U238, U235, Th232, Sm147, Ft238,
                     Ft235, Ft232, Ft147
    '''
    # if U235 is left with the default value, assume it must be
    # based on a constant ratio to U238
    if U235 is None:
        U235 = U238/137.818
    try:
        # First convert ints or floats to arrays of len(1) 
        # to accomodate array datatype
        single_date = False
        if type(t_guess) != np.ndarray:
            single_date = True
            t_guess = np.array([t_guess])
        # create list of current and previous time (with initial
        # value of 0) to track amount of change between iterations
        times = [np.zeros(len(t_guess)), t_guess]
        n = 0
        ignores = []
        # restricting to n<100 is really just a safeguard. In practice, most iterated
        # dates should require <10 cycles to calculate
        while np.nanmean(np.delete(abs(times[1]-times[0]),ignores)) > 1 and n<100:
            n+=1
            times[0] = times[1]
            f_t = (8*Ft238*U238 * (np.exp(times[1]*l_U238) - 1) +
                   7*Ft235*U235 * (np.exp(times[1]*l_U235) - 1) +
                   6*Ft232*Th232 * (np.exp(times[1]*l_Th232) - 1) +
                   Ft147*Sm147 * (np.e**(times[1]*l_Sm147) - 1))-He
            f_t_prime = (8*Ft238*l_U238*U238*np.exp(times[1]*l_U238) +
                         7*Ft235*l_U235*U235*np.exp(times[1]*l_U235) +
                         6*Ft232*l_Th232*Th232*np.exp(times[1]*l_Th232) +
                         Ft147*l_Sm147*Sm147*np.exp(times[1]*l_Sm147))
            new_date = times[1] - f_t/f_t_prime
            times[1] = new_date
            # give the method three cycles to converge, and then stop evaluating
            # samples that aren't converging, defined as changes of >0.1% of
            # value for each cycle after 3 cycles
            if n>3:
                rats = abs(times[1]/times[0]-1)
                ignores = np.where(rats>0.001)[0]
    
        final_date = np.delete(times[1], ignores)
        final_date = final_date[~np.isnan(final_date)]
        final_date = final_date[~np.isinf(final_date)]
        # If no dates converged, return NaN
        if len(final_date) == 0:
            return np.nan
        # Convert back to float if array was not passed originally
        if single_date:
            final_date = float(final_date[0])
        return final_date
    except:
        return np.nan

def meesters_dunai(He,
                   U238=0, U235=None, Th232=0, Sm147=0,
                   Ft238=1, Ft235=1, Ft232=1, Ft147=1):
    '''
    Calculate a (U-Th)/He date using the method of Meesters and Dunai (2005).
    Takes floats, ints, or numpy arrays for
    U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, and He.
    
    Users are referred to Meesters and Dunai (2005) for details and
    limitations on this method. In short, this method is very efficient
    and reasonably accurate for dates <500 Ma. Older dates (or dates
    where high accuracy is needed) should be calculated using the
    iterated_date() function instead.
                                                            
    If non-physical values are passed to this function (e.g., negative
    production rates of He), a linearized estimate of date is returned
    instead using the equation t = He/TotalProduction.
    
    Returns non-iterative date in years.
    
    Required inputs: He, and at least one radioniclide value.
    Optional inputs: U238, U235, Th232, Sm147, Ft238,
                     Ft235, Ft232, Ft147
    '''
    if U235 is None:
        U235 = U238/137.818
    # Calculate individual production values
    p238 = 8*l_U238*U238*Ft238
    p235 = 7*l_U235*U235*Ft235
    p232 = 6*l_Th232*Th232*Ft235
    p147 = l_Sm147*Sm147*Ft147
    # Generate total and weighted mean production
    p_tot = [p238, p235, p232, p147]
    P = p238 + p235 + p232 + p147
    l_wm = sum([l_tot[i]*p_tot[i]/P for i in range(len(p_tot))])
    try:
        return 1/l_wm*np.log((l_wm/P)*He+1)
    # If the Meesters and Dunai method fails, use
    # a linearized version of the age equaiton
    except ValueError:
        return He/sum(p_tot)
        
# Function to get age given He, radioisotopes, and Fts
def get_date(He,
             U238=0, U235=None, Th232=0, Sm147=0,
             Ft238=1, Ft235=1, Ft232=1, Ft147=1):
    '''
    Calculate raw (uncorrected) and alpha-ejection corrected
    dates. These are returned in a list with uncorrected
    date(s) first and corrected date(s) second.
    
    Returned dates are in years.
    
    Takes int, float, or numpy arrays as input.
    
    Required inputs: He, and at least one radioniclide value.
    Optional inputs: U238, U235, Th232, Sm147, Ft238,
                     Ft235, Ft232, Ft147
    '''
    if U235 is None:
        U235 = U238/137.818
    # Use M&D 2005 age as first guess for iterative age
    # start by calculating raw age using 1 for all Ft
    t_app_raw = meesters_dunai(He, U238, U235, Th232, Sm147)
    if type(U238) == np.ndarray:
        if len(U238) == 0:
            return [np.nan, np.nan]
    t_it = iterated_date(He, t_app_raw,
                         U238, U235, Th232, Sm147)
    
    # Repeat with Fts to get corrected age
    t_app_corr = meesters_dunai(He,
                                U238, U235, Th232, Sm147,
                                Ft238, Ft235, Ft232, Ft147)
    if type(U238) == np.ndarray:
        if len(U238) == 0:
            return [np.nan, np.nan]
    t_it_corr = iterated_date(He, t_app_corr,
                              U238, U235, Th232, Sm147,
                              Ft238, Ft235, Ft232, Ft147)
    return {'raw date': t_it, 'corrected date': t_it_corr}