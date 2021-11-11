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
This module contains the explicit math to calculate linear
uncertainty propagation for the (U-Th)/He age system.

Two functions are provided for this calculation:
"date_uncertainty" and "date_uncertainty_with235"

"date_uncertainty" should be used in all cases where 235U
was not measured directly but instead inferred from U or 238U
concentration. If 235U was measured directly, "date_uncertainty_with235"
should be used instead.

A number of helper functions are also present in this module.
These explicitly calculate the first derivatives of the age equation
with respect to date to propagate uncertainty. They are not imported
directly for the HeCalc package, but may be used independently if needed.

A large number of arguments are available for this module.
Required positional arguments:
    He, t
Optional arguments:
U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147,
He_s, U238_s, U235_s (only if U235 was measured) ,Th232_s, Sm147_s,
Ft238_s, Ft235_s, Ft232_s, Ft147_s
where U238_s is uncertainty in U238 for example.
All radionuclides and uncertainties default to a value of 0,
while all Fts default to a value of 1.

Note that at least one radionuclide value must be specified.
'''

import numpy as np

l_U238 = 1.55125e-10
l_U235 = 9.8485e-10
l_Th232 = 4.9475e-11
l_Sm147 = 6.539e-12

def _age_func_prime(U238, U235, Th232, Sm147,
                    Ft238, Ft235, Ft232, Ft147,
                    t):
    return ((8*Ft238*l_U238*U238*np.e**(t*l_U238) +
             7*Ft235*l_U235*U235*np.e**(t*l_U235) +
             6*Ft232*l_Th232*Th232*np.e**(t*l_Th232) +
             Ft147*l_Sm147*Sm147*np.e**(t*l_Sm147)))

def _He_prime(U238, U235, Th232, Sm147,
              Ft238, Ft235, Ft232, Ft147,
              He, t):
    return (1/_age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t))

def _U238_prime_noU235(U238, U235, Th232, Sm147,
                       Ft238, Ft235, Ft232, Ft147,
                       He, t):
    return (
        -(8*Ft238*(np.e**(l_U238*t)-1)+(7/137.818)*Ft235*(np.e**(l_U235*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _U238_prime_withU235(U238, U235, Th232, Sm147,
                         Ft238, Ft235, Ft232, Ft147,
                         He, t):
    return (
        -(8*Ft238*(np.e**(l_U238*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _U235_prime(U238, U235, Th232, Sm147,
                Ft238, Ft235, Ft232, Ft147,
                He, t):
    return (
        -(7*Ft235*(np.e**(l_U235*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _Th232_prime(U238, U235, Th232, Sm147,
                 Ft238, Ft235, Ft232, Ft147,
                 He, t):
    return (
        -(6*Ft232*(np.e**(l_Th232*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _Sm147_prime(U238, U235, Th232, Sm147,
                 Ft238, Ft235, Ft232, Ft147,
                 He, t):
    return (
        -(Ft147*(np.e**(l_Sm147*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _Ft238_prime(U238, U235, Th232, Sm147,
                 Ft238, Ft235, Ft232, Ft147,
                 He, t):
    return (
        -(8*U238*(np.e**(l_U238*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _Ft235_prime(U238, U235, Th232, Sm147,
                 Ft238, Ft235, Ft232, Ft147,
                 He, t):
    return (
        -(7*U235*(np.e**(l_U235*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _Ft232_prime(U238, U235, Th232, Sm147,
                 Ft238, Ft235, Ft232, Ft147,
                 He, t):
    return (
        -(6*Th232*(np.e**(l_Th232*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def _Ft147_prime(U238, U235, Th232, Sm147,
                 Ft238, Ft235, Ft232, Ft147,
                 He, t):
    return (
        -(Sm147*(np.e**(l_Sm147*t)-1))/
        _age_func_prime(U238, U235, Th232, Sm147, Ft238, Ft235, Ft232, Ft147, t)
        )

def date_uncertainty(He, t, He_s=0,
                     U238=0, U235=None, Th232=0, Sm147=0,
                     Ft238=1, Ft235=1, Ft232=1, Ft147=1,
                     U238_s=0, Th232_s=0, Sm147_s=0,
                     Ft238_s=0, Ft235_s=0, Ft232_s=0, Ft147_s=0):
    '''
    Returns symmetrical uncertainty via linear uncertainty propagation
    for a given date with uncertainty in He, radionuclides, and/or Ft values.
    
    Accepts integers, floats, or numpy arrays.
    
    WARNING:
    The exact date is required; this function does *not* calculate
    date explicitly. If an incorrect date is provided, the uncertainty
    calculation will also be incorrect.
    
    ~THIS FUNCTION SHOULD ONLY BE USED IF 235U IS INFERRED FROM 238U~
    If U235 was measured directly, the function date_uncertainty_with235()
    should be used instead.
    U235 is provided as an optional argument only if one wishes
    to use a non-standard value for the natural U ratio. If U235
    was measured directly, it will have its own uncertainty value
    and date_uncertainty_with235() should be used instead.
    
    Note that at least one radionuclide value must be specified.
    
    Parameters
    ----------
    He : float or array-like, optional
        The amount of helium measured. The exact units are unimportant
        as long as they are consistent with the other radionuclides
    
    t : float or array-like, optional
        The date for the uncertainty being calculated. This should be
        accurately calculated (e.g., using the get_date() function)
    
    He_s : float or array-like, optional
        The uncertainty in He, in the same units
        
    U238 : float or array-like, optional
        The amount of 238U measured
        
    U235 : float or array-like, optional
        The amount of 2385U measured. If no 235U is provided, it is assumed
        to be present at the standard ratio of 137.818 238U/235U
        
    Th232 : float or array-like, optional
        The amount of Th232 measured
        
    Sm147 : float or array-like, optional
        The amount of Sm147 measured
    
    Ft238 : float or array-like, optional
        The alpha ejection correction for 238U. Should be between 0-1.
        
    Ft235 : float or array-like, optional
        The alpha ejection correction for 235U. Should be between 0-1.
        
    Ft232 : float or array-like, optional
        The alpha ejection correction for 232Th. Should be between 0-1.
        
    Ft147 : float or array-like, optional
        The alpha ejection correction for 147Sm. Should be between 0-1.
    
    U238_s : float or array-like, optional
        The uncertainty in 238U

    Th232_s : float or array-like, optional
        The uncertainty in Th232
        
    Sm147_s : float or array-like, optional
        The uncertainty in Sm147
        
    Ft238_s : float or array-like, optional
        The uncertainty in Ft238
        
    Ft235_s : float or array-like, optional
        The uncertainty in Ft235
        
    Ft232_s : float or array-like, optional
        The uncertainty in Ft232
        
    Ft147_s : float or array-like, optional
        The uncertainty in Ft147
        
    Returns
    -------
    float or array-like of uncertainty values
    '''
    if U235 is None:
        U235 = U238/137.818
    return (
        (_He_prime(U238, U235, Th232, Sm147,
                   Ft238, Ft235, Ft232, Ft147,
                   He, t)*He_s)**2+
        (_U238_prime_noU235(U238, U235, Th232, Sm147,
                            Ft238, Ft235, Ft232, Ft147,
                            He, t)*U238_s)**2+
        (_Th232_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Th232_s)**2+
        (_Sm147_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Sm147_s)**2+
        (_Ft238_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft238_s)**2+
        (_Ft235_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft147_s)**2+
        (_Ft232_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft232_s)**2+
        (_Ft147_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft147_s)**2
        )**(1/2)

def date_uncertainty_with235(He, t, He_s=0,
                             U238=0, U235=0, Th232=0, Sm147=0,
                             Ft238=1, Ft235=1, Ft232=1, Ft147=1,
                             U238_s=0, U235_s=0, Th232_s=0, Sm147_s=0,
                             Ft238_s=0, Ft235_s=0, Ft232_s=0, Ft147_s=0):
    '''
    Returns symmetrical uncertainty for a given date with uncertainty
    in He, radionuclides, and/or all Ft values.
    
    Accepts integers, floats, or numpy arrays.
    
    WARNING:
    The exact date is required; this function does *not* calculate
    date explicitly. If an incorrect date is provided, the uncertainty
    calculation will also be incorrect.
    
        
    ~THIS FUNCTION SHOULD ONLY BE USED IF 235U HAS BEEN MEASURED DIRECTLY~
    If 235U is instead inferred from 238U measurement, the function
    "date_uncertainty" should be used instead.
    
    Note that at least one radionuclide value must be specified.
    
    Parameters
    ----------
    He : float or array-like, optional
        The amount of helium measured. The exact units are unimportant
        as long as they are consistent with the other radionuclides
    
    t : float or array-like, optional
        The date for the uncertainty being calculated. This should be
        accurately calculated (e.g., using the get_date() function)
    
    He_s : float or array-like, optional
        The uncertainty in He, in the same units
        
    U238 : float or array-like, optional
        The amount of 238U measured
        
    U235 : float or array-like, optional
        The amount of 2385U measured. If no 235U is provided, it is assumed
        to be present at the standard ratio of 137.818 238U/235U
        
    Th232 : float or array-like, optional
        The amount of Th232 measured
        
    Sm147 : float or array-like, optional
        The amount of Sm147 measured
    
    Ft238 : float or array-like, optional
        The alpha ejection correction for 238U. Should be between 0-1.
        
    Ft235 : float or array-like, optional
        The alpha ejection correction for 235U. Should be between 0-1.
        
    Ft232 : float or array-like, optional
        The alpha ejection correction for 232Th. Should be between 0-1.
        
    Ft147 : float or array-like, optional
        The alpha ejection correction for 147Sm. Should be between 0-1.
    
    U238_s : float or array-like, optional
        The uncertainty in 238U
        
    U235_s : float or array-like, optional
        The uncertainty in 235U, if measured.

    Th232_s : float or array-like, optional
        The uncertainty in Th232
        
    Sm147_s : float or array-like, optional
        The uncertainty in Sm147
        
    Ft238_s : float or array-like, optional
        The uncertainty in Ft238
        
    Ft235_s : float or array-like, optional
        The uncertainty in Ft235
        
    Ft232_s : float or array-like, optional
        The uncertainty in Ft232
        
    Ft147_s : float or array-like, optional
        The uncertainty in Ft147
        
    Returns
    -------
    float or array-like of uncertainty values
    '''
    return (
        (_He_prime(U238, U235, Th232, Sm147,
                   Ft238, Ft235, Ft232, Ft147,
                   He, t)*He_s)**2+
        (_U238_prime_withU235(U238, U235, Th232, Sm147,
                              Ft238, Ft235, Ft232, Ft147,
                              He, t)*U238_s)**2+
        (_U235_prime(U238, U235, Th232, Sm147,
                     Ft238, Ft235, Ft232, Ft147,
                     He, t)*U235_s)**2+
        (_Th232_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Th232_s)**2+
        (_Sm147_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Sm147_s)**2+
        (_Ft238_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft238_s)**2+
        (_Ft235_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft147_s)**2+
        (_Ft232_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft232_s)**2+
        (_Ft147_prime(U238, U235, Th232, Sm147,
                      Ft238, Ft235, Ft232, Ft147,
                      He, t)*Ft147_s)**2
        )**(1/2)