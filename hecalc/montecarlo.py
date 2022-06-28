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
This module takes a dataset and runs Monte Carlo uncertainty modeling.

It contains functions for the Monte Carlo modeling and to generate
histograms and parameterize those histograms if a user so desires.
"""

from scipy.stats import skewnorm as sk
import numpy as np
from .date_calculation import get_date

def make_histogram(mc, parameterize):
    '''
    Create histograms and potentially parameterize them as skew-normal
    probability distribution functions.
    
    Parameters
    ----------
    mc : Array-like
        The list of numbers to place in a histogram
        
    parameterize : bool
        whether to parameterize the histogram
        
    Returns
    -------
    Generates a dictionary with the histogram and potentially parameters.
    the histogram is returned as a list of the bin centers and a list of
    the number of samples within the bin. The parameters areassociated
    with the key "fit distribution", and are a list of the three parameters
    (skewness, location, scale) to calculate a skew-normal distribution
    '''
    # Generate histogram with bins centered
    bins = int(len(mc)/1000)
    if bins>1000:
        bins = 1000
    elif bins<10:
        bins = 10
    histo = list(np.histogram(mc, bins=bins))[::-1]
    histo[0] = (histo[0][1:] + histo[0][:-1]) / 2
    
    # Get stats for skewnormal function if the user called for parameterization
    if parameterize:
        # Subsample Monte Carlo results if lots of simulations were run.
        # skewnorm.fit is very slow with large samples, and the accuracy of
        # the fit is perfectly acceptable with 1e5 samples
        if len(mc) > int(1e5):
            sample = np.random.choice(mc, int(1e5), replace = False)
        else:
            sample = mc
        a, u, s = sk.fit(sample)
        # except:
        #     print('parameterization unsuccessful, returning only histogram')
        #     return {'histogram': histo,
        #             'fit distribution': [None, None, None]}
        return {'histogram': histo,
                'fit distribution': [a, u, s]}
    return {'histogram': histo}

def monte_carlo(mc_number, He, He_s=0,
                U238=0, U235=None, Th232=0, Sm147=0,
                Ft238=1, Ft235=1, Ft232=1, Ft147=1,
                U238_s=0, U235_s=None, Th232_s=0, Sm147_s=0,
                Ft238_s=0, Ft235_s=0, Ft232_s=0, Ft147_s=0,
                U238_U235_v=0, U238_Th232_v=0, U238_Sm147_v=0,
                U235_Th232_v=0, U235_Sm147_v=0, Th232_Sm147_v=0,
                Ft238_Ft235_v=0, Ft238_Ft232_v=0, Ft238_Ft147_v=0,
                Ft235_Ft232_v=0, Ft235_Ft147_v=0, Ft232_Ft147_v=0,
                histogram=False, parameterize=False):
    '''
    Performs Monte Carlo modeling for a (U-Th)/He dataset. Note that
    while all radionuclides are optional and default to 0, at least one
    *MUST* be >0.
    
    Parameters
    ----------
    mc_number : int
        The number of monte carlo cycles to perform
        
    He : float
        The amount of helium measured. The exact units are unimportant
        as long as they are consistent with the other radionuclides
    
    He_s : float, optional
        The uncertainty in He, in the same units
        
    U238 : float, optional
        The amount of 238U measured
        
    U238_s : float, optional
        The uncertainty in 238U
        
    U235 : float, optional
        The amount of 2385U measured. If no 235U is provided, it is assumed
        to be present at the standard ratio of 137.818 238U/235U
        
    U235_s : float, optional
        The uncertainty in 235U, if measured.
        
    Th232 : float, optional
        The amount of Th232 measured
        
    Th232_s : float, optional
        The uncertainty in Th232
        
    Sm147 : float, optional
        The amount of Sm147 measured
        
    Sm147_s : float, optional
        The uncertainty in Sm147
    
    Ft238 : float, optional
        The alpha ejection correction for 238U. Should be between 0-1.
        
    Ft238_s : float, optional
        The uncertainty in Ft238
        
    Ft235 : float, optional
        The alpha ejection correction for 235U. Should be between 0-1.
        
    Ft235_s : float, optional
        The uncertainty in Ft235
        
    Ft232 : float, optional
        The alpha ejection correction for 232Th. Should be between 0-1.
        
    Ft232_s : float, optional
        The uncertainty in Ft232
        
    Ft147 : float, optional
        The alpha ejection correction for 147Sm. Should be between 0-1.
        
    Ft147_s : float, optional
        The uncertainty in Ft147
        
    U238_U235_v : float, optional
        The covariance between U238 and U235
    
    U238_Th232_v : float, optional
        The covariance between U238 and Th232
    
    U238_Sm147_v : float, optional
        The covariance between U238 and Sm147
    
    U235_Th232_v : float, optional
        The covariance between U235 and Th232
    
    U235_Sm147_v : float, optional
        The covariance between U235 and Sm147
    
    Th232_Sm147_v : float, optional
        The covariance between Th232 and Sm147
    
    Ft238_Ft235_v : float, optional
        The covariance between Ft238 and Ft235
    
    Ft238_Ft232_v : float, optional
        The covariance between Ft238 and Ft232
    
    Ft238_Ft147_v : float, optional
        The covariance between Ft238 and Ft147
    
    Ft235_Ft232_v : float, optional
        The covariance between Ft235 and Ft232
    
    Ft235_Ft147_v : float, optional
        The covariance between Ft235 and Ft147
    
    Ft232_Ft147_v : float, optional
        The covariance between Ft232 and Ft147
        
    histogram : bool, optional
        Whether to generate a histogram for the Monte Carlo results
        
    parameterize : bool, optional
        Whether to parameterize the histogram with a skew-normal distribution.
        If histogram is False, parameterize will not be used
        
    Returns
    -------
    Produces a dictionary with entries for the raw and alpha-ejection corrected
    results. Each entry is a dictionary with mean, +/- 68% CI, and the number of
    successful cycles. If histogram is true, the histogram for each is returned 
    as a list of the bin centers and a list of the number of samples within the 
    bin. The parameters are associated with the key "fit distribution", and are 
    a list of the three parameters (skewness, location, scale) to calculate a 
    skew-normal distribution
    '''
    # Helium is always univariate, so generate normal sample here
    He = np.random.normal(He, He_s, mc_number)
    
    # Generate radionuclides covariance matrix (3x3 in most cases, 4x4 for U235 measured)
    # and corresponding multivariate normal arrays
    if U235 == None:
        rad_means = [U238, Th232, Sm147]
        rad_cov = [[U238_s**2, U238_Th232_v, U238_Sm147_v],
                   [U238_Th232_v, Th232_s**2, Th232_Sm147_v],
                   [U238_Sm147_v, Th232_Sm147_v, Sm147_s**2]]
        rads = np.random.multivariate_normal(rad_means, rad_cov, mc_number).T
        U238 = rads[0]
        Th232 = rads[1]
        Sm147 = rads[2]
        U235 = U238/137.818
    else:
        rad_means = [U238, U235, Th232, Sm147]
        rad_cov = [[U238_s**2, U238_U235_v, U238_Th232_v, U238_Sm147_v],
                   [U238_U235_v, U235_s**2, U235_Th232_v, U235_Sm147_v],
                   [U238_Th232_v, U235_Th232_v, Th232_s**2, Th232_Sm147_v],
                   [U238_Sm147_v, U235_Sm147_v, Th232_Sm147_v, Sm147_s**2]]
        rads = np.random.multivariate_normal(rad_means, rad_cov, mc_number).T
        U238 = rads[0]
        U235 = rads[1]
        Th232 = rads[2]
        Sm147 = rads[3]
    
    # Generate Ft covariance matrix and multivariate normal arrays
    Ft_means = [Ft238, Ft235, Ft232, Ft147]
    Ft_covs = [[Ft238_s**2, Ft238_Ft235_v, Ft238_Ft232_v, Ft238_Ft147_v],
               [Ft238_Ft235_v, Ft235_s**2, Ft235_Ft232_v, Ft235_Ft147_v],
               [Ft238_Ft232_v, Ft235_Ft232_v, Ft232_s**2, Ft232_Ft147_v],
               [Ft238_Ft147_v, Ft235_Ft147_v,Ft232_Ft147_v, Ft147_s**2]]
    Fts = np.random.multivariate_normal(Ft_means, Ft_covs, mc_number).T
    Ft238 = Fts[0]
    Ft235 = Fts[1]
    Ft232 = Fts[2]
    Ft147 = Fts[3]
    
    # Call the get_date function to calculate a date for each
    # of the randomly generated data
    MonteCarlo_t = get_date(He, U238, U235, Th232, Sm147,
                            Ft238, Ft235, Ft232, Ft147)
        
    # Create a dictionary to save results to
    mc_results = {'raw date': {},
                  'corrected date': {}}
    # Produce statistics and save
    for t in ['raw date', 'corrected date']:
        mean = np.mean(MonteCarlo_t[t])
        CIs_68 = np.percentile(MonteCarlo_t[t], [15.865, 84.135])
        CIs_95 = np.percentile(MonteCarlo_t[t], [15.865, 84.135])
        mc_results[t]['mean'] = mean
        mc_results[t]['+68% CI'] = CIs_68[1]
        mc_results[t]['-68% CI'] = CIs_68[0]
        mc_results[t]['+95% CI'] = CIs_95[1]
        mc_results[t]['-95% CI'] = CIs_95[0]
        # Report how many of the Monte Carlo simulations converged to
        # an age. For highly uncertain input data, this may be significantly
        # lower than the number of requested cycles.
        mc_results[t]['cycles'] = len(MonteCarlo_t[t])
        if histogram:
            histogram_out = make_histogram(MonteCarlo_t[t], parameterize)
            mc_results[t]['histogram'] = histogram_out['histogram']
            if parameterize:
                mc_results[t]['a'] = histogram_out['fit distribution'][0]
                mc_results[t]['u'] = histogram_out['fit distribution'][1]
                mc_results[t]['s'] = histogram_out['fit distribution'][2]
    
    return mc_results