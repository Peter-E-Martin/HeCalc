
[![DOI](https://zenodo.org/badge/407311456.svg)](https://zenodo.org/badge/latestdoi/407311456)

# HeCalc

This is a package written in Python 3.8 to reduce (U-Th)/He data, including statistically rigorous propagation of uncertainty. All uncertainties are input and output at the 1-sigma level.

## Installation

HeCalc is available through the python package index (PyPI), so you can install it using ```pip install hecalc```

### User Interface installation

HeCalc is also available as an executable which will run a self-contained Windows application (Mac version coming soon!) that does not require the use of Python. The executable is available with each HeCalc release, which can be found on the right-hand side of the GitHub repository page by clicking "Releases". Alternatively, the most recent download is available by [clicking here](https://github.com/Peter-E-Martin/HeCalc/releases/download/v0.2.1-alpha/HeCalc.exe).

If you instead wish to run the user interface through python, optional dependencies are available through the command ```pip install hecalc[gui]```. This optional installation mode allows hecalc to be run in the command line as described below. ***Warning***: ```pip install hecalc[gui]``` is incompadible with Anaconda. Anaconda's native PyQt version allows the HeCalc to be run, so if you wish to use HeCalc in an Anaconda environment ```pip install hecalc``` is the safest option.

## Running HeCalc

There are three main methods of running HeCalc:

- Command line
- Graphical User Interface (GUI)
- Through individual functions

The command line and GUI can both be triggered from the current codebase by running the hecalc_launcher.py script. These each will prompt the user for a series of decisions about the kind of data reduction to perform, the style of file to save out, and which dataset to read in. **All data uncertainties should be at the 1-sigma level.**

HeCalc functions both as an application and a Python package. This allows a user to write their own wrapper in Python to use the functions contained within HeCalc rather than using the main integrated HeCalc program (e.g., to avoid the necessity of reconfiguring data files for the input, or for smooth integration with existing lab data reduction schemes). Details on these functions are contained below.

### Running the GUI

Having installed HeCalc with ```pip install hecalc[gui]```, it is possible to run the GUI without having downloaded the executable. Simply run the commands:
```import hecalc.GUI as gui```
```gui.launch_GUI()```

Alternatively, an executable version of the GUI is available with each release. For those interested only in running the GUI, simply download the HeCalc.exe file from the most recent release. This will allow you to run HeCalc without any knowledge of Python (or even having Python installed on your computer).

## Input

Data input for HeCalc can be in .xlsx, .xls, .csv, or tab-delimited .txt form. The following columns **must** be present with these exact names:
Sample, mol 4He, mol 238U, mol 232Th, mol 147Sm,  238Ft, 235Ft, 232Ft, 147Ft
Each column must be followed by its 1-sigma uncertainty value, **even if that uncertainty is 0**. A typical header column will therefore look something like this:

Sample | mol 4He | ± | mol 238U | ± | mol 232Th | ± | mol 147Sm | ± | 238Ft | ± | 235Ft | ± | 232Ft | ± | 147Ft | ±

Correlated uncertainties in the radionuclides and Ft values can be accommodated by HeCalc. By default, uncertainties will be treated as independent. However, correlation coefficients for each known correlation can be included by adding columns with the naming convention "r value1-value2". For example, columns for known correlations between 238U and 232Th, 238Ft and 235Ft, and 238Ft and 232Ft would look like:

r 238U-232Th | r 238Ft-235Ft | r 238Ft-232Ft

Example files are included in the Test directory that can serve as template for data entry. These in includes options with both correlated and uncorrelated uncertainty.

## User Options

The user will be prompted to provide the following information for HeCalc:

 - The precision of the mean corrected Monte Carlo date, in percent relative to the date (i.e., for a sample with a mean age of 50 Ma, a 1% precision would be +/- 0.5 Ma). This determines the total number of Monte Carlo cycles to run.
 - The number of decimals to report in the output excel workbook. This simply dictates rounding of the output and does not affect the statistical results.
 - Whether 235U was directly measured, though this should be very rare given the most common lab setups
 - Whether to run Monte Carlo uncertainty propagation
 - Whether to run linear uncertainty propagation
 - Whether to generate histograms of the Monte Carlo simulations
 - Whether to parameterize the histograms (if generated) to produce an approximate skew-normal probability distribution for the data.

## Output

When HeCalc is fully run, it will produce an excel file with header for the user-requested precision (assuming that Monte Carlo uncertainty propagation was run) and the path for the source file. An example output file with all options selected is included in the Test directory. Each sample will then occupy a row with columns for the raw and alpha-ejection corrected values of:
date, mean date*, 1-sigma linear uncertainty**, +/- 68% confidence intervals*, % skewness*, fit a parameter***, fit u parameter***, fit s parameter***, and the number of Monte Carlo simulations*

\* Only if Monte Carlo uncertainty propagation was chosen
\** Only if Linear uncertainty propagation was chosen only
\*** Only if parameterization was chosen. a = shape (i.e., skewness), u = location, s = scale

If histogram generation was chosen with Monte Carlo uncertainty propagation, the histograms will be included in a second excel sheet labeled "Histogram Output", with raw and alpha ejection corrected histograms consisting of bin centers and the total number of samples in each bin.

## hecalc functions

If HeCalc is installed in site-packages (i.e., is downloaded as a Python package), several functions are available out of the box using the package name hecalc:

|Function|Use|
|--|--|
|hecalc.iterated_date()|Using input data and an estimate of date, uses the Newton-Raphson method to calculate exact dates|
|hecalc.meesters_dunai()|Generates a non-iterative date solution using the Meesters & Dunai 2005 method|
|hecalc.get_date|Combines the two above functions to provide raw and alpha ejection-corrected dates directly from a given dataset|
|hecalc.date_uncertainty()|Performs linear uncertainty propagation for a given dataset, assuming that 235U is calculated from the measurement of 238U. Uncertainties should be at the 1-sigma level.|
|hecalc.date_uncertainty_with235()|Performs the same linear uncertainty propagation but accounts for the fact that the 235U measurement is independent of the 238U measurement if it was measured directly. Uncertainties should be at the 1-sigma level.|
|hecalc.monte_carlo()|Runs Monte Carlo uncertainty propagation on a dataset assuming gaussian 1-sigma uncertainty, and outputs statistics and (if requested) the histogram and parameterized fit for the data|
|hecalc.hecalc_main()|This is a manual function to run exactly the same set of data reduction as the software version, though users can input specific options more flexibly|
