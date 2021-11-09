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
This module provides funcitonality for either the GUI (through a
threaded worker class using Qthread or through a simple function
to read in data and call the core functionality of HeCalc.

This module is not necessary to run HeCalc if a different data
read-in is desired; the core HeCalc functions may be called indepenently
by other scripts.

If main is used (through the GUI by running the main_GUI.py script
or command line by running this script), it will accept
an excel, csv, or tan-delimited text file. This file must have the headers:
Sample, 238U, 232Th, 147Sm, 4He, 238Ft, 235Ft, 232Ft, and 147Ft.
If 235U was measured directly, this column must also be present.

If run through the command line, each parameter will be set using
the input() built-in function. Files are loaded using tkinter.
The command line version of HeCalc does not support "manual input"
as the GUI does. This is because it is simpler to write a quick
script that calls the individual date calculation and uncertainty
propagation functions contained in this package than it is to include
this option every time HeCalc is run in the command line.

Uncertainty in each value is required (even if the uncertainty is 0).
The uncertainty must follow its respective column; the header for
these uncertainty columns have no labeling requirements.

The output from this module is an excel file with a header containing
the input file path, and (if Monte Carlo uncertainty popagation is used)
the user-requested precision.

The column headers for the output include sample name, time to run, and
corrected and uncorrected dates at a minimum. If linear uncertainty
propagation is requested, 1-s uncertainty columns are added.
If Monte Carlo uncertainty propagation is requested, +/- 68% confidence
intervals are added along with % skew.
'''

from PyQt5.QtCore import QObject, pyqtSignal
import time
import pandas as pd
from xlrd import XLRDError
import numpy as np
import re
from tabulate import tabulate
from openpyxl import styles, utils, Workbook#, formatting
# from sample_calc import sample_loop
import warnings
import tkinter
from tkinter.filedialog import askopenfilename, asksaveasfilename
from date_calculation import get_date
from linear_propagation import date_uncertainty, date_uncertainty_with235
from monte_carlo import monte_carlo

warnings.filterwarnings("ignore")

def _get_cols(linear, monteCarlo, parameterize):
    save_columns = ['Sample',
                    'Raw date',
                    'Mean raw date',
                    'Linear raw uncertainty',
                    ' +68% CI raw',
                    ' -68% CI raw',
                    '% Skewness of raw distribution',
                    'Corrected date',
                    'Mean corrected date',
                    'Linear corrected uncertainty',
                    ' +68% CI corrected',
                    ' -68% CI corrected',
                    '% Skewness of corrected distribution',
                    'Number of Monte Carlo simulations']
    if monteCarlo and not linear:
        save_columns = ['Sample',
                        'Raw date',
                        'Mean raw date',
                        ' +68% CI raw',
                        ' -68% CI raw',
                        '% Skewness of raw distribution',
                        'Corrected date',
                        'Mean corrected date',
                        ' +68% CI corrected',
                        ' -68% CI corrected',
                        '% Skewness of corrected distribution',
                        'Number of Monte Carlo simulations']
    elif linear and not monteCarlo:
        save_columns = ['Sample',
                        'Raw date',
                        'Linear raw uncertainty',
                        'Corrected date',
                        'Linear corrected uncertainty']
    elif not monteCarlo and not linear:
        save_columns = ['Sample',
                        'Raw date',
                        'Corrected date']
    if parameterize:
        ins_idx = save_columns.index('% Skewness of raw distribution')
        add_cols = ['raw fit a','raw fit u','raw fit s']
        save_columns = save_columns[:ins_idx+1]+add_cols+save_columns[ins_idx+1:]
        ins_idx = save_columns.index('% Skewness of corrected distribution')
        add_cols = ['corrected fit a','corrected fit u','corrected fit s']
        save_columns = save_columns[:ins_idx+1]+add_cols+save_columns[ins_idx+1:]
    return save_columns

def _load_file(file, measured_U235, sheet):
    cols = ['mol 238U',
            'mol 232Th',
            'mol 147Sm',
            'mol 4He',
            '238Ft',
            '235Ft',
            '232Ft',
            '147Ft']
    if measured_U235:
        cols.insert(1, 'mol 235U')
    
    if file[-5:] == '.xlsx' or file[-4:] == '.xls':
        if not sheet:
            try:
                data_load = pd.read_excel(file) # Load the data sheet
                assert set(cols).issubset(list(data_load.columns))
            except AssertionError:
                return None
        elif sheet:
            try:
                data_load = pd.read_excel(file, sheet_name = sheet) # Load the data sheet
                assert set(cols).issubset(list(data_load.columns))
            except:
                return None
    if file[-4:] == '.csv':
        data_load = pd.read_csv(file)
    if file[-4:] == '.txt':
        data_load = pd.read_csv(file, delimiter="\t")
    
    col_load = ['Sample']
    for c in cols:
        col_load.append(c)
        col_load.append(data_load.columns[data_load.columns.get_loc(c)+1])
    
    final_cols = ['Sample']
    for i in range(len(cols)):
        final_cols.append(cols[i])
        final_cols.append(u'\u00B1 '+cols[i])
    
    data = data_load[col_load]
    data.columns = final_cols
    short_names = {name:name.replace('mol ','') for name in data.columns if 'mol ' in name}
    data = data.rename(columns=short_names)
    
    return data

def _make_excel(save_out, save_columns, file, monteCarlo, precision_user, saveAs):
    # Generate the excel file
    book = Workbook()
    
    # Start by adding histograms to their own sheet and
    # then deleting them to clean up the save_out variable
    if 'raw histogram' in save_out:
        hist_sheet = book.active
        hist_sheet.title = 'Histogram Output'
        
        keys = {0: 'raw histogram',
                1: 'corrected histogram'}
        xy = {0: 'bin centers',
              1: 'values'}
        for col in range(int(len(save_out['Sample']))):
            for i in range(2):
                for n in range(2):
                    hist_sheet.cell(row=1, column=col*4+i*2+n+1,
                                    value= save_out['Sample'][col]+'\n'+keys[i]+'\n'+xy[n])
                    for row in range(len(save_out[keys[i]][col][n])):
                        hist_sheet.cell(row=row+2, column=col*4+i*2+n+1,
                                        value=save_out[keys[i]][col][n][row])
        del save_out['raw histogram']
        del save_out['corrected histogram']
        for header in hist_sheet["1:1"]:
            header.font = styles.Font(bold=True)
            header.alignment = styles.Alignment(wrapText=True)
            if 'raw' in header.value:
                hist_sheet.column_dimensions[utils.get_column_letter(header.column)].width = 13
            if 'corrected' in header.value:
                hist_sheet.column_dimensions[utils.get_column_letter(header.column)].width = 11
       
    # Create a list to format for Excel
    colHead_row = "2:2"
    save_final = []
    for s in save_out:
        save_final.append(save_out[s])
    save_final = list(zip(*save_final))
    save_final.insert(0, save_columns)
    save_final.insert(0, ['input = ' + file])
    if monteCarlo:
        save_final.insert(0, ['precision = ' + str(round(precision_user,3))+'%'])
        colHead_row = "3:3"
    
    # Make a new sheet for the uncertainty results
    sheet = book.active
    if sheet.title == 'Histogram Output':
        book.create_sheet('Uncertainty Output', 0)
        output = book.active
    else:
        output = book.active
        output.title = 'Uncertainty Output'
    
    # Write the list to the excel file
    for row in range(len(save_final)):
        for col in range(len(save_final[row])):
            output.cell(row=row+1, column=col+1, value=save_final[row][col])
    
    # Format the column headers
    for cell in output[colHead_row]:
        cell.font = styles.Font(bold=True)
        cell.alignment = styles.Alignment(wrapText=True)
        if '% Skewness' in cell.value:
            output.column_dimensions[utils.get_column_letter(cell.column)].width = 20
        if 'Linear' in cell.value:
            output.column_dimensions[utils.get_column_letter(cell.column)].width = 15
        if 'Corrected date' in cell.value:
            output.column_dimensions[utils.get_column_letter(cell.column)].width = 10
        if 'simulations' in cell.value:
            output.column_dimensions[utils.get_column_letter(cell.column)].width = 17
        if 'Mean corrected' in cell.value:
            output.column_dimensions[utils.get_column_letter(cell.column)].width = 14
    output.row_dimensions[3].height = 31
    
    try:
        save_out['histogram']
    except:
        pass
    
    return book

def _linear_propagation(nominal_t, data, measured_U235):
    if measured_U235:
        linear_raw_uncertainty = date_uncertainty_with235(
            data['4He'], nominal_t['raw date'], data[u'\u00B1 4He'],
            data['238U'], data['235U'], data['232Th'], data['147Sm'],
            U238_s=data[u'\u00B1 238U'], U235_s=data[u'\u00B1 235U'],
            Th232_s=data[u'\u00B1 232Th'], Sm147_s=data[u'\u00B1 147Sm']
                                                          )
        linear_corr_uncertainty = date_uncertainty_with235(
            data['4He'], nominal_t['corrected date'], data[u'\u00B1 4He'],
            data['238U'], data['235U'], data['232Th'], data['147Sm'],
            data['238Ft'], data['235Ft'], data['232Ft'], data['147Ft'],
            data[u'\u00B1 238U'], data[u'\u00B1 235U'],
            data[u'\u00B1 232Th'], data[u'\u00B1 147Sm'],
            data[u'\u00B1 238Ft'], data[u'\u00B1 235Ft'],
            data[u'\u00B1 232Ft'], data[u'\u00B1 147Ft']
                                                           )
    else:
        linear_raw_uncertainty = date_uncertainty(
            data['4He'], nominal_t['raw date'], data[u'\u00B1 4He'],
            data['238U'], data['238U']/137.818, data['232Th'], data['147Sm'],
            U238_s=data[u'\u00B1 238U'], Th232_s=data[u'\u00B1 232Th'],
            Sm147_s=data[u'\u00B1 147Sm']
                                                  )
        linear_corr_uncertainty = date_uncertainty(
            data['4He'], nominal_t['corrected date'], data[u'\u00B1 4He'],
            data['238U'], data['238U']/137.818, data['232Th'], data['147Sm'],
            data['238Ft'], data['235Ft'], data['232Ft'], data['147Ft'],
            data[u'\u00B1 238U'], data[u'\u00B1 232Th'], data[u'\u00B1 147Sm'],
            data[u'\u00B1 238Ft'], data[u'\u00B1 235Ft'],
            data[u'\u00B1 232Ft'], data[u'\u00B1 147Ft']
                                                   )
    return {'raw unc': linear_raw_uncertainty, 'corr unc': linear_corr_uncertainty}

def _sample_loop(save_out, sample_data, measured_U235, linear, monteCarlo,
                 histograms, parameterize, decimals, precision):
    
    U235 = sample_data['238U']/137.818
    if measured_U235:
        U235 = sample_data['235U']
    
    nominal_t = get_date(sample_data['4He'], sample_data['238U'], U235,
                         sample_data['232Th'], sample_data['147Sm'],
                         sample_data['238Ft'], sample_data['235Ft'],
                         sample_data['232Ft'], sample_data['147Ft'])
    save_out['Raw date'].append(round(nominal_t['raw date']/1e6, decimals))
    save_out['Corrected date'].append(round(nominal_t['corrected date']/1e6, decimals))
    
    # Calls helper function _linear_propagation, which unpacks relevant
    # data and calls the correct functions for linear uncertainty propagation
    linear_uncertainty = _linear_propagation(nominal_t, sample_data, measured_U235)
    if linear:
        save_out['Linear raw uncertainty'].append(round(linear_uncertainty['raw unc']/1e6,decimals))
        save_out['Linear corrected uncertainty'].append(round(linear_uncertainty['corr unc']/1e6,decimals))
    
    if monteCarlo:
        s_est = linear_uncertainty['raw unc']
        mean_est = nominal_t['corrected date']
        mc_number = s_est**2/(precision*mean_est)**2
        if mc_number < 5:
            mc_number = 5
        else:
            mc_number = int(mc_number)
        save_out['Number of Monte Carlo simulations'].append(mc_number)
            
        if not measured_U235:
            U235 = None
            U235_s = None
        elif measured_U235:
            U235 = sample_data['235U']
            U235_s = sample_data[u'\u00B1 235U']
            
        mc_results = monte_carlo(
            mc_number, sample_data['4He'], sample_data[u'\u00B1 4He'],
            sample_data['238U'], U235, sample_data['232Th'], sample_data['147Sm'],
            sample_data['238Ft'], sample_data['235Ft'], sample_data['232Ft'], sample_data['147Ft'],
            sample_data[u'\u00B1 238U'], U235_s,
            sample_data[u'\u00B1 232Th'], sample_data[u'\u00B1 147Sm'],
            sample_data[u'\u00B1 238Ft'], sample_data[u'\u00B1 235Ft'],
            sample_data[u'\u00B1 232Ft'], sample_data[u'\u00B1 147Ft'],
            histogram=histograms, parameterize=parameterize
                                 )
        
        for ft in ['raw', 'corrected']:
            save_out['Mean '+ft+' date'].append(round(mc_results[ft+' date']['mean']/1e6,decimals))
            save_out[' +68% CI '+ft].append(round(mc_results[ft+' date']['+68% CI']/1e6,decimals))
            save_out[' -68% CI '+ft].append(round(mc_results[ft+' date']['-68% CI']/1e6,decimals))
            save_out['% Skewness of '+ft+' distribution'].append(round(mc_results[ft+' date']['% Skew'],decimals))
            if histograms:
                mc_results[ft+' date']['histogram'][0] = np.around(mc_results[ft+' date']['histogram'][0]/1e6,decimals)
                mc_results[ft+' date']['histogram'][1] = mc_results[ft+' date']['histogram'][1]
                save_out[ft+' histogram'].append(mc_results[ft+' date']['histogram'])
                if parameterize:
                    save_out[ft+' fit a'].append(mc_results[ft+' date']['a'])
                    try:
                        save_out[ft+' fit u'].append(round(mc_results[ft+' date']['u']/1e6,decimals))
                        save_out[ft+' fit s'].append(round(mc_results[ft+' date']['s']/1e6,decimals))
                    except:
                        save_out[ft+' fit u'].append(mc_results[ft+' date']['u'])
                        save_out[ft+' fit s'].append(mc_results[ft+' date']['s'])
    return save_out

class WorkerClass(QObject):
    
    status = pyqtSignal(str)
    text = pyqtSignal(str)
    progress = pyqtSignal(float)
    finished = pyqtSignal()
    error = pyqtSignal()
    stop = pyqtSignal()
    loop_done = pyqtSignal()
    sheet_entry = pyqtSignal()

    def __init__(self,
                 file,
                 saveAs,
                 prec,
                 decimals,
                 histograms,
                 parameterize,
                 monteCarlo,
                 linear,
                 measured_U235,
                 manual_input,
                 manual_data,
                 sheet_name):
        super(WorkerClass, self).__init__()
        self.file = file
        self.saveAs = saveAs
        self.prec = prec
        self.decimals = decimals
        self.histograms = histograms
        self.parameterize = parameterize
        self.monteCarlo = monteCarlo
        self.linear = linear
        self.measured_U235 = measured_U235
        self.manual_input = manual_input
        self.manual_data = manual_data
        self.cancel = False
        self.sheet_name = sheet_name

    def mainLoop(self):
        save_columns = _get_cols(self.linear, self.monteCarlo, self.parameterize)
        save_out = {c: [] for c in save_columns}
        precision_user = self.prec*100

        if not self.manual_input:
            if self.histograms:
                save_out['raw histogram'] = []
                save_out['corrected histogram'] = []
            total_time = time.time()
            
            # Load data
            data = _load_file(self.file, self.measured_U235, self.sheet_name)
            if self.sheet_name is not None:
                self.status.emit('Sheet not found. Please try again.')
            if data is None:
                self.sheet_entry.emit()
                return
            
            try:  
                # This is the main loop that calculates progress and calls the data reduction code
                times = []
                for i in range(len(data)):
                    while not self.cancel:
                        # Update user of new sample being analyzed
                        self.progress.emit(100 * (len(save_out["Sample"])/len(data)))
                        self.status.emit('Now working on sample ' + str(data['Sample'].iloc[i]))
    
                        # Estimate time remaining and alert user
                        if len(times) == 0:
                            self.text.emit(
                                f'{len(save_out["Sample"])} of {len(data)} ' +
                                'samples analyzed\nEstimating time remaining')
                        elif round((len(data)-i)*np.mean(times)/60, 1) > 60:
                            t_left = round((len(data)-i)*np.mean(times)/3600, 1)
                            self.text.emit(
                                f'{len(save_out["Sample"])} of {len(data)} ' +
                                f'samples analyzed\nEst. {t_left} hr remaining')
                        else:
                            t_left = round((len(data)-i)*np.mean(times)/60, 1)
                            self.text.emit(
                                f'{len(save_out["Sample"])} of {len(data)} ' +
                                f'samples analyzed\nEst. {t_left} min remaining')
    
                        # Check that no illegal characters are being used and replace with underscores
                        if (re.sub('[^\w\-_]', '_', data['Sample'].iloc[i]) != 
                                data['Sample'].iloc[i]):
                            self.status.emit('sample name ' +
                                             str(data['Sample'].iloc[i]) +
                                             ' changed to ' +
                                             str(re.sub('[^\w\-_]', '_', data['Sample'].iloc[i])) +
                                             ' to avoid illegal characters in file name')
                            data['Sample'].iloc[i] = re.sub('[^\w\-_]', '_', data['Sample'].iloc[i])
                        save_out['Sample'].append(data['Sample'].iloc[i])
                        
                        sample_data = data.loc[[i]].to_dict(orient='record')[0]
                        save_out = _sample_loop(save_out,
                                                sample_data,
                                                self.measured_U235,
                                                self.linear,
                                                self.monteCarlo,
                                                self.histograms,
                                                self.parameterize,
                                                self.decimals,
                                                self.prec)
                        break
    
                book = _make_excel(save_out, save_columns, self.file,
                                   self.monteCarlo, precision_user, self.saveAs)
    
                # When user cancels, alert user to cancelled status and kill the worker thread
                if self.cancel:
                    self.progress.emit(0.)
                    self.text.emit('Data reduction canceled.\nNo data saved.')
                    self.status.emit('Data reduction canceled\n')
                    self.loop_done.emit()
    
                # Save the excel workbook, update user, and kill the worker thread
                if self.cancel is False:
                    try:
                        book.save(self.saveAs)
                        self.progress.emit(100.)
                        self.text.emit(f'{len(data)} samples analyzed.')
                        self.status.emit('Done! Total time to run: ' +
                                         str(round((time.time()-total_time)/60, 2)) +
                                         ' minutes\n')
                        self.finished.emit()
    
                    # If the user has a workbook with the same name open, alert and kill thread
                    except PermissionError:
                        self.progress.emit(100.)
                        self.text.emit(f'{len(data)} samples analyzed.')
                        self.status.emit('\nCould not overwrite existing file.'+
                                         '\nPlease close the open file and rerun.')
                        self.stop.emit()
                        self.error.emit()
    
            # If the excel sheet is formatted incorrectly, alert user and kill thread
            except (XLRDError, KeyError):
                self.status.emit('Please use correctly formatted data input.' +
                                  '\nSee the help tab or documentation for more detail')
                self.stop.emit()
                self.error.emit()
                
        if self.manual_input:
            self.histograms = False
            self.parameterize = False
            try:
                if self.measured_U235:
                    inputs = ['','4He', '238U', '235U', '232Th', '147Sm', '238Ft', '235Ft', '232Ft', '147Ft']
                    vals = ['value']+self.manual_data[::2]
                    uncs = [u'\u00B1']+self.manual_data[1::2]
                    table = [inputs,vals,uncs]
                    rotated = list(zip(*table))
                    self.status.emit(tabulate(rotated, tablefmt = 'presto', headers="firstrow"))
                else:
                    inputs = ['','4He', 'U', '232Th', '147Sm', '238Ft', '235Ft', '232Ft', '147Ft']
                    vals = ['value']+[self.manual_data[0]]+[self.manual_data[2]]+self.manual_data[6::2]
                    uncs = [u'\u00B1']+[self.manual_data[1]]+[self.manual_data[3]]+self.manual_data[7::2]
                    table = [inputs,vals,uncs]
                    rotated = list(zip(*table))
                    self.status.emit(tabulate(rotated, tablefmt = 'presto', headers="firstrow"))
                
                data = {'Sample': 'Manually Entered',
                        '238U': self.manual_data[2],
                        u'\u00B1 '+'238U': self.manual_data[3],
                        '235U': self.manual_data[4],
                        u'\u00B1 '+'235U': self.manual_data[5],
                        '232Th': self.manual_data[6],
                        u'\u00B1 '+'232Th': self.manual_data[7],
                        '147Sm': self.manual_data[8],
                        u'\u00B1 '+'147Sm': self.manual_data[9],
                        '4He': self.manual_data[0],
                        u'\u00B1 '+'4He': self.manual_data[1],
                        '238Ft': self.manual_data[10],
                        u'\u00B1 '+'238Ft': self.manual_data[11],
                        '235Ft': self.manual_data[12],
                        u'\u00B1 '+'235Ft': self.manual_data[13],
                        '232Ft': self.manual_data[14],
                        u'\u00B1 '+'232Ft': self.manual_data[15],
                        '147Ft': self.manual_data[16],
                        u'\u00B1 '+'147Ft': self.manual_data[17]}
                output_lists = _sample_loop(save_out,
                                            data,
                                            self.measured_U235,
                                            self.linear,
                                            self.monteCarlo,
                                            self.histograms,
                                            self.parameterize,
                                            self.decimals,
                                            self.prec)
                output = {}
                for d in output_lists:
                    if len(output_lists[d])>0:
                        output[d] = output_lists[d][0]
                key_change = ['Raw date', 'Raw\ndate',
                              'Mean raw date', 'Mean\nraw\ndate',
                              'Linear raw uncertainty', 'Linear\nraw\n'+u'\u00B1',
                              ' +68% CI raw', '+68%\nCI\nraw',
                              ' -68% CI raw', '-68%\nCI\nraw',
                              '% Skewness of raw distribution', '% raw\nskew',
                              'Corrected date', 'Corr.\ndate',
                              'Mean corrected date', 'Mean\ncorrected\ndate',
                              'Linear corrected uncertainty', 'Linear\ncorr.\n'+u'\u00B1',
                              ' +68% CI corrected', '+68%\nCI\ncorr.',
                              ' -68% CI corrected', '-68%\nCI\ncorr.',
                              '% Skewness of corrected distribution', '% corr.\nskew']
                for i in range(0,len(key_change),2):
                    if key_change[i] in output.keys():
                        output[key_change[i+1]] = output.pop(key_change[i])
                if self.monteCarlo:
                    self.status.emit('\n')
                    self.status.emit('Number of Monte Carlo simulations: '+str(output['Number of Monte Carlo simulations']))
                self.status.emit('\n')
                raws = [[k, v] for k, v in list(output.items()) if any(r in k for r in ['raw', 'Raw'])]
                rotate_raws = list(zip(*raws))
                self.status.emit(tabulate(rotate_raws))
                self.status.emit('\n')
                corrs = [[k, v] for k, v in list(output.items()) if any(r in k for r in ['corr', 'Corr'])]
                rotate_corrs = list(zip(*corrs))
                self.status.emit(tabulate(rotate_corrs))
                self.status.emit('\n')
                self.finished.emit()
            except:
                self.status.emit('Date cannot be calculated from these values\n')
                self.finished.emit()

def hecalc(file=None, saveAs=None, percent_precision=0.01, decimals=5, measured_U235=False,
           monteCarlo=True, linear=True, histograms=False, parameterize=False):
    '''
    This is the main module to run HeCalc as a standalone program
    through the command line interface. The user will first be
    prompted to select a number of options:
    1. the precision of the Monte Carlo mean determines the number
       of Monte Carlo simulations to perform; the number input here
       is the relative precision in percent (i.e., for a sample
       with a mean age of 50 Ma, a 1% precision would be +/- 0.5 Ma)
    2. the number of decimals to report simply dictates the
       level to which the code will round the result; this
       has no bearing on the statistical portions of the code
    3. Whether 235U was measured directly. This should only
       be selected if 235 is not assumed from 238U measurement
    4. Whether to run the Monte Carlo simulation.
    5. Whether to run linear uncertainty propagation
    6a. Whether to generate histograms. These are saved out as txt files
       with bin centers and relative height.
    6b. Whether to parameterize the histogram. In this case, the
        histogram will be fit to a skewnormal distribution as
        a first-order approximation of the date as a continuous
        random variable.
    
    The input file for this module must have the following headers:
    Sample, 238U, 232Th, 147Sm, 4He, 238Ft, 235Ft, 232Ft, and 147Ft
    Sample is the sample name
    238U, 232Th, 147Sm, and 4He are amounts of each nuclide and
    can be in any unit (e.g., mol, atoms, mol/g)
    238Ft, 235Ft, 232Ft, and 147Ft are the isotope-specific
    Ft values for each nuclide. see Ketcham et al., 2011 for
    how these may be calculated.
    
    required arguments: none
    optional arguments: file
    file: a full path to a correctly formatted excel, txt, or csv
          file on the local machine where HeCalc is being run.
    '''
    # Get user-defined variables through command line entry
    precision = percent_precision/100 # convert from % to proportion
    
    # create the list of columns to be used as headers in the output
    # these vary based on which (if any) uncertainty propagation is
    # chosen by the user.
    save_columns = _get_cols(linear,monteCarlo,parameterize)
    save_out = {c: [] for c in save_columns}
    
    if histograms:
        save_out['raw histogram'] = []
        save_out['corrected histogram'] = []
    
    # Get the file containing the data if no file was provided
    if file == None:
        root = tkinter.Tk()
        try:
            file = askopenfilename(title = 'Choose file to reduce') # brings up window to select data file
            root.wm_withdraw()
            assert any(ext in file[-4:] for ext in ['xlsx', '.xls', '.csv', '.txt'])
        except AssertionError:
            print('Input file must be excel, csv, or tab-delimited text file')
            file = askopenfilename(title = 'Choose file to reduce') # brings up window to select data file
            root.wm_withdraw()
    
    if saveAs == None:
        root = tkinter.Tk()
        saveAs = asksaveasfilename(filetypes=[("Excel files", "*.xlsx .*xls")],
                                     defaultextension=".xlsx")
        root.wm_withdraw()
    
    # Load in the data. The helper function _load_file
    # goes through a series of checks to make sure the 
    # file has the appropriate information. If not,
    # an input is requested in case an excel workbook
    # with multiple sheets has been provided
    data = _load_file(file, measured_U235, None)
    while data is None:
        sheet = input('Name of sheet to read data from: ')
        data = _load_file(file, measured_U235, sheet)
        if data is None:
            print('Ensure columns are correctly labeled')
    
    # This is the main loop that calculates progress
    # and calls the individual data reduction modules
    for i in range(len(data)):
        # Update user of new sample being analyzed
        print(('Now working on sample ' + str(data['Sample'].iloc[i])))
        save_out['Sample'].append(data['Sample'].iloc[i])

        # Check that no illegal characters are being used and replace with underscores
        if re.sub('[^\w\-_]', '_', data['Sample'].iloc[i]) != data['Sample'].iloc[i]:
            data['Sample'].iloc[i] = re.sub('[^\w\-_]', '_', data['Sample'].iloc[i])
        
        sample_data = data.loc[[i]].to_dict(orient='record')[0]
        
        save_out = _sample_loop(save_out, sample_data, measured_U235, linear, monteCarlo,
                                histograms, parameterize, decimals, precision)

    book = _make_excel(save_out, save_columns, file, monteCarlo, percent_precision, saveAs)
    
    # Save the excel workbook
    book.save(saveAs)
    
if __name__ == "__main__":
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
    hecalc(file=None, saveAs=None, percent_precision=percent_precision,
           decimals=decimals, measured_U235=measured_U235, monteCarlo=monteCarlo,
           linear=linear, histograms=histograms, parameterize=parameterize)