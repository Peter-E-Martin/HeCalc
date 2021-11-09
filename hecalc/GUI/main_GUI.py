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
Created on Sun Sep 26 08:42:43 2021

@author: Peter
"""

from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal
from .base_GUI import Ui_MainWindow
from os import path
import sys
import time
import re
import numpy as np
from ..main import _get_cols, _load_file, _sample_loop, _make_excel
from xlrd import XLRDError
from tabulate import tabulate

class GUI_App(Ui_MainWindow, QObject):

    cancel_signal = pyqtSignal(str)
    sheet_name = pyqtSignal(str)

    def __init__(self, window):
        super(GUI_App, self).__init__()
        self.setupUi(window)
        self.retranslateUi(window)
        self.file = None
        self.saveAs = None
        self.sheet_name = None
        self.tabWidget.setCurrentIndex(0)
        self.actionExit.triggered.connect(QtWidgets.QApplication.quit)#app.quit())
        self.actionHelp.triggered.connect(lambda: self.helpButton())
        self.actionAbout.triggered.connect(lambda: self.aboutButton())
        self.inputFileButton.clicked.connect(lambda: self.openFileNamesDialog())
        self.outputFileButton.clicked.connect(lambda: self.saveFileDialog())
        self.monteCarloCheckBox.toggled.connect(lambda: self.monteCarlo())
        self.U235_measuredCheck.toggled.connect(lambda: self.u235meas())
        self.outputHistCheck.toggled.connect(lambda: self.outputHistograms())
        self.cancelButton.clicked.connect(lambda: self.Canceled())
        self.goButton.clicked.connect(lambda: self.GO(self.file,
                                                      self.saveAs))

    def u235meas(self):
        if self.U235_measuredCheck.isChecked():
            self.u235inBox.setDisabled(False)
            self.u235sinBox.setDisabled(False)
        else:
            self.u235inBox.setDisabled(True)
            self.u235sinBox.setDisabled(True)

    def monteCarlo(self):
        if self.monteCarloCheckBox.isChecked():
            self.precisionEnter.setDisabled(False)
            self.outputHistCheck.setDisabled(False)
        else:
            self.precisionEnter.setDisabled(True)
            self.outputHistCheck.setDisabled(True)
            self.outputHistCheck.setChecked(False)
            
    def outputHistograms(self):
        if self.outputHistCheck.isChecked():
            self.pat.setEnabled(True)
        else:
            self.pat.setDisabled(True)
            self.pat.setChecked(False)

    def helpButton(self):
        helpBut = QtWidgets.QMessageBox()
        helpBut.setFixedSize(1000, 10000)
        helpBut.setWindowTitle("How to Use CU TRaIL Data Reduction")
        f = open('HTMLhelpmsg.txt', "r")
        msg = f.read()
        f.close()
        helpBut.setText(msg)
        helpBut.setStandardButtons(QtWidgets.QMessageBox.Ok)
        helpBut.exec_()

    def aboutButton(self):
        aboutBut = QtWidgets.QMessageBox()
        aboutBut.setWindowTitle("About CU TRaIL Data Reduction")
        aboutBut.setText(
            '''<html><head/><body><p align=\"center\" style=\"font-size:14px\">
            v0.0 | Dec. 15 2020<br/>Contact: Peter E. Martin
            (peter.martin-2@colorado.edu)<br/>Written in Python 3.8 using
            PyQt5</p></body></html>'''
            )
        aboutBut.setStandardButtons(QtWidgets.QMessageBox.Ok)
        aboutBut.exec_()

    def openFileNamesDialog(self):
        file, _filter = (QtWidgets
                         .QFileDialog
                         .getOpenFileName(None, "Choose File to Reduce", '.', "(*.xlsx *.xls *.csv *.txt)"))
        if file != '':
            self.userUpdate.append('Input: ' + path.split(file)[-1] + '\n')
            self.file = file

    def saveFileDialog(self):
        saveAs, _ = (QtWidgets
                      .QFileDialog
                      .getSaveFileName(None, "Save As", '.', "(*.xlsx)"))
        if saveAs != '':
            self.userUpdate.append('Output: ' + saveAs + '\n')
            self.saveAs = saveAs

    def GO(self, file, saveAs):
        if self.tabWidget.currentIndex() == 0:
            self.manual_input = False
            monteCarlo = self.monteCarloCheckBox.isChecked()
            linear = self.linearCheckBox.isChecked()
            measured_U235 = self.U235_measuredCheck.isChecked()
            outputHistogram = self.outputHistCheck.isChecked()
            parameterize = self.pat.isChecked()
            if file is not None and saveAs is not None:
                decimals = int(self.decBox.value())
                precision = float(self.precisionEnter.value()) / 100
                manual_data = None
                self.decBox.setDisabled(True)
                self.inputFileButton.setDisabled(True)
                self.outputFileButton.setDisabled(True)
                self.outputHistCheck.setDisabled(True)
                self.precisionEnter.setDisabled(True)
                self.monteCarloCheckBox.setDisabled(True)
                self.linearCheckBox.setDisabled(True)
                self.U235_measuredCheck.setDisabled(True)
                self.pat.setDisabled(True)
                self.label_3.setDisabled(True)
                self.cancelButton.setDisabled(False)
                self.goButton.setDisabled(True)
                self.mainCall(file,
                              saveAs,
                              precision,
                              decimals,
                              parameterize,
                              outputHistogram,
                              monteCarlo,
                              linear,
                              measured_U235,
                              manual_data)
                
            elif file is None:
                self.userUpdate.append('Please choose an input file')
            elif saveAs is None:
                self.userUpdate.append('Please choose an output folder')
        elif self.tabWidget.currentIndex() == 1:
            self.manual_input = True
            monteCarlo = self.monteCarloCheckBox.isChecked()
            linear = self.linearCheckBox.isChecked()
            measured_U235 = self.U235_measuredCheck.isChecked()
            decimals = int(self.decBox.value())
            precision = float(self.precisionEnter.value()) / 100
            saveAs = None
            outputHistogram = False
            parameterize = False
            manual_data = [float(self.heInBox.text()),
                           float(self.hesInBox.text()),
                           float(self.uInBox.text()),
                           float(self.usInBox.text()),
                           float(self.u235inBox.text()),
                           float(self.u235sinBox.text()),
                           float(self.thInBox.text()),
                           float(self.thsInBox.text()),
                           float(self.smInBox.text()),
                           float(self.smsInBox.text()),
                           float(self.ft238InBox.text()),
                           float(self.ft238sInBox.text()),
                           float(self.ft235InBox.text()),
                           float(self.ft235sInBox.text()),
                           float(self.ft232Inbox.text()),
                           float(self.ft232sInBox.text()),
                           float(self.ft147Inbox.text()),
                           float(self.ft147sInBox.text())]
            self.decBox.setDisabled(True)
            self.inputFileButton.setDisabled(True)
            self.outputFileButton.setDisabled(True)
            self.outputHistCheck.setDisabled(True)
            self.precisionEnter.setDisabled(True)
            self.monteCarloCheckBox.setDisabled(True)
            self.linearCheckBox.setDisabled(True)
            self.U235_measuredCheck.setDisabled(True)
            self.pat.setDisabled(True)
            self.label_3.setDisabled(True)
            self.goButton.setDisabled(True)
            self.mainCall(file,
                          saveAs,
                          precision,
                          decimals,
                          outputHistogram,
                          parameterize,
                          monteCarlo,
                          linear,
                          measured_U235,
                          manual_data)

                
    def progressUpdate(self, value):
        self.progressBar.setValue(value)

    def updateStatus(self, status):
        self.userUpdate.append(status)

    def updateText(self, text):
        self.progressLabel.setText(text)

    def Canceled(self):
        self.worker.cancel = True
        self.cancelButton.setText('Cancel')
        self.cancelButton.setDisabled(True)

    def loopDone(self):
        self.cancelButton.setDisabled(True)
        self.goButton.setDisabled(False)
        self.decBox.setDisabled(False)
        self.inputFileButton.setDisabled(False)
        self.outputFileButton.setDisabled(False)
        self.outputHistCheck.setDisabled(False)
        self.precisionEnter.setDisabled(False)
        self.monteCarloCheckBox.setDisabled(False)
        self.linearCheckBox.setDisabled(False)
        self.U235_measuredCheck.setDisabled(False)
        self.pat.setDisabled(False)
        self.label_3.setDisabled(False)

    def completed(self):
        self.cancelButton.setText('Reset')
        self.cancelButton.clicked.connect(lambda: self.reset())

    def reset(self):
        self.cancelButton.setText('Cancel')
        self.cancelButton.setDisabled(True)
        self.goButton.setDisabled(False)
        self.decBox.setDisabled(False)
        self.inputFileButton.setDisabled(False)
        self.outputFileButton.setDisabled(False)
        self.outputHistCheck.setDisabled(False)
        self.precisionEnter.setDisabled(False)
        self.monteCarloCheckBox.setDisabled(False)
        self.linearCheckBox.setDisabled(False)
        self.U235_measuredCheck.setDisabled(False)
        self.pat.setDisabled(False)
        self.label_3.setDisabled(False)

    def manual_completed(self):
        self.cancelButton.setText('Cancel')
        self.cancelButton.setDisabled(True)
        self.goButton.setDisabled(False)
        self.decBox.setDisabled(False)
        self.inputFileButton.setDisabled(False)
        self.outputFileButton.setDisabled(False)
        self.outputHistCheck.setDisabled(False)
        self.precisionEnter.setDisabled(False)
        self.monteCarloCheckBox.setDisabled(False)
        self.linearCheckBox.setDisabled(False)
        self.U235_measuredCheck.setDisabled(False)
        self.pat.setDisabled(False)
        self.label_3.setDisabled(False)

    def open_sheet_dialog(self):
        sheet , ok = QtWidgets.QInputDialog.getText(MainWindow,'Sheet Name Entry','<html style="font-size:12pt;text-align:center;">Please enter the name of the<br/>excel sheet containing data</html>')
        if ok:
            self.sheet_name = sheet
            self.GO(self.file, self.saveAs, self.outputHistogram)
        if not ok:
            self.reset()

    def mainCall(self,
                 file,
                 saveAs,
                 precision,
                 decimals,
                 parameterize,
                 outputHistogram,
                 monteCarlo,
                 linear,
                 measured_U235,
                 manual_data):
        if not self.manual_input:
            self.worker_thread = QThread()
            self.worker = WorkerClass(file,
                                      saveAs,
                                      precision,
                                      decimals,
                                      outputHistogram,
                                      parameterize,
                                      monteCarlo,
                                      linear,
                                      measured_U235,
                                      self.manual_input,
                                      manual_data,
                                      self.sheet_name)
            self.worker.moveToThread(self.worker_thread)
            self.worker.progress.connect(self.progressUpdate)
            self.worker.status.connect(self.updateStatus)
            self.worker.text.connect(self.updateText)
            self.worker_thread.started.connect(self.worker.mainLoop)
            self.worker.finished.connect(self.worker_thread.quit)
            self.worker.finished.connect(self.completed)
            self.worker.error.connect(self.worker_thread.quit)
            self.worker.stop.connect(self.reset)
            self.worker.sheet_entry.connect(self.open_sheet_dialog)
            self.worker.sheet_entry.connect(self.worker_thread.terminate)
            self.worker.loop_done.connect(self.loopDone)
            self.worker.loop_done.connect(self.worker_thread.quit)
            self.worker_thread.start()
        elif self.manual_input:
            self.worker_thread = QThread()
            self.worker = WorkerClass(file,
                                      saveAs,
                                      precision,
                                      decimals,
                                      parameterize,
                                      outputHistogram,
                                      monteCarlo,
                                      linear,
                                      measured_U235,
                                      self.manual_input,
                                      manual_data,
                                      self.sheet_name)
            self.worker.moveToThread(self.worker_thread)
            self.worker.progress.connect(self.progressUpdate)
            self.worker.status.connect(self.updateStatus)
            self.worker.text.connect(self.updateText)
            self.worker_thread.started.connect(self.worker.mainLoop)
            self.worker.finished.connect(self.worker_thread.quit)
            self.worker.finished.connect(self.manual_completed)
            self.worker.error.connect(self.worker_thread.quit)
            self.worker.stop.connect(self.reset)
            self.worker.loop_done.connect(self.loopDone)
            self.worker.loop_done.connect(self.worker_thread.quit)
            self.worker_thread.start()

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

def launch_GUI():
    try:
        import pyi_splash
        pyi_splash.close()
    except:
        pass
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('AlphaDecayIcon.ico'))
    MainWindow = QtWidgets.QMainWindow()
    ui = GUI_App(MainWindow)
    MainWindow.show()
    app.exec_()

if __name__ == "__main__":
    try:
        import pyi_splash
        pyi_splash.close()
    except:
        pass
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('AlphaDecayIcon.ico'))
    MainWindow = QtWidgets.QMainWindow()
    ui = GUI_App(MainWindow)
    MainWindow.show()
    app.exec_()