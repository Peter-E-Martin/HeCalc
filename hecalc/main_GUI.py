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
from base_GUI import Ui_MainWindow
from os import path
import sys
from main import WorkerClass

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
        self.actionExit.triggered.connect(app.quit)
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