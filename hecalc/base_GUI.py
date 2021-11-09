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

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(836, 793)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(14)
        MainWindow.setFont(font)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(1)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(28, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.goButton = QtWidgets.QPushButton(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(16)
        self.goButton.setFont(font)
        self.goButton.setObjectName("goButton")
        self.horizontalLayout_3.addWidget(self.goButton)
        spacerItem1 = QtWidgets.QSpacerItem(48, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.cancelButton = QtWidgets.QPushButton(self.centralwidget)
        self.cancelButton.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.cancelButton.setFont(font)
        self.cancelButton.setObjectName("cancelButton")
        self.horizontalLayout_3.addWidget(self.cancelButton)
        spacerItem2 = QtWidgets.QSpacerItem(28, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem2)
        self.gridLayout_2.addLayout(self.horizontalLayout_3, 3, 0, 2, 1)
        self.progressBar = QtWidgets.QProgressBar(self.centralwidget)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setTextVisible(False)
        self.progressBar.setObjectName("progressBar")
        self.gridLayout_2.addWidget(self.progressBar, 2, 1, 2, 1)
        self.statusMessenger = QtWidgets.QLabel(self.centralwidget)
        self.statusMessenger.setAlignment(QtCore.Qt.AlignCenter)
        self.statusMessenger.setObjectName("statusMessenger")
        self.gridLayout_2.addWidget(self.statusMessenger, 0, 1, 1, 1)
        self.userUpdate = QtWidgets.QTextBrowser(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Courier New")
        font.setPointSize(12)
        self.userUpdate.setFont(font)
        self.userUpdate.setObjectName("userUpdate")
        self.gridLayout_2.addWidget(self.userUpdate, 1, 1, 1, 1)
        self.progressLabel = QtWidgets.QLabel(self.centralwidget)
        self.progressLabel.setText("")
        self.progressLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.progressLabel.setObjectName("progressLabel")
        self.gridLayout_2.addWidget(self.progressLabel, 4, 1, 1, 1)
        self.scrollArea_3 = QtWidgets.QScrollArea(self.centralwidget)
        self.scrollArea_3.setMinimumSize(QtCore.QSize(0, 100))
        self.scrollArea_3.setWidgetResizable(True)
        self.scrollArea_3.setObjectName("scrollArea_3")
        self.scrollAreaWidgetContents_3 = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_3.setGeometry(QtCore.QRect(0, 0, 367, 674))
        self.scrollAreaWidgetContents_3.setObjectName("scrollAreaWidgetContents_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.scrollAreaWidgetContents_3)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.label_6 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_6.setFont(font)
        self.label_6.setAlignment(QtCore.Qt.AlignCenter)
        self.label_6.setObjectName("label_6")
        self.gridLayout_3.addWidget(self.label_6, 0, 0, 1, 1)
        self.tabWidget = QtWidgets.QTabWidget(self.scrollAreaWidgetContents_3)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.tabWidget.setFont(font)
        self.tabWidget.setUsesScrollButtons(False)
        self.tabWidget.setObjectName("tabWidget")
        self.fileTab = QtWidgets.QWidget()
        self.fileTab.setObjectName("fileTab")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.fileTab)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.scrollArea = QtWidgets.QScrollArea(self.fileTab)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 323, 327))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.gridLayout_7 = QtWidgets.QGridLayout()
        self.gridLayout_7.setObjectName("gridLayout_7")
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem3, 2, 1, 1, 1)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem4)
        self.outputFileButton = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        self.outputFileButton.setMaximumSize(QtCore.QSize(16777215, 16777214))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.outputFileButton.setFont(font)
        self.outputFileButton.setObjectName("outputFileButton")
        self.horizontalLayout_10.addWidget(self.outputFileButton)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem5)
        self.gridLayout_7.addLayout(self.horizontalLayout_10, 3, 1, 1, 1)
        self.pat = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        self.pat.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.pat.setFont(font)
        self.pat.setObjectName("pat")
        self.gridLayout_7.addWidget(self.pat, 6, 1, 1, 1)
        spacerItem6 = QtWidgets.QSpacerItem(13, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_7.addItem(spacerItem6, 1, 2, 1, 1)
        spacerItem7 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem7, 4, 1, 1, 1)
        self.outputHistCheck = QtWidgets.QCheckBox(self.scrollAreaWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.outputHistCheck.setFont(font)
        self.outputHistCheck.setObjectName("outputHistCheck")
        self.gridLayout_7.addWidget(self.outputHistCheck, 5, 1, 1, 1)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem8)
        self.inputFileButton = QtWidgets.QPushButton(self.scrollAreaWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.inputFileButton.setFont(font)
        self.inputFileButton.setObjectName("inputFileButton")
        self.horizontalLayout_9.addWidget(self.inputFileButton)
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_9.addItem(spacerItem9)
        self.gridLayout_7.addLayout(self.horizontalLayout_9, 1, 1, 1, 1)
        spacerItem10 = QtWidgets.QSpacerItem(13, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_7.addItem(spacerItem10, 1, 0, 1, 1)
        spacerItem11 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem11, 7, 1, 1, 1)
        spacerItem12 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem12, 0, 1, 1, 1)
        self.verticalLayout_7.addLayout(self.gridLayout_7)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.verticalLayout_4.addWidget(self.scrollArea)
        self.tabWidget.addTab(self.fileTab, "")
        self.manualTab = QtWidgets.QWidget()
        self.manualTab.setObjectName("manualTab")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.manualTab)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.scrollArea_2 = QtWidgets.QScrollArea(self.manualTab)
        self.scrollArea_2.setWidgetResizable(True)
        self.scrollArea_2.setObjectName("scrollArea_2")
        self.scrollAreaWidgetContents_2 = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_2.setGeometry(QtCore.QRect(0, 0, 323, 327))
        self.scrollAreaWidgetContents_2.setObjectName("scrollAreaWidgetContents_2")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout(self.scrollAreaWidgetContents_2)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.gridLayout_8 = QtWidgets.QGridLayout()
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.label_50 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_50.setFont(font)
        self.label_50.setAlignment(QtCore.Qt.AlignCenter)
        self.label_50.setObjectName("label_50")
        self.gridLayout_8.addWidget(self.label_50, 5, 0, 1, 1)
        self.ft232Inbox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft232Inbox.setFont(font)
        self.ft232Inbox.setObjectName("ft232Inbox")
        self.gridLayout_8.addWidget(self.ft232Inbox, 9, 2, 1, 1)
        self.ft235sInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft235sInBox.setFont(font)
        self.ft235sInBox.setObjectName("ft235sInBox")
        self.gridLayout_8.addWidget(self.ft235sInBox, 8, 4, 1, 1)
        self.line_16 = QtWidgets.QFrame(self.scrollAreaWidgetContents_2)
        self.line_16.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_16.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_16.setObjectName("line_16")
        self.gridLayout_8.addWidget(self.line_16, 1, 0, 1, 5)
        self.u235inBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        self.u235inBox.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.u235inBox.setFont(font)
        self.u235inBox.setObjectName("u235inBox")
        self.gridLayout_8.addWidget(self.u235inBox, 4, 2, 1, 1)
        self.label_44 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_44.setFont(font)
        self.label_44.setAlignment(QtCore.Qt.AlignCenter)
        self.label_44.setObjectName("label_44")
        self.gridLayout_8.addWidget(self.label_44, 10, 0, 1, 1)
        self.label_48 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_48.setFont(font)
        self.label_48.setAlignment(QtCore.Qt.AlignCenter)
        self.label_48.setObjectName("label_48")
        self.gridLayout_8.addWidget(self.label_48, 9, 0, 1, 1)
        self.label_45 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        self.label_45.setMinimumSize(QtCore.QSize(0, 38))
        self.label_45.setAlignment(QtCore.Qt.AlignCenter)
        self.label_45.setObjectName("label_45")
        self.gridLayout_8.addWidget(self.label_45, 0, 4, 1, 1)
        self.label_55 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_55.setFont(font)
        self.label_55.setAlignment(QtCore.Qt.AlignCenter)
        self.label_55.setObjectName("label_55")
        self.gridLayout_8.addWidget(self.label_55, 3, 0, 1, 1)
        self.label_47 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        self.label_47.setAlignment(QtCore.Qt.AlignCenter)
        self.label_47.setObjectName("label_47")
        self.gridLayout_8.addWidget(self.label_47, 0, 2, 1, 1)
        self.ft232sInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft232sInBox.setFont(font)
        self.ft232sInBox.setObjectName("ft232sInBox")
        self.gridLayout_8.addWidget(self.ft232sInBox, 9, 4, 1, 1)
        self.heInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.heInBox.setFont(font)
        self.heInBox.setObjectName("heInBox")
        self.gridLayout_8.addWidget(self.heInBox, 2, 2, 1, 1)
        self.label_53 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_53.setFont(font)
        self.label_53.setAlignment(QtCore.Qt.AlignCenter)
        self.label_53.setObjectName("label_53")
        self.gridLayout_8.addWidget(self.label_53, 4, 0, 1, 1)
        self.ft238InBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft238InBox.setFont(font)
        self.ft238InBox.setObjectName("ft238InBox")
        self.gridLayout_8.addWidget(self.ft238InBox, 7, 2, 1, 1)
        self.label_52 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_52.setFont(font)
        self.label_52.setAlignment(QtCore.Qt.AlignCenter)
        self.label_52.setObjectName("label_52")
        self.gridLayout_8.addWidget(self.label_52, 6, 0, 1, 1)
        self.ft238sInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft238sInBox.setFont(font)
        self.ft238sInBox.setObjectName("ft238sInBox")
        self.gridLayout_8.addWidget(self.ft238sInBox, 7, 4, 1, 1)
        self.label_46 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_46.setFont(font)
        self.label_46.setAlignment(QtCore.Qt.AlignCenter)
        self.label_46.setObjectName("label_46")
        self.gridLayout_8.addWidget(self.label_46, 2, 0, 1, 1)
        self.thInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.thInBox.setFont(font)
        self.thInBox.setObjectName("thInBox")
        self.gridLayout_8.addWidget(self.thInBox, 5, 2, 1, 1)
        self.ft147Inbox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft147Inbox.setFont(font)
        self.ft147Inbox.setObjectName("ft147Inbox")
        self.gridLayout_8.addWidget(self.ft147Inbox, 10, 2, 1, 1)
        self.uInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        self.uInBox.setEnabled(True)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.uInBox.setFont(font)
        self.uInBox.setObjectName("uInBox")
        self.gridLayout_8.addWidget(self.uInBox, 3, 2, 1, 1)
        self.line_18 = QtWidgets.QFrame(self.scrollAreaWidgetContents_2)
        self.line_18.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_18.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_18.setObjectName("line_18")
        self.gridLayout_8.addWidget(self.line_18, 2, 1, 9, 1)
        self.label_54 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        self.label_54.setAlignment(QtCore.Qt.AlignCenter)
        self.label_54.setObjectName("label_54")
        self.gridLayout_8.addWidget(self.label_54, 0, 0, 1, 1)
        self.u235sinBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        self.u235sinBox.setEnabled(False)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.u235sinBox.setFont(font)
        self.u235sinBox.setObjectName("u235sinBox")
        self.gridLayout_8.addWidget(self.u235sinBox, 4, 4, 1, 1)
        self.ft235InBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft235InBox.setFont(font)
        self.ft235InBox.setObjectName("ft235InBox")
        self.gridLayout_8.addWidget(self.ft235InBox, 8, 2, 1, 1)
        self.label_49 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_49.setFont(font)
        self.label_49.setAlignment(QtCore.Qt.AlignCenter)
        self.label_49.setObjectName("label_49")
        self.gridLayout_8.addWidget(self.label_49, 7, 0, 1, 1)
        self.hesInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.hesInBox.setFont(font)
        self.hesInBox.setObjectName("hesInBox")
        self.gridLayout_8.addWidget(self.hesInBox, 2, 4, 1, 1)
        self.smInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.smInBox.setFont(font)
        self.smInBox.setObjectName("smInBox")
        self.gridLayout_8.addWidget(self.smInBox, 6, 2, 1, 1)
        self.ft147sInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.ft147sInBox.setFont(font)
        self.ft147sInBox.setObjectName("ft147sInBox")
        self.gridLayout_8.addWidget(self.ft147sInBox, 10, 4, 1, 1)
        self.line_19 = QtWidgets.QFrame(self.scrollAreaWidgetContents_2)
        self.line_19.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_19.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_19.setObjectName("line_19")
        self.gridLayout_8.addWidget(self.line_19, 0, 3, 1, 1)
        self.line_20 = QtWidgets.QFrame(self.scrollAreaWidgetContents_2)
        self.line_20.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_20.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_20.setObjectName("line_20")
        self.gridLayout_8.addWidget(self.line_20, 2, 3, 9, 1)
        self.thsInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.thsInBox.setFont(font)
        self.thsInBox.setObjectName("thsInBox")
        self.gridLayout_8.addWidget(self.thsInBox, 5, 4, 1, 1)
        self.line_17 = QtWidgets.QFrame(self.scrollAreaWidgetContents_2)
        self.line_17.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_17.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_17.setObjectName("line_17")
        self.gridLayout_8.addWidget(self.line_17, 0, 1, 1, 1)
        self.usInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.usInBox.setFont(font)
        self.usInBox.setObjectName("usInBox")
        self.gridLayout_8.addWidget(self.usInBox, 3, 4, 1, 1)
        self.smsInBox = QtWidgets.QLineEdit(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.smsInBox.setFont(font)
        self.smsInBox.setObjectName("smsInBox")
        self.gridLayout_8.addWidget(self.smsInBox, 6, 4, 1, 1)
        self.label_51 = QtWidgets.QLabel(self.scrollAreaWidgetContents_2)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.label_51.setFont(font)
        self.label_51.setAlignment(QtCore.Qt.AlignCenter)
        self.label_51.setObjectName("label_51")
        self.gridLayout_8.addWidget(self.label_51, 8, 0, 1, 1)
        self.gridLayout_8.setColumnStretch(0, 20)
        self.gridLayout_8.setColumnStretch(1, 1)
        self.gridLayout_8.setColumnStretch(2, 30)
        self.gridLayout_8.setColumnStretch(3, 1)
        self.gridLayout_8.setColumnStretch(4, 40)
        self.verticalLayout_9.addLayout(self.gridLayout_8)
        self.scrollArea_2.setWidget(self.scrollAreaWidgetContents_2)
        self.verticalLayout_8.addWidget(self.scrollArea_2)
        self.tabWidget.addTab(self.manualTab, "")
        self.gridLayout_3.addWidget(self.tabWidget, 1, 0, 1, 1)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem13 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem13)
        self.U235_measuredCheck = QtWidgets.QCheckBox(self.scrollAreaWidgetContents_3)
        self.U235_measuredCheck.setText("")
        self.U235_measuredCheck.setObjectName("U235_measuredCheck")
        self.horizontalLayout_4.addWidget(self.U235_measuredCheck)
        self.label_3 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_4.addWidget(self.label_3)
        spacerItem14 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem14)
        self.gridLayout_3.addLayout(self.horizontalLayout_4, 2, 0, 1, 1)
        self.groupBox = QtWidgets.QGroupBox(self.scrollAreaWidgetContents_3)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.monteCarloCheckBox = QtWidgets.QCheckBox(self.groupBox)
        self.monteCarloCheckBox.setChecked(True)
        self.monteCarloCheckBox.setObjectName("monteCarloCheckBox")
        self.gridLayout.addWidget(self.monteCarloCheckBox, 0, 0, 1, 1)
        self.linearCheckBox = QtWidgets.QCheckBox(self.groupBox)
        self.linearCheckBox.setChecked(True)
        self.linearCheckBox.setObjectName("linearCheckBox")
        self.gridLayout.addWidget(self.linearCheckBox, 1, 0, 1, 1)
        self.gridLayout_3.addWidget(self.groupBox, 3, 0, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem15 = QtWidgets.QSpacerItem(13, 43, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem15)
        self.label_2 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_2.setFont(font)
        self.label_2.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.precisionEnter = QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_3)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.precisionEnter.setFont(font)
        self.precisionEnter.setDecimals(3)
        self.precisionEnter.setSingleStep(0.001)
        self.precisionEnter.setProperty("value", 0.01)
        self.precisionEnter.setObjectName("precisionEnter")
        self.horizontalLayout.addWidget(self.precisionEnter)
        spacerItem16 = QtWidgets.QSpacerItem(13, 43, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem16)
        self.gridLayout_3.addLayout(self.horizontalLayout, 4, 0, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem17 = QtWidgets.QSpacerItem(18, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem17)
        self.decLabel = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        font = QtGui.QFont()
        font.setPointSize(14)
        self.decLabel.setFont(font)
        self.decLabel.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.decLabel.setWordWrap(True)
        self.decLabel.setObjectName("decLabel")
        self.horizontalLayout_2.addWidget(self.decLabel)
        self.decBox = QtWidgets.QSpinBox(self.scrollAreaWidgetContents_3)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.decBox.setFont(font)
        self.decBox.setMaximum(10)
        self.decBox.setProperty("value", 5)
        self.decBox.setObjectName("decBox")
        self.horizontalLayout_2.addWidget(self.decBox)
        spacerItem18 = QtWidgets.QSpacerItem(17, 30, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem18)
        self.gridLayout_3.addLayout(self.horizontalLayout_2, 5, 0, 1, 1)
        self.scrollArea_3.setWidget(self.scrollAreaWidgetContents_3)
        self.gridLayout_2.addWidget(self.scrollArea_3, 0, 0, 3, 1)
        self.gridLayout_2.setColumnStretch(0, 5)
        self.gridLayout_2.setColumnStretch(1, 6)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 836, 22))
        font = QtGui.QFont()
        font.setFamily("MS Shell Dlg 2")
        font.setPointSize(10)
        self.menubar.setFont(font)
        self.menubar.setObjectName("menubar")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionAbout = QtWidgets.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.actionHelp = QtWidgets.QAction(MainWindow)
        self.actionHelp.setObjectName("actionHelp")
        self.actionExit = QtWidgets.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.menuHelp.addAction(self.actionHelp)
        self.menuHelp.addAction(self.actionAbout)
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "CU TRaIL Data Reducer"))
        self.goButton.setText(_translate("MainWindow", "Go!"))
        self.cancelButton.setText(_translate("MainWindow", "Cancel"))
        self.statusMessenger.setText(_translate("MainWindow", "Status"))
        self.label_6.setText(_translate("MainWindow", "Input Type"))
        self.outputFileButton.setText(_translate("MainWindow", "Save Output As"))
        self.pat.setText(_translate("MainWindow", "Parameterize Histograms"))
        self.outputHistCheck.setText(_translate("MainWindow", "Output Histograms"))
        self.inputFileButton.setText(_translate("MainWindow", "Choose Input File"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.fileTab), _translate("MainWindow", "From File"))
        self.label_50.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">232</span>Th</p></body></html>"))
        self.ft232Inbox.setText(_translate("MainWindow", "1.0"))
        self.ft235sInBox.setText(_translate("MainWindow", "0.0"))
        self.u235inBox.setText(_translate("MainWindow", "1.0"))
        self.label_44.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">147</span>F<span style=\" vertical-align:sub;\">T</span></p></body></html>"))
        self.label_48.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">232</span>F<span style=\" vertical-align:sub;\">T</span></p></body></html>"))
        self.label_45.setText(_translate("MainWindow", "Absolute\n"
"Uncertainty (1σ)"))
        self.label_55.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">238</span>U</p></body></html>"))
        self.label_47.setText(_translate("MainWindow", "Value"))
        self.ft232sInBox.setText(_translate("MainWindow", "0.0"))
        self.heInBox.setText(_translate("MainWindow", "1.0"))
        self.label_53.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">235</span>U</p></body></html>"))
        self.ft238InBox.setText(_translate("MainWindow", "1.0"))
        self.label_52.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">147</span>Sm</p></body></html>"))
        self.ft238sInBox.setText(_translate("MainWindow", "0.0"))
        self.label_46.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">4</span>He</p></body></html>"))
        self.thInBox.setText(_translate("MainWindow", "1.0"))
        self.ft147Inbox.setText(_translate("MainWindow", "1.0"))
        self.uInBox.setText(_translate("MainWindow", "1.0"))
        self.label_54.setText(_translate("MainWindow", "Variable"))
        self.u235sinBox.setText(_translate("MainWindow", "0.0"))
        self.ft235InBox.setText(_translate("MainWindow", "1.0"))
        self.label_49.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">238</span>F<span style=\" vertical-align:sub;\">T</span></p></body></html>"))
        self.hesInBox.setText(_translate("MainWindow", "0.0"))
        self.smInBox.setText(_translate("MainWindow", "1.0"))
        self.ft147sInBox.setText(_translate("MainWindow", "0.0"))
        self.thsInBox.setText(_translate("MainWindow", "0.0"))
        self.usInBox.setText(_translate("MainWindow", "0.0"))
        self.smsInBox.setText(_translate("MainWindow", "0.0"))
        self.label_51.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" vertical-align:super;\">235</span>F<span style=\" vertical-align:sub;\">T</span></p></body></html>"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.manualTab), _translate("MainWindow", "Manual"))
        self.label_3.setText(_translate("MainWindow", "Check if <sup>235</sup>U was measured directly"))
        self.groupBox.setTitle(_translate("MainWindow", "Modeling Type"))
        self.monteCarloCheckBox.setText(_translate("MainWindow", "Monte Carlo"))
        self.linearCheckBox.setText(_translate("MainWindow", "Linear Error Propagation"))
        self.label_2.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\">Enter desired precision<br/>of Monte Carlo mean in %</p></body></html>"))
        self.decLabel.setText(_translate("MainWindow", "<html><head/><body><p align=\"center\">Enter number of decimals to report</p></body></html>"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionAbout.setText(_translate("MainWindow", "About"))
        self.actionHelp.setText(_translate("MainWindow", "Help"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))