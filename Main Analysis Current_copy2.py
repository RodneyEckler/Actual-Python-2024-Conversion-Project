"""Data analysis for Grav Lab
Copied structure from APOLLO ATUI code 06/01/12
9/12/19 - fixed run titles from continuously appending - CDH
"""
import matplotlib
matplotlib.use('TkAgg')
from lmfit import minimize, Parameters, Parameter, report_fit


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
#import matplotlib.pyplot as plt

import tkinter as Tk
import tkinter.messagebox as msg
import numpy as num
import scipy
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import sys
from RO import Wdg
import Pmw
from RO import Alg
import csv

import os


# globals
labelFont = ('Times', 10, 'bold')

_HelpPrefix = "GravLab analysis software"


class NoteBook(Pmw.NoteBook):
    def __init__(self, master, **kargs):

            ## set up the notebook
            Pmw.NoteBook.__init__(self, master, **kargs)

            self.FuncCall = Alg.GenericCallback

            self.nb = Pmw.NoteBook(master)
            self.nb.pack(fill='both', expand=1)

            self.p1 = self.nb.add(' Main ')
            self.gr_p1 = Wdg.Gridder(self.p1)

            self.p2 = self.nb.add(' Fit ')
            self.gr_p2 = Wdg.Gridder(self.p2)


            # initialize values
            self.runNum = 00000
            self.nSensors = 0  #number of sensors read in (columns in data file)
            self.nPoints = 0  #number of data poijts taken
            self.nreads = 0  #number of reads in run
            self.sampleInterval = 1.0  #time bewtween reads
            self.overSamplingRate = 10.0  #oversamplig rate in Hz
            self.sensors = []  #list of sensor properties TBD
#            self.dataDir = os.path.join(os.path.expanduser("~"), "Documents", "Gravity Lab", "Python", "Data", "") #a platform and user independent solution
#            self.hdrDir = os.path.join(os.path.expanduser("~"), "Documents", "Gravity Lab", "Python", "Headers", "") 
##            self.dataDir = "C:\\Users\\Jare\\Desktop\\Gravity Lab\\Python\\Data\\"  #data file directory Jeremy
##            self.hdrDir = "C:\\Users\\Jare\\Desktop\\Gravity Lab\\Python\\Headers\\"  #header file directory Jeremy
            #self.dataDir ="C:\\Users\\cdh33\\Documents\\My Dropbox\\LabVIEW Data\\Data\\"         #data file directory office
            #self.hdrDir = "C:\\Users\\cdh33\\Documents\\My Dropbox\\LabVIEW Data\\Headers\\"      #header file directory office
            #self.dataDir="D:\\LabVIEW Data\\data\\" #Acquisition computer
            #self.hdrDir="D:\\LabVIEW Data\\headers\\" #Acquisition computer
            #self.dataDir="L:\\phyx-hoyle-research\LabVIEW Data\\data\\" #Network #TEMPORARY FIX 02.27.2020
            self.dataDir="C:\\LabVIEW Data\\Data\\"
            #self.hdrDir="L:\\phyx-hoyle-research\LabVIEW Data\\headers\\" #Network
            self.hdrDir="C:\\LabVIEW Data\\Headers\\"
            self.chanNamesDict = {}  #names of Data Channels taken from header file. format example: {0:'Sum, etc.}
            self.cal1Dict = {}  #dictionary of linear calibration coefficients
            self.cal2Dict = {}  #dictionary of quadratic calibration coefficients
            self.pageList = []  #list of notebook pages
            self.gridderList = []  #list of gridder objects
            self.pageNameList = []  #list of page names (strings)
            self.buttonList = []  #list of buttons to be enabled/disbled
            self.entryList = []  #list of entrywidgets to be enabled/disbled
            self.tInit = 0.  #initial time for analysis
            self.tFin = 0.
            self.initial = 0
            self.final = 0
            self.runTitle = ''  #run title taken from header
            self.thetaX = []  #diffX/sum
            self.thetaXcal1 = .9945*num.pi/180. #in rad, 1.0747 degrees  #linear calibration coefficient for thetaX
            self.thetaXcal2 = .091*num.pi/180. #0.2070  #quadratic " "
            self.thetaYcal1 = 1.0  #linear calibration coefficient for thetaY
            self.thetaYcal2 = 0.0  #quadratic " "
            self.thetaY = []  #diffY/sum
            self.polyOrder = -1  #polynhomial order (if -1 no polynomial is subtracted)
            self.thetaXcal1Name = 'rad'  #units
            self.thetaXcal2Name = 'rad'  #units
            self.thetaYcal1Name = 'rad'  #units
            self.thetaYcal2Name = 'rad'  #units
            self.polyOrderCut = -1  #polynomial order for cut data
            self.cutLength = 0.0  #length of cts for data
            self.harmList = []
            self.LockinTimeConstant = 0.0
            self.ttor = 0.0 # Torsion Period from Header file
            
            #####
            # Main Page Widgets

            #read data button
            self.readDataButton = Wdg.Button(
                    master=self.p1,
                    text='Read Data',
                    state='disabled',
                    width=20,
                    command=self.readData,
                    helpText='read data and header info'
            )
            self.buttonList.append(self.readDataButton)

            #Run number to be analyzed
            self.runNumWdg = Wdg.IntEntry(self.p1,
                                                                             defValue=None,
                                                                             minValue=00000,
                                                                             maxValue=99999,
                                                                             autoIsCurrent=False,
                                                                             isCurrent=True,
                                                                             helpText='Run # to analyze'
            )
            self.runNumWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.runNumWdg))
            self.runNumWdg.bind('<Return>', self.FuncCall(self.setRunNum, ))

            #data directory
            self.dataDirWdg = Wdg.StrEntry(self.p1,
                                                                              defValue=self.dataDir,
                                                                              autoIsCurrent=False,
                                                                              isCurrent=True,
                                                                              helpText='Data file directory',
                                                                              width=40
            )
            self.dataDirWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.dataDirWdg))
            self.dataDirWdg.bind('<Return>', self.FuncCall(self.setDataDir, ))

            #header directory
            self.hdrDirWdg = Wdg.StrEntry(self.p1,
                                                                             defValue=self.hdrDir,
                                                                             autoIsCurrent=False,
                                                                             isCurrent=True,
                                                                             helpText='Header file directory',
                                                                             width=40
            )
            self.hdrDirWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.hdrDirWdg))
            self.hdrDirWdg.bind('<Return>', self.FuncCall(self.setHdrDir, ))

            # run title
            self.runTitleWdg = Wdg.StrLabel(self.p1,
                                                                               text='',
                                                                               helpText='run title'
            )

            # sampling rate
            self.samplingIntervalWdg = Wdg.FloatLabel(self.p1,
                                                                                                     defValue=None,
                                                                                                     helpText=' sampling rate'
            )

            # oversampling rate
            self.oversamplingRateWdg = Wdg.FloatLabel(self.p1,
                                                                                                     defValue=None,
                                                                                                     helpText=' oversampling rate'
            )

            # oversamples per sample
            self.oversamplesPerSampleWdg = Wdg.FloatLabel(self.p1,
                                                                                                             defValue=None,
                                                                                                             helpText=' oversamples per sample')
            #Laser Frequency (Hz)
            self.LaserFrequencyWdg = Wdg.Label(self.p1,
                                                                                               defValue=None,
                                                                                               helpText=' Laser Frequency')

            #Lockin Time Constant
            self.LockinTimeConstantWdg = Wdg.Label(self.p1,
                                                                                                    defValue=None,
                                                                                                    helpText='Lockin Time Constant')
            
            # oversampling rate
            self.oversamplingRateWdg = Wdg.FloatLabel(self.p1,
                                                                                                     defValue=None,
                                                                                                     helpText=' Oversampling Rate'
            )

            #Sum Lockin Scale
            self.SumLockinScaleWdg = Wdg.Label(self.p1,
                                                                                             defValue=None,
                                                                                             helpText=' Sum Lockin Scale')

            #DiffX Lockin Scale
            self.DiffxLockinScaleWdg = Wdg.Label(self.p1,
                                                                                                     defValue=None,
                                                                                                     helpText=' Diffx Lockin Scale')

            #DiffY Lockin Scale
            self.DiffyLockinScaleWdg = Wdg.Label(self.p1,
                                                                                                     defValue=None,
                                                                                                     helpText=' Diffy Lockin Scale')
            #Torsion Period
            self.TorsionPeriodWdg = Wdg.Label(self.p1,
                                                                                                     defValue=None,
                                                                                                     helpText =' Torsion Period')

            # initial time to analyze
            self.timeInitWdg = Wdg.FloatEntry(self.p1,
                                                                                     defValue=self.tInit,
                                                                                     autoIsCurrent=False,
                                                                                     isCurrent=True,
                                                                                     state='disabled',
                                                                                     helpText=' first time to analyze',
                                                                                     width=35

            )
            self.timeInitWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.timeInitWdg))
            self.timeInitWdg.bind('<Return>', self.FuncCall(self.setTInit))
            self.entryList.append(self.timeInitWdg)

            # final time to analyze
            self.timeFinWdg = Wdg.FloatEntry(self.p1,
                                                                                    autoIsCurrent=False,
                                                                                    isCurrent=True,
                                                                                    defValue=self.tInit,
                                                                                    state='disabled',
                                                                                    helpText=' last time to analyze',
                                                                                    width=35
            )
            self.timeFinWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.timeFinWdg))
            self.timeFinWdg.bind('<Return>', self.FuncCall(self.setTFin))
            self.entryList.append(self.timeFinWdg)

            # Quit Button
            self.quitButton = Wdg.Button(
                    master=self.p1,
                    text='Quit',
                    width=20,
                    command=self._quit,
                    helpText='Exit the program'
            )

            #plot data button
            self.rawPlotButton = Wdg.Button(
                    master=self.p1,
                    text='Plot Data',
                    width=20,
                    command=self.FuncCall(self.rawPlot, 'All'),
                    state='disabled',
                    helpText='plot data'
            )
            self.buttonList.append(self.rawPlotButton)

            #plot fft data button
            self.fftPlotButton = Wdg.Button(
                    master=self.p1,
                    text='FFT Plot',
                    width=20,
                    command=self.FuncCall(self.fftPlot, 'All'),
                    state='disabled',
                    helpText='plot fft data'
            )
            self.buttonList.append(self.fftPlotButton)

            #plot torsion filter button
            self.TorsionFilterPlotButton = Wdg.Button(
                    master=self.p1,
                    text='Torsion Filter Plot',
                    width=20,
                    command=self.FuncCall(self.TorsionFilterPlot, 'All'),
                    state='disabled',
                    helpText='plot Torsion Filter data'
                    )
            self.buttonList.append(self.TorsionFilterPlotButton)

            #torsion filtered fft plot button
            self.TorsionFftPlotButton = Wdg.Button(
                    master=self.p1,
                    text='Torsion fft plot',
                    width=20,
                    command=self.FuncCall(self.TorsionFftPlot, 'All'),
                    state='disabled',
                    helpText='plot Torsion Filtered fft data'
                    )
            self.buttonList.append(self.TorsionFftPlotButton)

            #calibration checkbutton
            self.calOn = Wdg.Checkbutton(
                    master=self.p1,
                    text="Apply Calibration Coefficients",
                    selectcolor='white',
                    defValue=False,
                    state='disabled',
                    helpText='Check to apply calibration coefficients to data'
            )
            self.calOn.bind('<ButtonRelease-1>', self.FuncCall(self.setCal))

            #subtract polynomial fit
            self.polyOrderWdg = Wdg.IntEntry(self.p1,
                                                                                    autoIsCurrent=False,
                                                                                    isCurrent=True,
                                                                                    defValue=-1,
                                                                                    state='disabled',
                                                                                    helpText=' polynomial order',
                                                                                    width=10
            )
            self.polyOrderWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.polyOrderWdg))
            self.polyOrderWdg.bind('<Return>', self.FuncCall(self.setPolyOrder))
            self.entryList.append(self.polyOrderWdg)




            # place widgets
            self.gr_p1.gridWdg('Data File Directory: ', self.dataDirWdg, sticky='w')
            self.gr_p1.gridWdg('Header File Directory: ', self.hdrDirWdg, sticky='w')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg("Run Number: ", self.runNumWdg, sticky='w')
            self.gr_p1.gridWdg('Run Title: ', self.runTitleWdg, sticky='w')
            self.gr_p1.gridWdg("Sample Interval: ", self.samplingIntervalWdg, sticky='w')
            self.gr_p1.gridWdg("Oversampling Rate: ", self.oversamplingRateWdg, sticky='w')
            self.gr_p1.gridWdg("Laser Frequency: ", self.LaserFrequencyWdg,sticky='w')
            self.gr_p1.gridWdg("Lockin Time Constant: ", self.LockinTimeConstantWdg,sticky='w')
            self.gr_p1.gridWdg("Sum Lockin Scale: ",self.SumLockinScaleWdg,sticky='w')
            self.gr_p1.gridWdg("DiffX Lockin Scale: ",self.DiffxLockinScaleWdg,sticky='w')
            self.gr_p1.gridWdg("DiffY Lockin Scale: ",self.DiffyLockinScaleWdg,sticky='w')
            self.gr_p1.gridWdg("Torsion Period: ",self.TorsionPeriodWdg,sticky='w')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg("Initial Time to Analyze [s]: ", self.timeInitWdg, sticky='w')
            self.gr_p1.gridWdg("Final Time to Analyze [s]: ", self.timeFinWdg, 'Enter "0" for last point', sticky='w')
            self.gr_p1.gridWdg("Polynomial Order (-1 for none): ", self.polyOrderWdg, sticky='w')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.calOn, sticky='WNS')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.rawPlotButton, sticky='W')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.fftPlotButton, sticky='W')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.TorsionFilterPlotButton, sticky='W')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.TorsionFftPlotButton, sticky='W')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.readDataButton, row=20, col=2, sticky='W')
            self.gr_p1.gridWdg('', sticky='W')
            self.gr_p1.gridWdg(None, self.quitButton, row=22, col=2, sticky='W')

            ##Fit page widgets
            #
            #set length of cuts in seconds
            self.cutLengthWdg = Wdg.FloatEntry(self.p2,
                                                                                      autoIsCurrent=False,
                                                                                      isCurrent=True,
                                                                                      defValue=0.00,
                                                                                      state='disabled',
                                                                                      helpText=' length of cuts in seconds',
                                                                                      width=10
            )
            self.cutLengthWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.cutLengthWdg))
            self.cutLengthWdg.bind('<Return>', self.FuncCall(self.setCutLength))
            self.entryList.append(self.cutLengthWdg)

            #set polynomial order for cuts
            self.polyOrderCutWdg = Wdg.IntEntry(self.p2,
                                                                                       autoIsCurrent=False,
                                                                                       isCurrent=True,
                                                                                       defValue=1,
                                                                                       state='disabled',
                                                                                       helpText=' polynomial order for cuts',
                                                                                       width=10
            )
            self.polyOrderCutWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.polyOrderCutWdg))
            self.polyOrderCutWdg.bind('<Return>', self.FuncCall(self.setPolyOrderCut))
            self.entryList.append(self.polyOrderCutWdg)

            #cut and fit data
            self.cutPlotButton = Wdg.Button(
                    master=self.p2,
                    text='Cut and Fit',
                    width=20,
                    command=self.FuncCall(self.cutFit),
                    state='disabled',
                    helpText='cut data into sections and fit each section'
            )
            self.buttonList.append(self.cutPlotButton)

            #Torsion filter fit checkbox
            self.TorsionOn = Wdg.Checkbutton(
                    master=self.p2,
                    text="Apply Torsion Filter",
                    selectcolor='white',
                    defValue=False,
                    state='disabled',
                    helpText='Check to apply Torsion Filter to data'
            )
            self.calOn.bind('<ButtonRelease-1>', self.FuncCall(self.setCal))

            #harmonic entry
            self.harmWdg = Wdg.StrEntry(self.p2,
                                                                       defValue='',
                                                                       autoIsCurrent=False,
                                                                       isCurrent=True,
                                                                       helpText='period of harmonics to analyze',
                                                                       width=20
            )
            self.harmWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.harmWdg))
            self.harmWdg.bind('<Return>', self.FuncCall(self.setHarm, ))

             # initial time to analyze
            self.timeInitCutWdg = Wdg.FloatEntry(self.p2,
                                                                                     defValue=self.tInit,
                                                                                     autoIsCurrent=False,
                                                                                     isCurrent=True,
                                                                                     state='disabled',
                                                                                     helpText=' first time to analyze',
                                                                                     width=35

            )
            self.timeInitCutWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.timeInitCutWdg))
            self.timeInitCutWdg.bind('<Return>', self.FuncCall(self.setTInitCut))
            self.entryList.append(self.timeInitCutWdg)

            # final time to analyze
            self.timeFinCutWdg = Wdg.FloatEntry(self.p2,
                                                                                    autoIsCurrent=False,
                                                                                    isCurrent=True,
                                                                                    defValue=self.tInit,
                                                                                    state='disabled',
                                                                                    helpText=' last time to analyze',
                                                                                    width=35
            )
            self.timeFinCutWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=self.timeFinCutWdg))
            self.timeFinCutWdg.bind('<Return>', self.FuncCall(self.setTFinCut))
            self.entryList.append(self.timeFinCutWdg)

            self.gr_p2.gridWdg('Cut Length in seconds: ', self.cutLengthWdg, sticky='w')
            self.gr_p2.gridWdg('Polynomial order for cuts (fixed to 1): ', self.polyOrderCutWdg, sticky='w')
            self.gr_p2.gridWdg('Period list of harmonics to fit (in s, space separated): ', self.harmWdg, sticky='w')
            self.gr_p2.gridWdg('Insert start time: ',self.timeInitCutWdg, sticky='w')
            self.gr_p2.gridWdg('Insert end time: ',self.timeFinCutWdg, sticky='w')
            self.gr_p2.gridWdg(None, self.cutPlotButton, sticky='w')
            self.gr_p2.gridWdg(None, self.TorsionOn, sticky='WNS')

            self.runNumWdg.focus_force()
    

    # Callback Functions

    def readHdr(self):
            filename = self.hdrDir + 'Run' + str(self.runNum) + '.hdr'

            #test to see if file exists
            try:
                    fp = open(filename)
                    self.nohdr = 0
            except:
                    self.nohdr = msg.showwarning("Open file", "Cannot open this file\n(%s)" % filename)
                    return

            hdrReader = csv.reader(open(filename, 'rb'), delimiter=' ')

            next(hdrReader)
            next(hdrReader)
            next(hdrReader)

            self.sampleInterval = float(hdrReader.next()[0])

            next(hdrReader)
            next(hdrReader)
            next(hdrReader)

            a = next(hdrReader)

            self.runTitle = ''

            for i in a:
                    self.runTitle += ' ' + i

            self.runTitleWdg.set(self.runTitle)

            b = next(hdrReader)

            #if len(b) <5:
            #hdrReader.next()

            self.oversampleInterval = float(hdrReader.next()[-1])
            self.oversamplingRateWdg.set(num.divide(1., self.oversampleInterval))

            self.oversamplesPerSample = float(hdrReader.next()[-1])
            self.sampleInterval = self.oversamplesPerSample * self.oversampleInterval
            self.samplingIntervalWdg.set(self.sampleInterval)

            self.LockinTC = float(hdrReader.next()[-1])

            next(hdrReader)
            next(hdrReader)
            next(hdrReader)
            next(hdrReader)

            sens = 0

            for row in hdrReader:
                    if len(row) > 0:
                            row = [''.join(row[:])]
                            row = row[0].strip().split()

                            cal1 = row[-2].replace('(', ' ').replace(')', ' ').split()
                            cal2 = row[-1].replace('(', ' ').replace(')', ' ').split()

                            self.chanNamesDict[sens] = row[0]

                            if len(cal1) > 1:
                                    cal1 = [float(cal1[0]), cal1[1]]
                            else:
                                    cal1 = [float(cal1[0]), 'None']
                            if len(cal2) > 1:
                                    cal2 = [float(cal2[0]), cal2[1]]
                            else:
                                    cal2 = [float(cal2[0]), 'None']

                            self.cal1Dict[sens] = cal1
                            self.cal2Dict[sens] = cal2
                            sens += 1

    def readHdrNew(self):   #if run number greater than equal to 685
            filename = self.hdrDir + 'Run' + str(self.runNum) + '.hdr'

            #test to see if file exists
            try:
                    fp = open(filename)
                    self.nohdr = 0
            except:
                    self.nohdr = msg.showwarning("Open file", "Cannot open this file\n(%s)" % filename)
                    return

            hdrReader = csv.reader(open(filename, 'rb'), delimiter=' ')

            a = next(hdrReader)
            self.runTitle = ''
            for i in a:
                    self.runTitle += ' ' + i

            self.runTitleWdg.set(self.runTitle)

            #b = hdrReader.next() ##needed??
            next(hdrReader) # delete if above line is needed

            self.sampleInterval = float(hdrReader.next()[-1])
            self.samplingIntervalWdg.set(self.sampleInterval) #sets sampling interval to the Wdg

            next(hdrReader) #skips oversamples per sample in header file

            self.oversampleInterval = float(hdrReader.next()[-1])
            self.oversamplingRateWdg.set(self.oversampleInterval)  #oversampling rate

            if int(self.runNum)<726:
                self.LockinTimeConstant = float(hdrReader.next()[-1])
                self.LaserFrequency = float(hdrReader.next()[-1])
                self.SumLockinScale = float(hdrReader.next()[-1])
                self.DiffxLockinScale = float(hdrReader.next()[-1])
                self.DiffyLockinScale = float(hdrReader.next()[-1])
                next(hdrReader)
            elif int(self.runNum)>735:
                self.LockinTimeConstant = ' '.join(hdrReader.next()[3:5])
                self.LaserFrequency = ' '.join(hdrReader.next() [2:5])
                self.SumLockinScale=' '.join(hdrReader.next() [3:5])
                self.DiffxLockinScale=' '.join(hdrReader.next() [3:5])
                self.DiffyLockinScale=' '.join(hdrReader.next() [3:5])
                self.ttor = float(hdrReader.next()[-1])
                self.thetaCal = float(hdrReader.next()[-1])
            else:
                self.LockinTimeConstant = ' '.join(hdrReader.next()[3:5])
                self.LaserFrequency = ' '.join(hdrReader.next() [2:5])
                self.SumLockinScale=' '.join(hdrReader.next() [3:5])
                self.DiffxLockinScale=' '.join(hdrReader.next() [3:5])
                self.DiffyLockinScale=' '.join(hdrReader.next() [3:5])
                next(hdrReader)
            self.LockinTimeConstantWdg.set(self.LockinTimeConstant)
            self.LaserFrequencyWdg.set(self.LaserFrequency)
            self.SumLockinScaleWdg.set(self.SumLockinScale)
            self.DiffxLockinScaleWdg.set(self.DiffxLockinScale)
            self.DiffyLockinScaleWdg.set(self.DiffyLockinScale)
            self.TorsionPeriodWdg.set(self.ttor)

            self.ttor = self.ttor/self.sampleInterval

            next(hdrReader)
            next(hdrReader)
            next(hdrReader)

            sens = 0

            for row in hdrReader:
                    if len(row) > 0:
                            row = [''.join(row[:])]
                            row = row[0].strip().split()

                            cal1 = row[-2].replace('(', ' ').replace(')', ' ').split()
                            cal2 = row[-1].replace('(', ' ').replace(')', ' ').split()

                            self.chanNamesDict[sens] = row[0]

                            if len(cal1) > 1:
                                    cal1 = [float(cal1[0]), cal1[1]]
                            else:
                                    cal1 = [float(cal1[0]), 'None']
                            if len(cal2) > 1:
                                    cal2 = [float(cal2[0]), cal2[1]]
                            else:
                                    cal2 = [float(cal2[0]), 'None']

                            self.cal1Dict[sens] = cal1
                            self.cal2Dict[sens] = cal2
                            sens += 1

    def readData(self):

            #read header file
            if int(self.runNum)<685:
                    self.readHdr()            #reads data from readHdr
            if int(self.runNum)>=685:
                    self.readHdrNew()         #reads data from readHdrNew

            #raise error and exit if no header found
            if self.nohdr:
                    msg.showerror("Error", "No Header File Found")
                    return

            filename = self.dataDir + 'Run' + str(self.runNum) + '.dat'

            #test to see if file exists
            try:
                    fp = open(filename)
            except:
                    msg.showwarning("Error", "Cannot open this file\n(%s)" % filename)
                    return

            self.timeInitWdg.set(0.0)
            self.timeFinWdg.set(0.0)
            self.polyOrder = -1
            self.polyOrderWdg.set(self.polyOrder)

            datReader = csv.reader(open(filename, 'rb'), delimiter='\t')
            firstRow = next(datReader)
            self.nSensors = len(firstRow)  # number of sensors read in
            self.sensDat = []
            floatRow0 = []

            cal1 = []
            cal2 = []

            if self.calOn.getBool(): 
                    for i in range(self.nSensors): #takes cal coefficient from each sens
                            cal1.append(list(self.cal1Dict.values())[i][0])
                            cal2.append(list(self.cal2Dict.values())[i][0])
            else:
                    cal1 = num.ones(self.nSensors, float)
                    cal2 = num.zeros(self.nSensors, float)

            for i in range(self.nSensors):
                    floatRow0.append(cal1[i] * float(firstRow[i]) + cal2[i] * float(firstRow[i]) ** 2)

            self.sensDat.append(floatRow0)

            for row in datReader:
                    floatRow = []

                    for i in range(self.nSensors):
                            val = float(row[i])
                            floatRow.append(cal1[i] * val + cal2[i] * val ** 2)

                    self.sensDat.append(floatRow)

            self.datArray = num.array(self.sensDat)
            self.datArray = self.datArray.transpose()
            
        
            # number of data points read in
            self.nreads = len(self.datArray[0])

            # make array of time values using sample interval
            self.timeArray = num.arange(0., self.nreads * self.sampleInterval, self.sampleInterval)

            # useful indices to have for later
            for i in range(len(self.chanNamesDict)):
                    if self.chanNamesDict[i] == 'Sum':
                            self.sumIndex = list(self.chanNamesDict.keys())[i]
                    if self.chanNamesDict[i] == 'DiffX':
                            self.diffXIndex = list(self.chanNamesDict.keys())[i]
                    if self.chanNamesDict[i] == 'DiffY':
                        self.diffYIndex = list(self.chanNamesDict.keys())[i]

            # calculate thetaX and append to data array
            if 'DiffX' and 'Sum' in list(self.chanNamesDict.values()):
                    DiffOverSumX = num.divide(self.datArray[self.diffXIndex], self.datArray[self.sumIndex])
                    if int(self.runNum)<685:
                        c1=1
                        c2=1
                    
                    if self.calOn.getBool():
                        if int(self.runNum)>=685:
                            
                            dscale = float(self.DiffxLockinScale.split(' ')[0])
                            sscale = float(self.SumLockinScale.split(' ')[0])
                            
                            c1 = self.thetaXcal1*(dscale/100)/(sscale/100)      #100 and 100 refer to the scale used for the specific calibration
                            c2 = self.thetaXcal2*(dscale/100)/(sscale/100)*(dscale/100)/(sscale/100)  

                        self.thetaX = c1 * DiffOverSumX + c2 * DiffOverSumX * DiffOverSumX
                    else:
                            self.thetaX = DiffOverSumX

                    self.datArray = num.vstack((self.datArray, self.thetaX))
                    self.chanNamesDict.update({self.nSensors: 'ThetaX'})
                    self.cal1Dict.update({self.nSensors: [self.thetaXcal1, self.thetaXcal1Name]})
                    self.cal2Dict.update({self.nSensors: [self.thetaXcal2, self.thetaXcal2Name]})
                    self.nSensors += 1

            # calculate thetaY and append to data array
            if 'DiffY' and 'Sum' in list(self.chanNamesDict.values()):

                    DiffOverSumY = num.divide(self.datArray[self.diffYIndex], self.datArray[self.sumIndex])
                    c1 = self.thetaYcal1
                    c2 = self.thetaYcal2

                    if self.calOn.getBool():
                            self.thetaY = c1 * DiffOverSumY + c2 * DiffOverSumY * DiffOverSumY
                    else:
                            self.thetaY = DiffOverSumY

                    self.datArray = num.vstack((self.datArray, self.thetaY))
                    self.chanNamesDict.update({self.nSensors: 'ThetaY'})
                    self.cal1Dict.update({self.nSensors: [self.thetaYcal1, self.thetaYcal1Name]})
                    self.cal2Dict.update({self.nSensors: [self.thetaYcal2, self.thetaYcal2Name]})
                    self.nSensors += 1

            # these are useful indices to have for later
            for i in range(len(self.chanNamesDict)):
                    if self.chanNamesDict[i] == 'ThetaX':
                            self.thetaXIndex = list(self.chanNamesDict.keys())[i]
                    if self.chanNamesDict[i] == 'ThetaY':
                            self.thetaYIndex = list(self.chanNamesDict.keys())[i]

            self.tInit = self.timeArray[0]
            self.tFin = self.timeArray[-1]

            self.initial = 0
            self.final = self.nreads ## - 1 

            for wdg in self.buttonList:
                    wdg.configure(state='active')
            for wdg in self.entryList:
                    wdg.configure(state='normal')
            self.polyOrderCutWdg.configure(state='disabled')
                 
            self.calOn.configure(state='disabled')
            self.readDataButton.configure(state='disabled')
            self.polyOrderWdg.configure(state='normal')

            self.rawPlot('All')

            
            
    def rawPlot(self, sensList):

            initial = int(self.initial)
            final = int(self.final)

            self.delPlots()

            if 'DiffY' and 'Sum' in list(self.chanNamesDict.values()):

                    ## make scatter plot of XY/Sum
                    xy = self.nb.add(' XY ')
                    self.pageList.append(xy)
                    self.gridderList.append(Wdg.Gridder(xy))

                    xyf = Figure(figsize=(4, 4), dpi=100)
                    plt = xyf.add_subplot(111)
                    
                    if self.runNum == 8300: ## delete, keep else portion.
                        plt.plot(self.thetaX[initial:initial+10000], self.thetaY[initial:initial+10000]) ## 2/25 get rid of this if /else keep the else
                    else:
                        plt.plot(self.thetaX[initial:final], self.thetaY[initial:final])
                        
                    if self.calOn.getBool():
                            plt.set_title('Scatter Plot')
                    else:
                            plt.set_title(' Raw Scatter Plot')

                    plt.set_xlabel('Theta X')
                    plt.set_ylabel('Theta Y')

                    canvas = FigureCanvasTkAgg(xyf, xy)
                    try:
                        canvas.show()
                    except OverflowError: # if data is too large
                        print( "Data are too large to plot.")
                    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                    toolbar = NavigationToolbar2TkAgg(canvas, xy)
                    toolbar.update()
                    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                                
            if sensList == 'All':
                    for i in range(self.nSensors):
                            a = self.nb.add(self.chanNamesDict[i])
                            self.pageList.append(a)
                            self.gridderList.append(Wdg.Gridder(a))

                            f = Figure(figsize=(4, 4), dpi=100)
                            plt = f.add_subplot(111)
                           
                            if self.runNum == 8300:
                                x= self.timeArray[initial:initial+10000] ## 2/25, above comment
                                y= self.datArray[i][initial:initial+10000]
                            else:
                                x = self.timeArray[initial:final]
                                y = self.datArray[i][initial:final]
                                
                            plt.plot(x, y)

                            if self.calOn.getBool():
                                    plt.set_title('Calibrated Data')
                                    plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + ']')
                            else:
                                    plt.set_title('Raw Data')
                                    plt.set_ylabel(self.chanNamesDict[i] + ' [Dimensionless]')

                            plt.set_xlabel('Time [s]')
                            plt.set_xlim(min(x), max(x))
                            plt.set_ylim(min(y) - 0.05 * abs(min(y)), max(y) + 0.05 * abs(max(y)))

                            canvas = FigureCanvasTkAgg(f, a)
                            canvas.show()
                            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                            toolbar = NavigationToolbar2TkAgg(canvas, a)
                            toolbar.update()
                            canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
                       
                    print(('thetay',self.datArray[-1]))
                    print(('thetax',self.datArray[-2]))
                    print(('?',self.datArray[-3]))
    def TorsionFilterPlot(self, sensList):
        initial = int(self.initial)
        final = int(self.final)
        self.delPlots()


        torFilter = []
        TorsionList=[self.sumIndex,self.diffXIndex,self.diffYIndex,self.thetaXIndex,self.thetaYIndex]
       
        
        
        
        torStart = int(num.floor(self.ttor*0.25) + 1)
        torEnd = int(num.floor(len(self.datArray[0])-(self.ttor*0.25)-1))
        print(('torsionlist',TorsionList))
        print(('torStart',torStart))
        print(('torEnd',torEnd))
        print(('ln of dict',len(self.chanNamesDict)))
        for i in range(len(self.chanNamesDict)):
           tor = []
           if i in TorsionList:
##               corFactor = 1/(num.cos((num.pi * self.ttor)/(2 * 201.8/self.sampleInterval))) #have python find attractor period
               """corFactor = 1/(num.cos((num.pi * self.ttor)/(2 * 100.0/self.sampleInterval))) #I forgot what the 100 is. Based on the above, I assume it is meant to be the attractor period.
               corfactor should appear in the fitting, not in the plot"""
               
               for j in range(torStart, torEnd):
                   qb0 = j - self.ttor*0.25 # quarter of period before,float
                   qa0 = j + self.ttor*0.25 # quarter of period after, float
                   qb1 = int(num.floor(qb0) + 1)  # discrete quarter of period before + 1, int
                   qa1 = int(num.floor(qa0) + 1)  # discrete quarter of period after + 1, int
                   #This works by finding the slope of the line connecting datArray[i][qb0] and datArray[i][qb1]
                   #and interpolating data between them
                   interp0 = (self.datArray[i][qb1] - self.datArray[i][int(num.floor(qb0))]) * (qb0 - num.floor(qb0)) + \
                                                self.datArray[i][int(num.floor(qb0))] #interpolated point quarter before, assumes linear
                   interp1 = (self.datArray[i][qa1] - self.datArray[i][int(num.floor(qa0))]) * (qa0 - num.floor(qa0)) + \
                                                self.datArray[i][int(num.floor(qa0))] #interpolated point quarter after, assumes linear
                   val = (interp0 + interp1)/2
#                   val = val* corFactor
                   tor.append(val)
               torFilter.append(tor)
           else:
               for j in range(int(torStart), int(torEnd)):
                    val = self.datArray[i][j]
                    tor.append(val)
               torFilter.append(tor)
        '''       
        print('qb0',qb0)
        print('qa0',qa0)
        print('qb1',qb1)
        print('qa1',qa1)
        print('interp0',interp0)
        print('interp1',interp1)
        print('data[i]',self.datArray[0])
        '''
        print(('torfilter',torFilter))
        #print("self.nsensors",self.nSensors)
        #print("chanNamesDic",self.chanNamesDict[i])
        
        torFilter = num.array(torFilter)

        if 'DiffY' and 'Sum' in list(self.chanNamesDict.values()):

            ## make scatter plot of XY/Sum
            xy = self.nb.add(' XY ')
            self.pageList.append(xy)
            self.gridderList.append(Wdg.Gridder(xy))

            xyf = Figure(figsize=(4,4), dpi=100)
            plt = xyf.add_subplot(111)
            plt.plot(self.thetaX[initial:final], self.thetaY[initial:final])
            #print self.thetaY.shape

            if self.calOn.getBool():
                plt.set_title('Scatter Plot')
            else:
                plt.set_title(' Torsion filter Scatter Plot')

            plt.set_xlabel('Theta X')
            plt.set_ylabel('Theta Y')

            canvas = FigureCanvasTkAgg(xyf, xy)
            canvas.show()
            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, xy)
            toolbar.update()
            canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
            # make diff & sum plot
            '''ds=self.nb.add(' DiffSum ')
                        self.pageList.append(ds)
                        self.gridderList.append(Wdg.Gridder(ds))

                        dsf = Figure(figsize=(4,4), dpi=100)
                        plt = dsf.add_subplot(111)
                        plt.plot(self.timeArray[initial:final],arrn[self.diffXIndex][initial:final])
                        plt.plot(self.timeArray[initial:final],arrn[self.sumIndex][initial:final])

                        canvas = FigureCanvasTkAgg(dsf, ds)
                        canvas.show()
                        canvas.get_tk_widged().pack(side=Tk.TOP,fill=Tk.BOTH, expand=1)

                        toolbar = NavigationToolbar2TkAgg(canvas,ds)
                        toolbar.update()
                        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1) '''

        if sensList == 'All':
            for i in range(self.nSensors):
                a = self.nb.add(self.chanNamesDict[i])
                self.pageList.append(a)
                self.gridderList.append(Wdg.Gridder(a))

                f = Figure(figsize=(4, 4), dpi=100)
                plt = f.add_subplot(111)

                # plots both filtered data and raw data,
                # only plots filtered data on data in TorsionList
              

                x = self.timeArray[initial:final]
                y = self.datArray[i][initial:final]
                plt.plot(x, y, 'b-')
                
                if i in TorsionList:
                    xTorsion = self.timeArray[torStart:torEnd]
                    yTorsion = torFilter[i]
                    plt.plot(xTorsion, yTorsion, 'g-')
                else:
                    pass

                if self.calOn.getBool():
                    plt.set_title('Calibrated Data')
                    plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + ']')
                else:
                    plt.set_title('Torsion Filtered Data')
                    plt.set_ylabel(self.chanNamesDict[i] + ' [Dimensionless]')

                plt.set_xlabel('Time [s]')
                plt.set_xlim(min(x), max(x))
                plt.set_ylim(min(y) - 0.05 * abs(min(y)), max(y) + 0.05 * abs(max(y)))

                canvas = FigureCanvasTkAgg(f, a)
                canvas.show()
                canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                toolbar = NavigationToolbar2TkAgg(canvas, a)
                toolbar.update()
                canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)


    def TorsionFftPlot(self,sensList):

        initial = int(self.initial)
        final = int(self.final)

        self.delPlots()

        plotDir = self.dataDir[0:-6] + "\\fft\\"
        
        torStart = int(num.floor(self.ttor*0.25) + 1)
        torEnd = int(num.floor(len(self.datArray[0])-(self.ttor*0.25)-1))
        torFilter = []
        TorsionList=[self.sumIndex,self.diffXIndex,self.diffYIndex,self.thetaXIndex,self.thetaYIndex]
        for i in range(len(self.chanNamesDict)):
           tor = []
           if i in TorsionList:
##               corFactor = 1/(num.cos((num.pi * self.ttor)/(2 * (float(self.harmList[0])/self.sampleInterval))))
 #              corFactor = 1/(num.cos((num.pi * self.ttor)/(2 * 201.8/self.sampleInterval))) #better determine attractor mass period, place instead of 201.8
               for j in range(torStart, torEnd):
                   qb0 = j - self.ttor*0.25 # quarter of period before,float
                   qa0 = j + self.ttor*0.25 # quarter of period after, float
                   qb1 = int(num.floor(qb0) + 1)  # discrete quarter of period before + 1, int
                   qa1 = int(num.floor(qa0) + 1) # discrete quarter of period after + 1, int
                   interp0 = (self.datArray[i][qb1] - self.datArray[i][int(num.floor(qb0))])*(qb0 - num.floor(qb0)) + \
                                                self.datArray[i][int(num.floor(qb0))] #interpolated point quarter before, assumes linear
                   interp1 = (self.datArray[i][qa1] - self.datArray[i][int(num.floor(qa0))])*(qa0 - num.floor(qa0)) + \
                                                self.datArray[i][int(num.floor(qa0))] #interpolated point quarter after, assumes linear
                   val = (interp0 + interp1)/2
#                   val = val * corFactor
                   tor.append(val)
               torFilter.append(tor)
           else:
               for j in range(torStart, torEnd):
                    val = self.datArray[i][j]
                    tor.append(val)
               torFilter.append(tor)
               
        torFilter = num.array(torFilter)
        
        freqLength = len(num.fft.rfft(torFilter[0][initial:final]))
        freqMin = num.divide(1.0, freqLength * self.sampleInterval)
        freqMax = num.divide(0.5, self.sampleInterval)
        freqStep = num.divide(freqMax, float(freqLength))

        self.freqArray = num.arange(0, freqMax, freqStep)

        if len(self.freqArray) > freqLength:
                self.freqArray = self.freqArray[0:-1]

        self.delPlots()

        if sensList == 'All':
                for i in range(self.nSensors):
                        a = self.nb.add(self.chanNamesDict[i])
                        self.pageList.append(a)
                        self.gridderList.append(Wdg.Gridder(a))

                        f = Figure(figsize=(4, 4), dpi=100)
                        plt = f.add_subplot(111)
                        x = self.freqArray
                        y = num.divide(1., freqLength) * num.abs(num.fft.rfft(torFilter[i][initial:final]))
                        plt.loglog(x, y)

                        if self.calOn.getBool():
                                plt.set_title('Amplitude Spectrum (Calibrated)')
                                plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + '/Sqrt(Hz)]')
                        else:
                                plt.set_title('Amplitude Spectrum (Raw)')
                                plt.set_ylabel(self.chanNamesDict[i] + ' [V/Sqrt(Hz)]')

                        plt.set_xlabel('Frequency [Hz]')

                        plt.set_xlim(min(x), max(x))
                        plt.set_ylim(0.8 * min(y[1:]), 1.2 * max(y[1:]))

                        canvas = FigureCanvasTkAgg(f, a)
                        canvas.show()
                        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                        toolbar = NavigationToolbar2TkAgg(canvas, a)
                        toolbar.update()
                        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                        fname = plotDir+'Run'+str(self.runNum)+self.chanNamesDict[i]+'.fft'
                        plotPoints = num.append(x,y).reshape(2,len(x)).transpose()
                        num.savetxt(fname, plotPoints)


    def fftPlot(self, sensList):

            initial = int(self.initial)
            final = int(self.final)

            freqLength = 0

            freqLength = len(num.fft.rfft(self.datArray[0][initial:final]))
            freqMin = num.divide(1.0, freqLength * self.sampleInterval)
            freqMax = num.divide(0.5, self.sampleInterval)
            freqStep = num.divide(freqMax, float(freqLength))

            self.freqArray = num.arange(0, freqMax, freqStep)

            if len(self.freqArray) > freqLength:
                    self.freqArray = self.freqArray[0:-1]

            self.delPlots()

            plotDir = self.dataDir[0:-6] + "\\fft\\"

            if sensList == 'All':
                    for i in range(self.nSensors):
                            a = self.nb.add(self.chanNamesDict[i])
                            self.pageList.append(a)
                            self.gridderList.append(Wdg.Gridder(a))

                            f = Figure(figsize=(4, 4), dpi=100)
                            plt = f.add_subplot(111)
                            x = self.freqArray
                            y = num.divide(1., freqLength) * num.abs(num.fft.rfft(self.datArray[i][initial:final]))
                            plt.loglog(x, y)

                            if self.calOn.getBool():
                                    plt.set_title('Amplitude Spectrum (Calibrated)')
                                    plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + '/Sqrt(Hz)]')
                            else:
                                    plt.set_title('Amplitude Spectrum (Raw)')
                                    plt.set_ylabel(self.chanNamesDict[i] + ' [V/Sqrt(Hz)]')

                            plt.set_xlabel('Frequency [Hz]')

                            plt.set_xlim(min(x), max(x))
                            plt.set_ylim(0.8 * min(y[1:]), 1.2 * max(y[1:]))

                            canvas = FigureCanvasTkAgg(f, a)
                            canvas.show()
                            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                            toolbar = NavigationToolbar2TkAgg(canvas, a)
                            toolbar.update()
                            canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                            fname = plotDir+'Run'+str(self.runNum)+self.chanNamesDict[i]+'.fft'
                            plotPoints = num.append(x,y).reshape(2,len(x)).transpose()
                            num.savetxt(fname, plotPoints)

    def delPlots(self):
            if len(self.nb.pagenames()) > 2:
                    for page in self.nb.pagenames()[2:]:
                            self.nb.delete(page)

    def setTInit(self, evt):
            """set the initial time to analyze"""
            val = float(self.timeInitWdg.get())

            if val > float(self.tFin) and float(self.tFin) > 0:
                    msg.showwarning("Error", "Invalid Time Interval")
                    return
            else:
                    self.tInit = val
                    self.initial = num.floor(num.divide(float(self.tInit), self.sampleInterval))
                    self.timeInitWdg.setIsCurrent(True)

    def setTInitCut(self, evt):
            """set the initial time to analyze"""
            val = float(self.timeInitCutWdg.get())

            if val > float(self.tFin) and float(self.tFin) > 0:
                    msg.showwarning("Error", "Invalid Time Interval")
                    return
            else:
                    self.tInit = val
                    self.initial = num.floor(num.divide(float(self.tInit), self.sampleInterval))
                    self.timeInitCutWdg.setIsCurrent(True)

    def setTFin(self, evt):
            """set the final time to analyze"""
            val = float(self.timeFinWdg.get())

            if val < float(self.tInit) and val != 0.0:
                    msg.showwarning("Error", "Invalid Time Interval")
                    return
            elif val == 0:
                    self.tFin = self.timeArray[-1]
                    self.final = self.nreads
                    self.timeFinWdg.setIsCurrent(True)
            else:
                    self.tFin = val
                    self.final = num.floor(num.divide(float(self.tFin), self.sampleInterval))
                    self.timeFinWdg.setIsCurrent(True)

    def setTFinCut(self, evt):
            """set the final time to analyze"""
            val = float(self.timeFinCutWdg.get())

            if val < float(self.tInit) and val != 0.0:
                    msg.showwarning("Error", "Invalid Time Interval")
                    return
            elif val == 0:
                    self.tFin = self.timeArray[-1]
                    self.final = self.nreads
                    self.timeFinCutWdg.setIsCurrent(True)
            else:
                    self.tFin = val
                    self.final = num.floor(num.divide(float(self.tFin), self.sampleInterval))
                    self.timeFinCutWdg.setIsCurrent(True)

    def setCal(self, evt):
            """choose to apply calibrations or not"""
            state = self.calOn.getBool()

            self.runNumWdg.setIsCurrent(True)

    def setRunNum(self, evt):
            """set the run number to be analyzed"""
            self.runNum = self.runNumWdg.get()
            self.calOn.configure(state='active')
            self.readDataButton.configure(state='active')
            self.rawPlotButton.configure(state='disabled')
            self.fftPlotButton.configure(state='disabled')
            self.TorsionFilterPlotButton.configure(state='disabled')
            self.TorsionFftPlotButton.configure(state='disabled')
            self.TorsionOn.configure(state='active')
            self.runNumWdg.setIsCurrent(True)
            self.readDataButton.focus_force()

    def setDataDir(self, evt):
            """set the run number to be analyzed"""
            self.dataDir = self.dataDirWdg.get()
            self.dataDirWdg.setIsCurrent(True)

    def setHdrDir(self, evt):
            """set the run number to be analyzed"""
            self.hdrDir = self.hdrDirWdg.get()
            self.hdrDirWdg.setIsCurrent(True)

    def setPolyOrder(self, evt):
            """set the polynomial order"""
            """subtracts polynomial from data. Does not subtract a0 (offset) term """

            initial = int(self.initial)
            final = int(self.final)
            y = []
            self.polyOrder = int(self.polyOrderWdg.get())
            self.polyOrderWdg.setIsCurrent(True)
            x = self.timeArray[initial:final]

            self.delPlots()

            self.polyOrderWdg.configure(state='disabled')

            for i in range(self.nSensors):
                    y.append(self.datArray[i][initial:final])
                    fit = num.polyfit(x, y[i], self.polyOrder)
                    term = num.zeros(len(y[i]))

                    for j in range(len(fit)):
                            porder = len(fit) - j - 1
                            term += -fit[j] * x ** (len(fit) - j - 1)
                            print((self.chanNamesDict[i] + '  polycoeff  ' + str(porder) + ': ' + str(fit[j])))
                    print (' ')

                    for k in range(len(y[i])):
                            self.datArray[i][k] = y[i][k] + term[k]
                            
                    self.datArray[i]+= fit[-1]

                    self.initial = 0
                    self.final = len(self.timeArray[initial:final])

                    a = self.nb.add(self.chanNamesDict[i])
                    self.pageList.append(a)
                    self.gridderList.append(Wdg.Gridder(a))

                    f = Figure(figsize=(4, 4), dpi=100)
                    plt = f.add_subplot(111)
                    
                    if self.calOn.getBool():
                            plt.set_title('Polynomial Subtracted (Calibrated)')
                            plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + ']')
                    else:
                            plt.set_title('Polynomial Subtracted (Raw)')
                            plt.set_ylabel(self.chanNamesDict[i] + ' [Dimensionless]')

                    plt.set_xlabel('Time [s]')
                    x1 = self.timeArray[initial:final]
                    y1 = self.datArray[i][self.initial:self.final]
                    plt.plot(x1, y1)

                    plt.set_xlim(min(x1), max(x1))
                   
                    if min(y1[1:])<0:
                        ylow = 1.1*min(y1[1:])
                    else: 
                        ylow = 0.9 * min(y1[1:])
                        
                    if max(y1[1:])<0:
                        yhi = 0.9*max(y1[1:])
                    else: 
                        yhi = 1.1*max(y1[1:])
                        
                    plt.set_ylim(ylow, yhi)

                    canvas = FigureCanvasTkAgg(f, a)
                    canvas.show()
                    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

                    toolbar = NavigationToolbar2TkAgg(canvas, a)
                    toolbar.update()
                    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    def setPolyOrderCut(self, evt):
            self.polyOrderCut = int(self.polyOrderCutWdg.get())
            self.polyOrderCutWdg.setIsCurrent(True)

    def setCutLength(self, evt):
            self.cutLength = float(self.cutLengthWdg.get())
            self.cutLengthWdg.setIsCurrent(True)

    def cutFit(self):
        torStart = int(num.floor(self.ttor*0.25) + 1)
        torEnd = int(num.floor(len(self.datArray[0])-(self.ttor*0.25)-1))                       
        torFilter = []
        TorsionList=[self.sumIndex,self.diffXIndex,self.diffYIndex,self.thetaXIndex,self.thetaYIndex]
        for i in range(len(self.chanNamesDict)):
           tor = []
           if i in TorsionList:
               corFactor = 1/(num.cos((num.pi * self.ttor)/(2 * (float(self.harmList[0])/self.sampleInterval))))
               for j in range(torStart, torEnd):
                   qb0 = j - self.ttor*0.25 # quarter of period before,float
                   qa0 = j + self.ttor*0.25 # quarter of period after, float
                   qb1 = int(num.floor(qb0) + 1) # discrete quarter of period before + 1, int
                   qa1 = int(num.floor(qa0) + 1) # discrete quarter of period after + 1, int
                   #This works by finding the slope of the line connecting datArray[i][qb0] and datArray[i][qb1]
                   #and interpolating data between them
                   interp0 = (self.datArray[i][qb1] - self.datArray[i][int(num.floor(qb0))])*(qb0 - num.floor(qb0)) + \
                                                self.datArray[i][int(num.floor(qb0))] #interpolated point quarter before, assumes linear
                   interp1 = (self.datArray[i][qa1] - self.datArray[i][int(num.floor(qa0))])*(qa0 - num.floor(qa0)) + \
                                                self.datArray[i][int(num.floor(qa0))] #interpolated point quarter after, assumes linear
                   val = (interp0 + interp1)/2
                   val = val * corFactor
                   tor.append(val)
               torFilter.append(tor)
           else:
               for j in range(torStart, torEnd):
                    val = self.datArray[i][j]
                    tor.append(val)
               torFilter.append(tor)     
        torFilter = num.array(torFilter)
        torLength = len(torFilter)

        
        self.cutPoints = num.int(num.floor(num.divide(self.cutLength, self.sampleInterval)))
        cutPoints = self.cutPoints
        initial=int(self.initial)
        final=int(self.final)
        self.newdat=[]
        for i in range(self.nSensors):
            if self.TorsionOn.getBool():
                self.newdat.append(torFilter[i][initial:final-cutPoints])
            else:
                self.newdat.append(self.datArray[i][initial:final])
    
    
        self.ncuts = int(num.floor(num.divide(len(self.newdat[0]), cutPoints)))
        self.polyOrderCut = 1                                                       #needs to change

        self.delPlots()

        self.cutAttempt = 0
        self.removeList = []

        self.senssub = num.zeros(shape=(self.nSensors, self.ncuts, cutPoints))

        self.fitList = num.zeros(shape=(self.nSensors, self.ncuts, 2))               # 2 needs to change
        self.errList = num.zeros(shape=(self.nSensors, self.ncuts, 2))               # 2 needs to change

        self.fitWdgList = []
        self.fitFigList = []
        self.fitCanList = []
        
        self.cutLevels = num.zeros(self.nSensors)

        self.aList = num.zeros(shape=(self.nSensors, self.ncuts))
        self.aErrList = num.zeros(shape=(self.nSensors, self.ncuts))
        self.bList = num.zeros(shape=(self.nSensors, self.ncuts))
        self.bErrList = num.zeros(shape=(self.nSensors, self.ncuts))

        # (More detailed definitions can be found throughout the function)
        # List of a and b amplitudes for each cut
        self.fits = num.zeros(shape=(self.nSensors, self.ncuts, 2, len(self.harmList))).tolist()
        # Averaged a and b amplitudes
        self.avgs = [[0 for m in range(len(self.harmList))] for n in range(2)]
        # Standard deviations of a and b
        self.stds = [[0 for m in range(len(self.harmList))] for n in range(2)]
        # Average overall amplitude (0) and error (1)
        self.ampavg = [[0 for m in range(len(self.harmList))] for n in range(2)]
        
        # Overall amplitude per cut
        self.avgcut = num.zeros(shape=(self.ncuts, len(self.harmList))).tolist()
        # Errors of a and b amplitudes for each cut and harmonic (for writestats function)
        self.errs = num.zeros(shape=(self.nSensors, self.ncuts, 2, len(self.harmList))).tolist()
        # polyfit coefficients
        self.polyfits = num.zeros(shape=(self.nSensors, self.ncuts, self.polyOrderCut + 1))
        self.polyfits2 = num.zeros(shape=(self.nSensors, self.ncuts, self.polyOrderCut + 1))
        # function fitted per cut, fittedFunc2 is from writestats()
        self.fittedFunc = num.zeros(shape=(self.nSensors, self.ncuts, cutPoints))
        self.fittedFunc2 = num.zeros(shape=(self.nSensors, self.ncuts, cutPoints))
        # residuals
        self.residuals = num.zeros(shape=(self.nSensors, self.ncuts, cutPoints))
        # chi square, chisq2 is used for writestats function
        self.chisq = num.zeros(shape=(self.nSensors, self.ncuts))
        # functions made from polyfit coefficients, accomodates drift
        self.polyfunc = num.zeros(shape=(self.nSensors, self.ncuts, cutPoints))
        self.polyfunc2 = num.zeros(shape=(self.nSensors, self.ncuts, cutPoints)) # not currently used for anything
        # cuts to keep after rejection
        self.goodCut = num.zeros(shape=(self.nSensors, 0)).tolist()
        
        
        
        for i in range(self.nSensors):
                # This sets up the function for only the first harmonic entered and will plot
                # using this function. The other harmonics and all stats for each are written
                # in the function writestats (found below)

            for j in range(self.ncuts):
                #This section divides the data into cuts
                beg = j * cutPoints #beg is for begin, as in, where in the data does this cut begin.
                end = (j + 1) * cutPoints

                x = self.timeArray[0:cutPoints]
                self.cutTime = x[:] #I don't recall why self.cutTime is created, it seems to be the same as x.

                self.data = self.newdat[i][beg:end]
                self.senssub[i][j] = self.data # self.senssub is data for sensor i, and cut j. (I assume I wrote this to divide the data into cuts)

            for j in range(self.ncuts):
                m = float(self.harmList[0])
                mHalf = int(m/2)
                mQuarter = int(m/4)

                # This section attempts to make each cut meet the next continuously. Varied amplitude makes it slightly discountinuous.
                if j == 0: # first cut
                    point0 = (self.senssub[i][j][0] + self.senssub[i][j][int(mHalf/self.sampleInterval)])/2
                    point1 = (self.senssub[i][j][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j+1][int(mQuarter/self.sampleInterval)])/2
                elif j == self.ncuts-1: # last cut
                    point0 = (self.senssub[i][j-1][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j][int(mQuarter/self.sampleInterval)])/2
                    point1 = (self.senssub[i][j][int(cutPoints-1 - (mHalf/self.sampleInterval))] + self.senssub[i][j][cutPoints-1])/2
                else: # intermediate cuts
                    point0 = (self.senssub[i][j-1][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j][int(mQuarter/self.sampleInterval)])/2
                    point1 = (self.senssub[i][j][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j+1][int(mQuarter/self.sampleInterval)])/2
                slope = (point1 - point0)/self.cutLength
            
                func=lambda t, a1, b1: slope * t + point0 + a1 * num.sin(2 * t * num.pi / float(m)) + b1 * num.cos(2 * t * num.pi / float(m))

                pmax = num.max(self.senssub[i][j])
                pmin = num.min(self.senssub[i][j])
                self.pinit = pmax - pmin
                p0 = num.array([self.pinit, self.pinit]) #why is this an array? what is this for?

                # Here the function is fit to the data using scipy curve_fit
                p, cov = curve_fit(func, x, self.senssub[i][j])
                
                    
                #p[0] and p[1] correspond to amplitudes a and b respectively
                #self.fitList creates a list of all a and b amplitudes for each cut and for each sensor
                #self.errList creates a list of a and b errors for each cut and each sensor
                self.polyfits[i][j] = [slope, point0]
                self.fitList[i][j] = [p[0], p[1]]
                try:
                    self.errList[i][j] = [num.sqrt(cov[0][0]), num.sqrt(cov[1][1])]
                except TypeError:
                    pass
                
                self.polyfunc[i][j] = self.polyfits[i][j][0] * x + self.polyfits[i][j][1]
                self.fittedFunc[i][j] = p[0] * num.sin(2 * x * num.pi / float(m)) + p[1] * num.cos(2 * x * num.pi / float(m)) + self.polyfunc[i][j]
                
                for k in range(cutPoints):
                    self.residuals[i][j][k] = self.fittedFunc[i][j][k] - self.senssub[i][j][k]

                for k in range(cutPoints):
                    self.chisq[i][j] += num.sqrt((self.fittedFunc[i][j][k] - self.senssub[i][j][k])**2 / (k+1))
                
                # Here the amp and error lists are being compiled into a form so that num.mean and num.std can be found

                self.aList[i][j] = self.fitList[i][j][0]
                self.aErrList[i][j] = self.errList[i][j][0]
                self.bList[i][j] = self.fitList[i][j][1]
                self.bErrList[i][j] = self.errList[i][j][1]

                # This is the old way to print out/calculate the amps and errors. Keeping just in case

##                '''print num.mean(self.aList[i]), num.std(self.aList[i]), num.mean(self.bList[i]), num.std(self.bList[i])
##                print self.chanNamesDict[i] + ' sin=', num.mean(self.aList[i]), 'sinErr=', num.std(self.aList[i]), 'cos=', num.mean(self.bList[i]), 'cosErr=', num.std(self.bList[i])'''

        # writestats function calculates for each frequency the mean a and b amplitudes, the standard
        # deviations for each, the overall average amplitude and the overall amplitude
        # errors and writes this to a file called statdata.
        # writestats also uses the lmfit package to fit the function instead of curve_fit
        
        

        for i in range(self.nSensors):
            a = self.nb.add(self.chanNamesDict[i])
            self.pageList.append(a)

            f = Figure(figsize=(4, 4), dpi=100)
            self.fitFigList.append(f)
            
            plt = f.add_subplot(311)
            x = self.timeArray[0:end]+self.initial*self.sampleInterval
            y = self.newdat[i][0:end]
            y2 = self.fittedFunc[i].flatten()
            y3 = self.polyfunc[i].flatten()

            plt.plot(x, y, 'b-')
            plt.plot(x, y2, 'y-')
            plt.plot(x, y3, 'g--')

            x2 = num.array(list(range(self.ncuts))) * self.cutLength + num.divide(self.cutLength, 2)+self.initial*self.sampleInterval
            plt.errorbar(x2, self.aList[i], yerr=self.aErrList[i], fmt='ro')

            x3 = num.array(list(range(self.ncuts))) * self.cutLength + num.divide(self.cutLength, 2)+self.initial*self.sampleInterval
            plt.errorbar(x2, self.bList[i], yerr=self.bErrList[i], fmt='go')

            if self.calOn.getBool():
                plt.set_title('Fitted Data')
                plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + ']')
            else:
                plt.set_title('Raw Fitted Data')
                plt.set_ylabel(self.chanNamesDict[i] + ' [Dimensionless]')

##            plt.set_xlabel('Time [s]')
            plt.set_xlim(min(x), max(x))
##            plt.set_ylim(min(y) - 0.05 * abs(min(y)), max(y) + 0.05 * abs(max(y)))
            
            plt2 = f.add_subplot(312)
                        
            y4 = self.residuals[i].flatten()
            plt2.plot(x, y4, 'r')
            plt2.set_xlabel('Time [s]')
            plt2.set_xlim(min(x), max(x))
            plt2.set_title('Residuals')

            plt3 = f.add_subplot(313)
            plt3.hist(self.chisq[i], color='r')
            plt3.set_xlabel('Chi Square')

            canvas = FigureCanvasTkAgg(f, a)
            canvas.show()
            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

            self.fitCanList.append(canvas)

            toolbar = NavigationToolbar2TkAgg(canvas, a)
            toolbar.update()
            canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

            cutEntryWdg = Wdg.FloatEntry(a,
                                            autoIsCurrent=False,
                                            isCurrent=True,
                                            defValue=None,
                                            width=35
            )
            cutEntryWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=cutEntryWdg))
            cutEntryWdg.bind('<Return>', self.FuncCall(self.setCutLevel, var=i))
            cutEntryWdg.pack()

            cutButton = Wdg.Button(
                master=a,
                text='Remove Cuts',
                width=20,
                command=self.FuncCall(self.cut, var=i),
            )
            cutButton.pack()

            self.fitWdgList.append([cutEntryWdg, cutButton])
            
        self.writestats() 

    def writestats(self):
        FitDir = "\\".join(self.dataDir.split('\\', -1)[:-2]+['\\fit\\']) #works on windows
  #      FitDir = os.path.join(*([os.sep]+self.dataDir.split(os.sep)[:-2]+['fit',""])) #platform-agnostic
        FileName = FitDir + 'Run' + str(self.runNum) + '.fit'
        stats = open(FileName, 'w')

        cutPoints = self.cutPoints
        
 #       pinit = self.pinit
        x = self.cutTime

        self.chisq2 = num.zeros(shape=(self.nSensors, self.ncuts))
        self.totals = [[[0 for m in range(len(self.harmList))] for n in range(2)] for i in range(self.nSensors)]
        self.amptotals = [[0 for m in range(len(self.harmList))] for i in range(self.nSensors)]
        self.totstd = num.zeros(shape=(self.nSensors, 2, len(self.harmList))).tolist()
        self.amptotstd = num.zeros(shape=(self.nSensors, len(self.harmList))).tolist()

        ## If I recall properly, the tcA, tcPhi, tcFactor are used to compute the timeconstant correction factor (hence tc).
        tc = float(self.LockinTimeConstant[0])
##        tcA = (1 + (self.attFreq**2)*(tc**2))
        tcA = 1#(1 + (((2*num.pi)/(100))**2)*(tc**2))
##        sA = (1+(((2*num.pi)/(100))**2)*(2**2))
##        tcPhi = 2/(num.tan(self.attFreq*tc))
##        tcFactor = tcA * (num.cos(tcPhi) + 1j * num.sin(tcPhi))

        for i in range(self.nSensors):
    ##                stats.write(str(self.chanNamesDict[i])+'\n\n')
            for j in range(self.ncuts):
                # Sets up the parameters to be used in calculating the amplitudes and errors.

                # Parameters is an object in the lmfit package. For this function, in some
                # instances parameters are used instead of variables. Further documentation can be found here:
                # cars9.uchicago.edu/software/python/lmfit/parameters.html
                
                pmax = num.max(self.senssub[i][j])
                pmin = num.min(self.senssub[i][j])
                self.pinit = pmax - pmin
                
                params=Parameters()
                for m in range(len(self.harmList)):
                    params.add('a'+str(m+1),value=self.pinit,vary=True)
                    params.add('b'+str(m+1),value=self.pinit,vary=True)
##                params.add('c',value=pinit,vary=True)
##                params.add('d',value=pinit,vary=True)

                # Sets up the function for lmfit

                mHarm = float(self.harmList[0])
                mHalf = mHarm/2
                mQuarter = mHarm/4

                if j == 0: # first cut
                    point0 = (self.senssub[i][j][0] + self.senssub[i][j][int(mHalf/self.sampleInterval)])/2
                    point1 = (self.senssub[i][j][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j+1][int(mQuarter/self.sampleInterval)])/2
                elif j == self.ncuts-1: # last cut
                    point0 = (self.senssub[i][j-1][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j][int(mQuarter/self.sampleInterval)])/2
                    point1 = (self.senssub[i][j][int(cutPoints-1 - (mHalf/self.sampleInterval))] + self.senssub[i][j][cutPoints-1])/2
                else: # intermediate cuts
                    point0 = (self.senssub[i][j-1][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j][int(mQuarter/self.sampleInterval)])/2
                    point1 = (self.senssub[i][j][int(cutPoints - mQuarter/self.sampleInterval)] + self.senssub[i][j+1][int(mQuarter/self.sampleInterval)])/2
                slope = (point1 - point0)/self.cutLength
##
##          I think these few lines have been replaced by the if/else above                
##                point0 = (self.senssub[i][j][0] + self.senssub[i][j][mHalf/self.sampleInterval])/2
##                point1 = (self.senssub[i][j][self.cutPoints-1] + self.senssub[i][j][self.cutPoints-1 - (mHalf/self.sampleInterval)])/2
##                slope = (point1 - point0)/self.cutLength

                def pfunc(params, t, ydata):
                    a=[]
                    b=[]
##                    c=[]
##                    d=[]
                    
                    for m in range(len(self.harmList)):
                        a.append(params['a'+str(m+1)].value)
                        b.append(params['b'+str(m+1)].value)
                        
##                    c.append(params['c'].value)
##                    d.append(params['d'].value)
                    
                    k=0
                            
                    model = a[0] * num.sin(2 * t * num.pi / float(mHarm)) + b[0] * num.cos(2 * t * num.pi / float(mHarm)) + slope * t + point0
                    if len(self.harmList)>1:
                        for m in self.harmList[1:]: ### needs to be fixed, harmonics are different ############
                            k+=1
                            model = model + a[k] * num.sin(2 * t * num.pi / float(mHarm)) + b[k] * num.cos(2 * t * num.pi / float(mHarm))
                    return model - self.senssub[i][j]

                # result, final and report_fit are all part of the lmfit package.
                # Further documentation can be found at:
                # cars9.uchicago.edu/software/python/lmfit/

                result = minimize(pfunc,params, args=(x,self.senssub[i][j]))
                final = self.senssub[i][j] + result.residual
##                report_fit(params)


                # am and bm are the a and b amplitudes respectively
                for m in range(len(self.harmList)):
                    am = num.asscalar(result.params['a'+str(m+1)].value)
                    bm = num.asscalar(result.params['b'+str(m+1)].value)

                    # Here each amplitude a (0) and b (1) is written to a list for each sensor(i), cut(j) and frequency (m).
                    self.fits[i][j][0][m] = am# * tcA # tcA is the timeConstant amplitude correction
                    self.fits[i][j][1][m] = bm# * tcA
                    a = num.array(self.fits)
##                self.polyfits2[i][j][0] = num.asscalar(params['c'].value)
##                self.polyfits2[i][j][1] = num.asscalar(params['d'].value)

                self.polyfunc2[i][j] = slope * x + point0
                
                fittedHarmonics = num.zeros(len(x))
                for m in range(len(self.harmList)):
                    harm = self.harmList[m]
                    fittedHarmonics += self.fits[i][j][0][m] * num.sin(2 * x * num.pi / float(harm)) + self.fits[i][j][1][m] * num.cos(2 * x * num.pi / float(harm))
                self.fittedFunc2[i][j] = slope * x + point0 + fittedHarmonics

                # Here we calculate the chi square for each cut
                for k in range(self.cutPoints):
                    self.chisq2[i][j] += num.sqrt((self.fittedFunc2[i][j][k] - self.senssub[i][j][k])**2 / (k+1))

                # avgcut[j][m] is the overall amplitude of each cut
                for m in range(len(self.harmList)):
                    self.avgcut[j][m] = num.sqrt(self.fits[i][j][0][m]**2 + self.fits[i][j][1][m]**2)

                # totals[i][n][m] is the totaled a(0) and b(1) amplitudes for each sensor(i) and frequency(m)
                # amptotals[i][m] is the sum of the total amplitudes for each sensor(i) and frequency(m)
                # if/else statement excludes rejected cuts from totals
                if j in self.goodCut[i] or self.goodCut[i] == []:
                    for m in range(len(self.harmList)):
                         self.amptotals[i][m] += self.avgcut[j][m]
                         for n in range(2):
                            self.totals[i][n][m] += self.fits[i][j][n][m]
                else:
                    continue


                
            # if/else statment is to account for rejected cuts
            if len(self.goodCut[i]) > 0:
                ncutsKept = len(self.goodCut[i])
            else:
                ncutsKept = self.ncuts
                
            # avgs[0][m] is the average a value for each harm m averaged over all cuts
            # avgs[1][m] is the average b value for each harm m averaged over all cuts
            # ampavg[0][m] is the overall amplitude for each harm m averaged over all cuts
                
            for m in range(len(self.harmList)):
                self.ampavg[0][m] = self.amptotals[i][m]/ncutsKept
                for n in range(2):
                    self.avgs[n][m] = self.totals[i][n][m]/ncutsKept

            # Here is where the standard deviation of the mean for a and b is calculated for each harm m
            ## old ## ampavg[1][m] is the error of the overall amp for each m averaged over all cuts
            # ampavg[1][m] is the standared deviation of the overall amp for each m averaged over all cuts
            
            for j in range(self.ncuts): #if/else accounts for rejected data, i.e., if the jth cut is 'good' then it is included, or if no data is rejected (self.goodCut[i] == []) then all are included.
                if j in self.goodCut[i] or self.goodCut[i] == []:
                    for m in range(len(self.harmList)):
                        self.amptotstd[i][m] += abs(self.avgcut[j][m] - self.ampavg[0][m])**2
                        self.ampavg[1][m] = num.sqrt(self.amptotstd[i][m]/((ncutsKept-1)))
                        for n in range(2):
                            self.totstd[i][n][m] += abs(self.fits[i][j][n][m] - self.avgs[n][m])**2
                            self.stds[n][m] = num.sqrt(self.totstd[i][n][m]/((ncutsKept-1)))
                else:
                    continue

            # Writes the file _[first harmonic] in the current directory

            #for m in range(len(self.harmList)):
                # Un-quote below to write the values from curve_fit package (as opposed to lmfit package) for comparison
    ##                    ''' stats.write('sin'+str(m+1)+'=' + str(num.mean(self.aList[i])) + ' sin'+str(m+1)+'Err=' + str(num.std(self.aList[i]))+'\n')
    ##                   stats.write('cos'+str(m+1)+'=' + str(num.mean(self.bList[i])) + ' cos'+str(m+1)+'Err=' + str(num.std(self.bList[i]))+'\n')'''
    ##                    '''stats.write('Period: '+str(self.harmList[m])+'\n')
    ##                    stats.write('sin' + str(m+1) + '=' + str(self.avgs[0][m]) + ' sin'+ str(m+1) + 'Err=' + str(self.stds[0][m]) + '\n')
    ##                    for j in range(self.ncuts):
    ##                        stats.write('sin' + str(m+1) + '=' + str(self.fits[i][j][0][m]) + '\n')
    ##                    stats.write('cos' + str(m+1) + '=' + str(self.avgs[1][m]) + ' cos'+ str(m+1) + 'Err=' + str(self.stds[1][m]) + '\n')
    ##                    for j in range(self.ncuts):
    ##                        stats.write('cos' + str(m+1) + '=' + str(self.fits[i][j][1][m]) + '\n')
    ##                    stats.write('amp' + str(m+1) + '=' + str(self.ampavg[0][m]) + ' amp' + str(m+1) + 'Err=' + str(self.ampavg[1][m]) + '\n\n')
    ##'''
            line = ''
            for m in range(len(self.harmList)):
                line += str(self.chanNamesDict[i]) + ': Period %s' % str(self.harmList[m]) + ', Torsion Filter: ' + str(self.TorsionOn.getBool()) + '\n\n'
    ##                    line += 'Torsion Filter: ' + str(self.TorsionOn.getBool()) + '\n\n'
                line += '[a amplitudes]' + '    \t' + '[b amplitudes]' + '  \t' + '[total amplitudes]' + '    \t' + '   [chi square]' + '\n'
                for j in range(self.ncuts):
                    if j in self.goodCut[i] or self.goodCut[i] == []: #this if/else statement gives us a way to distinguish rejected data in the .fit file. Rejected data does not contribute to stats.
                        line += str(self.fits[i][j][0][m]) + '    \t' + str(self.fits[i][j][1][m]) + '    \t' + str(self.avgcut[j][m]) + '    \t' \
                                + str(self.chisq2[i][j]) + '\n'
                    else:
                        line += '(' + str(self.fits[i][j][0][m]) + ')   \t(' + str(self.fits[i][j][1][m]) + ')   \t(' + str(self.avgcut[j][m]) + ')   \t(' \
                                + str(self.chisq2[i][j]) + ')' + '\n'
                line += '[Average amplitudes]' + '\n'
                line += str(self.avgs[0][m]) + '    \t' + str(self.avgs[1][m]) + '    \t' + str(self.ampavg[0][m]) + '\n'
                line += '[Standard deviations]' + '\n'
                line += str(self.stds[0][m]) + '    \t' + str(self.stds[1][m]) + '    \t' + str(self.ampavg[1][m]) + '\n\n'
                        
            stats.write(line)
            
        stats.close()
        return

    def setCutLevel(self, evt, var):
        i = var
        self.cutLevels[i] = self.fitWdgList[i][0].get()
        self.fitWdgList[i][0].setIsCurrent(True)

    def cut(self, var):
        i = var

        aList = list(self.aList[i])
        aErrList = list(self.aErrList[i])

        bList = list(self.bList[i])
        bErrList = list(self.bErrList[i])

        chisqList = list(self.chisq[i])

        if self.cutLevels[i] == 0:
            return

        for j in range(self.ncuts):
            if chisqList[j] > self.cutLevels[i]:
                self.removeList.append([i, j])


        a2List = [] # a amplitudes to keep
        a2ErrList = [] # a erros to keep
        b2List = [] # b amplitudes to keep
        b2ErrList = [] # b errors to keep
        chisq2List = []

        rList = [] # cuts to be removed

        for r in self.removeList:
            if r[0] == i:
                rList.append(r[1])

        for val in range(len(aList)):
            if val in rList:
                continue
            #else: This happens
            a2List.append(aList[val])
            a2ErrList.append(aErrList[val])
            b2List.append(bList[val])
            b2ErrList.append(bErrList[val])
            chisq2List.append(chisqList[val])
            
        page = self.nb.pagenames()[2 + i]
        # references page, rather than deleting it
        a = self.nb.page(2 + i)

        # references figure, rather than recreating it
        f = self.fitFigList[i]
        plt = f.add_subplot(311)
        plt.clear()

        y2 = self.fittedFunc[i].flatten()
        y3 = self.fittedFunc2[i].flatten()
        y4 = self.residuals[i].flatten()

        yfit = self.polyfunc[i].flatten()
        yfit2 = self.polyfunc2[i].flatten()
        y = self.newdat[i].flatten()[0:len(y2)]
        x = self.timeArray[0:len(y2)]
        plt.plot(x, y, 'b-')
        plt.plot(x, y2, 'y-')
        plt.plot(x, yfit, 'g--')
        
        self.goodCut[i] = []
        for val in range(self.ncuts):
            if val not in rList:
                self.goodCut[i].append(val)
        x2 = num.array(self.goodCut[i]) * self.cutLength + num.divide(self.cutLength, 2)
        plt.errorbar(x2, a2List, yerr=a2ErrList, fmt='ro')
        plt.errorbar(x2, b2List, yerr=b2ErrList, fmt='go')

        if self.calOn.getBool():
            plt.set_title('Fitted Data')
            plt.set_ylabel(self.chanNamesDict[i] + ' ' + '[' + self.cal1Dict[i][1] + ']')
        else:
            plt.set_title('Raw Fitted Data')
            plt.set_ylabel(self.chanNamesDict[i] + ' [V]')

        plt.set_xlabel('Time [s]')
        plt.set_xlim(min(x), max(x))
##        plt.set_ylim(min(y) - 0.05 * abs(min(y)), max(y) + 0.05 * abs(max(y)))
        
        plt2 = f.add_subplot(312)
        plt2.plot(x, y4, 'r')
        plt2.set_xlabel('Time [s]')
        plt2.set_xlim(min(x), max(x))
        plt2.set_title('Residuals')
    
        plt3 = f.add_subplot(313)
        plt3.clear()
        plt3.hist(b2ErrList, color='g')
        plt3.hist(a2ErrList, color='r')
        plt3.hist(chisq2List, color='r')
        
        plt3.set_xlabel('Chi Square')
        plt3.set_ylabel('# of Cuts')

        # references canvas, rather than recreating it
        canvas = self.fitCanList[i]
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

##        toolbar = NavigationToolbar2TkAgg(canvas, a)
##        toolbar.update()
##        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
##
##        cutEntryWdg = Wdg.FloatEntry(a,
##                                        autoIsCurrent=False,
##                                        isCurrent=True,
##                                        defValue=self.cutLevels[i],
##                                        width=35
##        )
##        cutEntryWdg.bind('<KeyPress>', self.FuncCall(self.setBackground, var=cutEntryWdg))
##        cutEntryWdg.bind('<Return>', self.FuncCall(self.setCutLevel, var=i))
##        cutEntryWdg.pack()
##
##        cutButton = Wdg.Button(
##            master=a,
##            text='Remove Cuts',
##            width=20,
##            command=self.FuncCall(self.cut, var=i),
##        )
##        cutButton.pack()
##
##        self.fitWdgList[i] = [cutEntryWdg, cutButton]
##
        print((page, ' sin=', num.mean(a2List), 'sinErr=', num.std(a2List), 'cos=', num.mean(
            b2List), 'cosErr=', num.std(b2List)))

        self.cutAttempt = 1

        self.writestats()
        
##        FitDir = "\\".join(self.dataDir.split('\\', -1)[:-2]+['\\fit\\'])
##        FileName = FitDir + 'Run' + str(self.runNum) + '.fit'
##        stats = open(FileName, 'w')
        

    def setHarm(self, evt):

        s = self.harmWdg.get()

        if len(s) > 0:
            self.harmList = list(map(float, s.split()))
            self.nHarms = len(self.harmList)

        else:
            msg.showerror("Error", "Invalid Harmonics Entered")
            return

        self.harmWdg.setIsCurrent(True)

    def setBackground(self, evt, var):
        var.setIsCurrent(False)

    def _quit(self):
        root.quit()
        root.destroy()


def vp_start_gui():
    global root
    root = Tk.Tk()
    root.geometry("1200x900+200+50")
    root.title("Run Analysis")
    w = NoteBook(root)
    root.mainloop()


if __name__ == '__main__':
    def vp_start_gui():
        global root
        root = Tk.Tk()
        root.geometry("1200x900+200+50")
        root.title("Run Analysis")
        w = NoteBook(root)
        root.mainloop()

    vp_start_gui()

################ inertia correction, insert in code later #########################
##        eTau = 216792
##        self.attFreq = (2*num.pi)/(201.8)
##        self.torFreq = (2*num.pi)/self.ttor
##        qualFac = (eTau * num.pi)/(self.ttor * self.sampleInterval)
##        iA = num.sqrt((1-(self.attFreq/self.torFreq)**2)**2 + (1/qualFac)**2) # check syntax
##        iPhi = 1/(num.tan((self.torFreq**2)/(qualFac*(self.torFreq**2 - self.attFreq**2)))) # check syntax
##        iFactor = iA * (num.cos(iPhi) + 1j * num.sin(iPhi))
################ time constant correction #########################################
##    tc = self.LockinTimeConstant
##    tcA = (1 + (self.attFreq**2)*(tc**2)
##    tcPhi = 2/(num.tan(self.attFreq*tc))
##    tcFactor = tcA * (num.cos(tcPhi) + 1j * num.sin(tcPhi)
