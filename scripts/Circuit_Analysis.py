import numpy as np
import matplotlib.pyplot as plt
import csv
from lmfit import Parameters, minimize
import pandas as pd
from tkinter import filedialog, Tk

# dataframe print formatting
pd.options.display.float_format = '{:.2f}'.format

# can uncomment this to force inline plotting if desired. ignore error message
# %matplotlib inline

parameters = Parameters()

# this script assumes that you are going to supply it with amplitude vs frequency (Hz) data collected using the same circuit setup as us (https://doi.org/10.1002/pol.20210763)

# make sure that all imported csv files are the basic comma delimited file format, and not another csv format


# TODO
###############################
# True - get values of Δf and ΔΓ at each harmonic given bare QCM and with-film data; False - get values of fn and Γn at each harmonic given either bare QCM or with-film data
getShifts = False
# if you're not collecting shifts, specify here whether the system is bare or has a film on it
isBare = True
# the highest harmonic included in the data set; odd numbers only; current plotting setup restricts highestHarmonic to a max of 15 - feel free to go into the plotting setup and change this if needed
highestHarmonic = 11
# tune this to optimize fitting range for the data set. it is the percentage of the trace from left to right to include for fitting. you want to avoid influence of side peaks in your data to fit
fitInclusion = 100
# use this when using a universal fitInclusion is not working well. each value will default to the above but you can change them each individually
cutIndividually = False

# parameter guesses and bounds are specified below; they have been fairly well optimized for most data, but you may need to adjust them occasionally to get good fits

# the resonance frequency at harmonic n; input guess does not matter as the guess is handled automatically later
parameters.add('fn', value = 0, min = 0, max = 1e8, vary = True)
# dissipation at the above resonance; lower bound must be above zero to prevent division by zero
parameters.add('dissip', value = 20, min = 0.0001, max = 1e4, vary = True)
# background signal, which creates a varied small downward shift in all of the data
parameters.add('background', value = -0.3, min = -1, max = 0, vary = True)
# R_2/R_m in the circuit; R_2 = 17.61 ohms; R_m seems to increase with higher harmonics
parameters.add('An', value = 17.61/40, min = 0, max = 17.61/5, vary = True)
# 2*pi*R_2*Z_0 in the circuit; should be a very small number, typically negative
parameters.add('Bn', value = -1e-9, min = -1e-8, max = 1e-9, vary = True)

# these variables tune the way that parameter guesses and bounds change with increasing harmonics, as their behavior is relatively predictable and it helps with fitting; as with above, these should already be well optimized
dissipGuessCoef = 1.5 # a given harmonic gets a dissipation initial guess of dissipGuessCoef times the previous harmonic's fit dissipation
dissipMinCoef = 0.05 # a given harmonic is assumed to have dissipation at least dissipMinCoef times the previous harmonic's fit dissipation
AnGuessCoef = 0.7 # a given An value gets an An initial guess of AnGuessCoef times the previous harmonic's fit An
fnBoundingScale = 100 # (Hz) the resonance frequency is assumed to be within fnBoundingScale Hz of the frequency at which the highest amplitude data point occurs
###############################


# array of harmonics
harmonics = []
for i in range(1, highestHarmonic + 1, 2):
    harmonics.append(i)
numHarmonics = len(harmonics)
    
# setting the chosen default fit inclusion percentage for each data set
fitInclusionIndivBare, fitInclusionIndivFilm = [], []
for i in range(0, numHarmonics):
    fitInclusionIndivBare.append(fitInclusion)
    fitInclusionIndivFilm.append(fitInclusion)
    
# TODO
# set any desired individual fitInclusion values here if cutIndividually is chosen
###############################
if cutIndividually:
    # for example: 
    fitInclusionIndivFilm[harmonics.index(11)] = 55
    # specify more values as needed here
###############################

harmonicNames = []
for i in range(numHarmonics):
    harmonicNames.append('n = ' + str(harmonics[i]))

# constant in circuit equation
C = 1.41

# circuit equation; input: frequencies, fitting parameters; output: amplitude
# n subscripts are used here to distinguish values which have different values at each resonance; f is the independent frequency variable while fn is a resonance frequency at harmonic n
# note that these are not the same A and B constants seen in associated modulus modeling
def circuitEqtn(f, fn, dissip, background, An, Bn):   
    return background + 10*np.log10(C**2 + (Bn**2 * f**2) + ((An**2 + (2*C*An) - ((Bn*An*(f**2 - fn**2))/dissip))/
                                              (1 + ((f**2 - fn**2)**2)/(4*dissip**2 * f**2))))

# a function that takes a set of parameters and returns a list of all of the residuals associated with each guess; this is called repeatedly by the fitting algorithm as it varies parameter values
def residuals(paramGuesses):        
    paramGuessVals = list(paramGuesses.valuesdict().values())
    
    # guess values are sent through the circuit equation; fList and ampList will always already be filled with the appropriate data whenever this function is called
    residualList = []
    for i in range(len(fList)):
        residualList.append(float(circuitEqtn(fList[i], paramGuessVals[0], paramGuessVals[1], paramGuessVals[2], paramGuessVals[3], paramGuessVals[4]) - ampList[i]))

    return residualList

# bare QCM
if getShifts or isBare:
    
    # reading data files
    
    # select the file for the first harmonic; the files should only be named '1.csv', '3.csv', and so on; it is recommended to have a folder for each set of resonance traces
    print('\nSelect the file for the first harmonic of your bare QCM data. The rest will be read automatically.\n(files should be named \'1.csv\', \'3.csv\', ...)\n')
    
    # open tkinter window
    root = Tk()
    # setting tkinter windows to appear at the front
    root.attributes("-topmost", True)
    # open file dialog
    firstFilePath = filedialog.askopenfilename(parent = root, title = 'Open Bare QCM Resonance Trace Data')
    # remove tkinter window
    root.destroy()
    
    if firstFilePath == '':
        print('\nFile selection canceled.')
        raise SystemExit()
    # full resonance traces for all harmonics, to compare to fit
    freqFullBareQCM, ampFullBareQCM = [], []
    # resonance trace portions used for fitting for all harmonics
    freqToFitBareQCM, ampToFitBareQCM = [], []
     
    for i in range(0, numHarmonics):
        
        # temporary arrays holding individual harmonic's data each loop
        freqFull, ampFull, freqToFit, ampToFit = [], [], [], []
        
        # reads the file corresponding to the current harmonic
        filePath = firstFilePath.replace('1.csv', str(harmonics[i]) + '.csv')
        file = csv.reader(open(filePath), delimiter = ",")
        
        # extracting current harmonic's trace from file
        for row in file:
            if row[0] != '':
                freqFull.append(float(row[0]))
                ampFull.append(float(row[1]))
                
        # a simple means of not including side peaks in fit
        cutoff = int(len(freqFull)*(fitInclusionIndivBare[i]/100))
        for j in range(0, cutoff):
            freqToFit.append(freqFull[j])
            ampToFit.append(ampFull[j])
        
        freqFullBareQCM.append(freqFull)
        ampFullBareQCM.append(ampFull)
        freqToFitBareQCM.append(freqToFit)
        ampToFitBareQCM.append(ampToFit)
        
    # fitting
    
    # holds full set of fit parameter values
    fitTableBareQCM = []

    # holds fit fn and Γn values
    fitFreqBareQCM, fitDissipBareQCM = [], []
    
    # the figure that will hold plots of each harmonic's fit
    figBare = plt.figure(figsize = (16, 10))
    figBare.suptitle('Bare QCM')
    
    print('Running fitting algorithm...\n')
    for i in range(0, numHarmonics):  
        
        fList = []
        ampList = []
        fnGuess = 0
        maxAmp = max(ampToFitBareQCM[i])
               
        for j in range(0, len(freqToFitBareQCM[i])):
            
            # fList and ampList are filled with the supplied data
            fList.append(freqToFitBareQCM[i][j])
            ampList.append(ampToFitBareQCM[i][j])
            
            # the guess for fn is set to the frequency of the data point at which amplitude is highest; this will be quite close to the fit value assuming you collect a lot of data points
            if ampToFitBareQCM[i][j] == maxAmp:
                fnGuess = freqToFitBareQCM[i][j]
        
        # parameter values and bounds are adjusted according to the input adjustment variables here
        parameters['fn'].min = fnGuess - fnBoundingScale
        parameters['fn'].max = fnGuess + fnBoundingScale
        parameters['fn'].value = fnGuess
        if not i == 0:
            parameters['dissip'].value = fitDissipBareQCM[i-1] * dissipGuessCoef
            parameters['dissip'].min = fitDissipBareQCM[i-1] * dissipMinCoef
            parameters['An'].value = fitTableBareQCM[i-1][3] * AnGuessCoef
        
        # this is where the fitting algorithm is run; it minimizes the sum of the residuals squared
        resultBareQCM = minimize(residuals, parameters)
        
        # fit parameter values
        fitTableBareQCM.append(list(resultBareQCM.params.valuesdict().values()))
        
        fitFreqBareQCM.append(fitTableBareQCM[i][0])
        fitDissipBareQCM.append(fitTableBareQCM[i][1])
         
        # plotting resonance trace fits
    
        axBare = plt.subplot(2, 4, i + 1)
        
        # f(x)_theoretical using best fit parameter values in circuit equation
        solvedAmp = []
        for j in range(0, len(freqFullBareQCM[i])):
            solvedAmp.append(circuitEqtn(freqFullBareQCM[i][j], fitTableBareQCM[i][0], fitTableBareQCM[i][1],
                                                 fitTableBareQCM[i][2], fitTableBareQCM[i][3], fitTableBareQCM[i][4]))
        
        axBare.scatter(freqFullBareQCM[i], ampFullBareQCM[i], color = 'skyblue', label = 'full trace')
        axBare.scatter(freqToFitBareQCM[i], ampToFitBareQCM[i], color = 'black', label = 'points used to fit')
        axBare.plot(freqFullBareQCM[i], solvedAmp, label='best fit', color = 'red')
        axBare.set(xlabel = 'Frequency', ylabel = 'Amplitude', title = 'n = '+str(harmonics[i]))
    
    # make legend, show figure
    handles, labels = axBare.get_legend_handles_labels()
    plt.figlegend(handles, labels, loc = 'lower right', fontsize = 15)
    plt.tight_layout()
    plt.show(); 

    #print fundamental frequency of your QCM
    print('\nFundamental frequency f_f = %.6f MHz' %(fitFreqBareQCM[0]*1e-6)) 

    # proceeding with final plotting and option to save fn and Γn, given choice to only analyze bare QCM
    if not getShifts:
        
        # plotting fn and Γn
        plt.scatter(harmonics, fitFreqBareQCM, color = 'black')
        plt.xlabel('Harmonic n')
        plt.ylabel('fn')
        plt.show();
        plt.scatter(harmonics, fitDissipBareQCM, color = 'black')
        plt.xlabel('Harmonic n')
        plt.ylabel('Γn')
        plt.show();
            
        # initialize data matrix with zeros; you will have the option to save this as a csv file after it is filled in
        dataToSaveBareQCM = [[0 for i in range(2)] for j in range(numHarmonics)]
        
        # fill data matrix with frequency and dissipation values, print them
        for i in range(0, numHarmonics):
            dataToSaveBareQCM[i][0] = fitFreqBareQCM[i]
            dataToSaveBareQCM[i][1] = fitDissipBareQCM[i]
        print('\n')
        print(pd.DataFrame(dataToSaveBareQCM, columns=['fn (Hz)', 'Γn (Hz)'], index = harmonicNames))
        
        # save data to csv
        print('\nYou may save this data table as a csv file, or just hit \'Cancel\'.')
         
        # open tkinter window
        root = Tk()
        # setting tkinter windows to appear at the front
        root.attributes("-topmost", True)
        # open file dialog
        fileName = filedialog.asksaveasfilename(parent = root, defaultextension = ".csv", title = 'Save f_n, Gamma_n if Needed')
        # remove tkinter window
        root.destroy()
        
        if not fileName == '':
            pd.DataFrame(dataToSaveBareQCM).to_csv(fileName, header=['Frequency (Hz)', 'Dissipation (Hz)'], index = None)

# for measurements with a film atop the QCM
if getShifts or not isBare:
    
    # reading data files
    
    # select the file for the first harmonic; the files should only be named '1.csv', '3.csv', and so on; it is recommended to have a folder for each set of resonance traces
    print('\nSelect the file for the first harmonic of your film-included data. The rest will be read automatically.\n(files should be named \'1.csv\', \'3.csv\', ...)\n')
    
    # open tkinter window
    root = Tk()
    # setting tkinter windows to appear at the front
    root.attributes("-topmost", True)
    # open file dialog
    firstFilePath = filedialog.askopenfilename(parent = root, title = 'Open QCM+Film Resonance Trace Data')
    # remove tkinter window
    root.destroy()
    
    if firstFilePath == '':
        print('\nFile selection canceled.')
        raise SystemExit()
    # full resonance traces for all harmonics, to compare to fit
    freqFullWithFilm, ampFullWithFilm = [], []
    # resonance trace portions used for fitting for all harmonics
    freqToFitWithFilm, ampToFitWithFilm = [], []
     
    for i in range(0, numHarmonics):
        # temporary arrays holding individual harmonic's data each loop
        freqFull, ampFull, freqToFit, ampToFit = [], [], [], []
        
        # reads the file corresponding to the current harmonic
        filePath = firstFilePath.replace('1.csv', str(harmonics[i]) + '.csv')
        file = csv.reader(open(filePath), delimiter = ",")
        
        # extracting current harmonic's trace from file
        for row in file:
            if row[0] != '':
                freqFull.append(float(row[0]))
                ampFull.append(float(row[1]))
                
        # a simple means of not including side peaks in fit
        cutoff = int(len(freqFull)*(fitInclusionIndivFilm[i]/100))
        for j in range(0, cutoff):
            freqToFit.append(freqFull[j])
            ampToFit.append(ampFull[j])
        
        freqFullWithFilm.append(freqFull)
        ampFullWithFilm.append(ampFull)
        freqToFitWithFilm.append(freqToFit)
        ampToFitWithFilm.append(ampToFit)
        
    # fitting
    
    # holds full set of fit parameter values
    fitTableWithFilm = []

    # holds fit fn and Γn values
    fitFreqWithFilm, fitDissipWithFilm = [], []
    
    # the figure that will hold plots of each harmonic's fit
    figFilm = plt.figure(figsize = (16, 10))
    figFilm.suptitle('With Film')
    
    print('Running fitting algorithm...\n')
    for i in range(0, numHarmonics):  
        
        fList = []
        ampList = []
        fnGuess = 0
        maxAmp = max(ampToFitWithFilm[i])
        
        for j in range(0, len(freqToFitWithFilm[i])):
            
            # fList and ampList are filled with the supplied data
            fList.append(freqToFitWithFilm[i][j])
            ampList.append(ampToFitWithFilm[i][j])
            
            # the guess for fn is set to the frequency of the data point at which amplitude is highest; this will be quite close to the fit value assuming you collect a lot of data points
            if ampToFitWithFilm[i][j] == maxAmp:
                fnGuess = freqToFitWithFilm[i][j]
                
        # parameter values and bounds are adjusted according to the input adjustment variables here
        parameters['fn'].min = fnGuess - fnBoundingScale
        parameters['fn'].max = fnGuess + fnBoundingScale
        parameters['fn'].value = fnGuess
        if not i == 0:
            parameters['dissip'].value = fitDissipWithFilm[i-1] * dissipGuessCoef
            parameters['dissip'].min = fitDissipWithFilm[i-1] * dissipMinCoef
            parameters['An'].value = fitTableWithFilm[i-1][3] * AnGuessCoef
        # resetting these attributes for the first harmonic in case shifts are being obtained and thus the code for the bare QCM was run first
        else:
            parameters['dissip'].value = 20
            parameters['dissip'].min = 0.0001
            parameters['An'].value = 17.61/40

        # this is where the fitting algorithm is run; it minimizes the sum of the residuals squared
        resultWithFilm = minimize(residuals, parameters)
        
        # fit parameter values
        fitTableWithFilm.append(list(resultWithFilm.params.valuesdict().values()))
        
        fitFreqWithFilm.append(fitTableWithFilm[i][0])
        fitDissipWithFilm.append(fitTableWithFilm[i][1])
         
        # plotting resonance trace fits
    
        axFilm = plt.subplot(2, 4, i + 1)
        
        # f(x)_theoretical using best fit parameter values in circuit equation
        solvedAmp = []
        for j in range(0, len(freqFullWithFilm[i])):
            solvedAmp.append(circuitEqtn(freqFullWithFilm[i][j], fitTableWithFilm[i][0], fitTableWithFilm[i][1],
                                                 fitTableWithFilm[i][2], fitTableWithFilm[i][3], fitTableWithFilm[i][4]))
        
        axFilm.scatter(freqFullWithFilm[i], ampFullWithFilm[i], color = 'skyblue', label = 'full trace')
        axFilm.scatter(freqToFitWithFilm[i], ampToFitWithFilm[i], color = 'black', label = 'points used to fit')
        axFilm.plot(freqFullWithFilm[i], solvedAmp, label='best fit', color = 'red')
        axFilm.set(xlabel = 'Frequency', ylabel = 'Amplitude', title = 'n = '+str(harmonics[i]))
        
    # make legend, show figure
    handles, labels = axFilm.get_legend_handles_labels()
    plt.figlegend(handles, labels, loc = 'lower right', fontsize = 15)
    plt.tight_layout()
    plt.show(); 

    # proceeding with final plotting and option to save fn and Γn, given choice to only analyze QCM with film
    if not getShifts:
        
        # plotting fn and Γn
        plt.scatter(harmonics, fitFreqWithFilm, color = 'black')
        plt.xlabel('Harmonic n')
        plt.ylabel('fn')
        plt.show();
        plt.scatter(harmonics, fitDissipWithFilm, color = 'black')
        plt.xlabel('Harmonic n')
        plt.ylabel('Γn')
        plt.show();
            
        # initialize data matrix with zeros; you will have the option to save this as a csv file after it is filled in
        dataToSaveWithFilm = [[0 for i in range(2)] for j in range(numHarmonics)]
        
        # fill data matrix with frequency and dissipation values, print them
        for i in range(0, numHarmonics):
            dataToSaveWithFilm[i][0] = fitFreqWithFilm[i]
            dataToSaveWithFilm[i][1] = fitDissipWithFilm[i]
        print('\n')
        print(pd.DataFrame(dataToSaveWithFilm, columns=['fn (Hz)', 'Γn (Hz)'], index = harmonicNames))
        
        # save data to csv
        print('\nYou may save this data table as a csv file, or just hit \'Cancel\'.')
        fileName = ''
        
        # open tkinter window
        root = Tk()
        # setting tkinter windows to appear at the front
        root.attributes("-topmost", True)
        # open file dialog
        fileName = filedialog.asksaveasfilename(parent = root, defaultextension = ".csv", title = 'Save f_n, Gamma_n if Needed')
        # remove tkinter window
        root.destroy()
        
        if not fileName == '':
            pd.DataFrame(dataToSaveWithFilm).to_csv(fileName, header=['Frequency (Hz)', 'Dissipation (Hz)'], index = None)

# for analysis of the shifts in data from a bare QCM to a QCM with a film on it
if getShifts:
    
    # Δf and ΔΓ at each harmonic; we drop the n subscripts at this point as it is implied that these variables describe harmonic-dependent resonances
    freqShifts, dissipShifts = [], []
    
    # initialize data matrix with zeros; you will have the option to save this as a csv file after it is filled in
    shiftsToSave = [[0 for i in range(2)] for j in range(numHarmonics)]
    
    # subtracting bare QCM fn and Γn from the values with a film, filling these shifts into the above data structures
    for i in range(0, numHarmonics):
        freqShifts.append(fitFreqWithFilm[i] - fitFreqBareQCM[i])
        dissipShifts.append(fitDissipWithFilm[i] - fitDissipBareQCM[i])
        shiftsToSave[i][0] = freqShifts[i]
        shiftsToSave[i][1] = dissipShifts[i]
        
    # plotting and printing Δf and ΔΓ
    plt.scatter(harmonics, freqShifts, color = 'black')
    plt.xlabel('Harmonic n')
    plt.ylabel('Δf')
    plt.show();
    plt.scatter(harmonics, dissipShifts, color = 'black')
    plt.xlabel('Harmonic n')
    plt.ylabel('ΔΓ')
    plt.show();
    print('\n')
    print(pd.DataFrame(shiftsToSave, columns=['Δf (Hz)', 'ΔΓ (Hz)'], index = harmonicNames))
    
    # save data to csv
    print('\nYou may save this data table as a csv file, or just hit \'Cancel\'.')
    fileName = ''
    
    # open tkinter window
    root = Tk()
    # setting tkinter windows to appear at the front
    root.attributes("-topmost", True)
    # open file dialog
    fileName = filedialog.asksaveasfilename(parent = root, defaultextension = ".csv", title = 'Save Delta f, Delta Gamma if Needed')
    # remove tkinter window
    root.destroy()
        
    if not fileName == '':
        pd.DataFrame(shiftsToSave).to_csv(fileName, header=['Frequency Shift (Hz)', 'Dissipation Shift (Hz)'], index = None)