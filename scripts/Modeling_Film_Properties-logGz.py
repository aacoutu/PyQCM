from sympy.interactive import printing
printing.init_printing(use_latex=True)
import sympy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import filedialog, Tk
import csv
import mpmath as mp
from lmfit import Parameters, minimize
from IPython.display import display # can use to display the symbolic math

# dataframe print formatting
pd.options.display.float_format = '{:.2f}'.format

# can uncomment this to force inline plotting if desired. ignore error message
# %matplotlib inline

# make sure that all imported csv files are the basic comma delimited file format, and not another csv format

parameters = Parameters()


# TODO
###############################
# number of layers to model the film as; if using m > 2, it is best to stick with odd numbers if you want a layer cutoff to exist at the compositional interface
m = 2
# the highest harmonic included in the data set; odd numbers only
highestHarmonic = 7
# fundamental frequency of the bare QCM used
fundFreq = 4.999*1e6 # Hz
# you can use a Levenberg-Marquardt (LM) fitting algorithm to find best fit parameter values; if set to False, you will just get the result of the input parameter value guesses
performFitting = False

# here you can input guesses for parameter values; if you are performing fitting, then here is also where you will set which parameters should vary and any bounds you'd like to impose
# if you are doing a single-layer analysis (m = 1) then use the parameters associated with the top of the film
# allowing several varied parameters while modeling with several layers is not advised
# do not alter the order of these parameter initializations, as the order they are passed to later functions will matter and be affected by initialization order

parameters.add('hTop', value = 0.380, min = 0.01, max = None, vary = False) # (m^-6) thickness of the top polymer layer
parameters.add('rhoTop', value = 895, min = 0.01, max = None, vary = False) # (kg m^-3) density of the top polymer
parameters.add('GfundStorTop', value = 2, min = 0.01, max = None, vary = False) # (MPa) storage modulus of the top polymer at the first harmonic (near the fundamental frequency of the QCM) 
parameters.add('GfundLossTop', value = 3.3, min = 0, max = None, vary = False) # (MPa) loss modulus of the top polymer at the first harmonic (near the fundamental frequency of the QCM) 
parameters.add('betaStorTop', value = 0.5, min = -2, max = 2, vary = False) # frequency dependence of the storage modulus of the top polymer
parameters.add('betaLossTop', value = 0.74, min = -2, max = 2, vary = False) # frequency dependence of the loss modulus of the top polymer

parameters.add('hBot', value = 1.360, min = 0.01, max = None, vary = False) # (m^-6) thickness of the bottom polymer layer
parameters.add('rhoBot', value = 1040, min = 0.01, max = None, vary = False) # (kg m^-3) density of the bottom polymer
parameters.add('GfundStorBot', value = 1900, min = 0.01, max = None, vary = False) # (MPa) storage modulus of the bottom polymer at the first harmonic (near the fundamental frequency of the QCM) 
parameters.add('GfundLossBot', value = 70, min = 0, max = None, vary = False) # (MPa) loss modulus of the bottom polymer at the first harmonic (near the fundamental frequency of the QCM) 
parameters.add('betaStorBot', value = 0, min = -2, max = 2, vary = False) # frequency dependence of the storage modulus of the bottom polymer
parameters.add('betaLossBot', value = 0, min = -2, max = 2, vary = False) # frequency dependence of the loss modulus of the bottom polymer

# our focus is on the complex shear modulus profile, however one could theoretically use these methods to study local density or independently study different local compenents of the shear modulus in a range of systems
parameters.add('wG', value = 0.1, min = 0.1, max = 500, vary = True) # (nm) shear modulus gradient width; must be greater than zero
parameters.add('muG', value = 0, min = -500, max = 500, vary = False) # (nm) asymmetry; i.e. offset the gradient's center has from the polymer-polymer interface; positive is farther from the top of the sample

parameters.add('wRho', value = 0.1, min = 0.1, max = 20, vary = False) # (nm) density gradient width
parameters.add('muRho', value = 0, min = -5, max = 5, vary = False) # (nm) density gradient asymmetry

# how much of the bottom side polymer to allow sliced into the gradient model for m > 2; this optimizes resolution where it matters and makes it easier to ensure a layer cutoff exists at the compositional interface
# if ensuring such a layer cutoff is important, be sure to set m to an odd number; from there, it's typically easiest to set botSideAmount equal to parameters['hTop'].value, but alternatives can achieve the same
# example alternative: m = 5, botSideAmount = 1/3 or 3 times parameters['hTop'].value; 1 layer is used for the bottom, 4 remain, 1 is used for one side of the interface, 3 for the other side
# if you allow hTop to vary, botSideAmount will be imposed later on to always equal hTop in order to preserve the aforementioned expected behavior
botSideAmount = parameters['hTop'].value # (m^-6)

# if your input experimental data is of the form Δf - Δf_t=0 instead of just Δf, set analyzeChange to True and enter the modeled Δf_t=0 and ΔΓ_t=0 values corresponding to your model of the initial state of the system
analyzeChange = False
freqShiftInitList = [-9252.06, -27922.43, -46945.33, -66447.06]
dissipShiftInitList = [18.54, 262.07, 912.92, 2112.36]
###############################


# array of harmonics
harmonics = []
for i in range(1, highestHarmonic + 1, 2):
    harmonics.append(i)
numHarmonics = len(harmonics)

# select the csv file for a single sample output by the Circuit Analysis code, or a similarly formatted csv file using averages from a set of samples and uncertainties included in the 3rd and 4th columns
print('\nSelect the csv file.\n\nFormat: rows - optional single header row, then harmonic order; columns - Δf (Hz), ΔΓ (Hz), uncertainty in Δf, uncertainty in ΔΓ\n\nIf analyzing data from a single sample, just leave the last two columns blank.\n')

# open tkinter window
root = Tk()
# setting tkinter windows to appear at the front
root.attributes("-topmost", True)
# open file dialog
filePath = filedialog.askopenfilename(parent = root, title = 'Open Delta f, Delta Gamma data')
# remove tkinter window
root.destroy()

if filePath == '':
    print('\nFile selection canceled.')
    raise SystemExit()
file = csv.reader(open(filePath), delimiter = ",")
# arrays to hold the data to be extracted from the csv file
freqShiftList, dissipShiftList, freqShiftUncList, dissipShiftUncList = [], [], [], []
 
# function to identify if a string is representable as a float
def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

# reading the selected file
rowIndex = 0
hasHeader = False
for row in file:
    # skip header row if there is one
    if rowIndex == 0 and not isFloat(row[0]):
        hasHeader = True
        rowIndex = rowIndex + 1
        continue
    
    if (hasHeader and rowIndex > numHarmonics) or (not hasHeader and rowIndex >= numHarmonics):
        # extra data exists but is being skipped due to highestHarmonic setting
        if row[0] != '':
            break
        # highestHarmonic setting is demanding data not included
        else:
            print('Error: insufficient data in csv file for highestHarmonic value entered.')
            raise SystemExit()
    
    freqShiftList.append(float(row[0]))
    dissipShiftList.append(float(row[1]))
    
    # if only two columns are present, a value of 1 will be filled in for every standard deviation; this means that every Δf and ΔΓ value will be weighted equally
    if len(row) == 2:
        freqShiftUncList.append(1)
        dissipShiftUncList.append(1)
    else:
        freqShiftUncList.append(float(row[2]))
        dissipShiftUncList.append(float(row[3]))
        
    rowIndex = rowIndex + 1

# many variables will have a symbol and a value; the value is a number, while the symbol is used for working with mathematical relations; here a few symbols are created
omega, n, resShift = sp.symbols('ω,n,\widetilde{Δf}')

# setting the mathematical relation of symbols that produces the value of ω; if your input experimental data is of the form (Δf - Δf_t=0) instead of just Δf, it is here that the code takes that into account
if analyzeChange:
    freqShiftInit, dissipShiftInit = sp.symbols('Δf_{t=0}, ΔΓ_{t=0}')
    # under these circumstances the variable, resShift, simply takes on the meaning of (Δf - Δf_t=0) + i(ΔΓ - ΔΓ_t=0) here instead of Δf + iΔΓ
    omegaVal = 2*np.pi*(n*fundFreq + freqShiftInit + dissipShiftInit*1j + resShift) 
else:
    omegaVal = 2*np.pi*(n*fundFreq + resShift)

# symbols which take on different values at different layers are initialized as arrays containing zeros here; polymer film: m layers (0 through m-1), QCM: layer m, total layer number: m+1
rho = [0 for i in range(m+1)] # density
G = [0 for i in range(m+1)] # complex shear modulus
GfundStor = [0 for i in range(m+1)] # storage modulus (real part) at first harmonic; only exists for convenient extraction of local storage modulus later
GfundLoss = [0 for i in range(m+1)] # loss modulus (imaginary part) at first harmonic; only exists for convenient extraction of local loss modulus later

z = [0 for i in range(m+1)] # position at top of layer (z-axis starts at top of sample, positive direction going down into the sample)
k = [0 for i in range(m+1)] # wavenumber
imped = [0 for i in range(m+1)] # acoustic impedance

# initializing arrays to hold the values associated with each symbol
rhoVals = [0 for i in range(m+1)]
GVals = [0 for i in range(m+1)]
GfundStorVals = [0 for i in range(m+1)]
GfundLossVals = [0 for i in range(m+1)]

zVals = [0 for i in range(m+1)]
kVals = [0 for i in range(m+1)]
impedVals = [0 for i in range(m+1)]

# filling QCM layer m values with AT cut quartz constants
rhoVals[m] = 2648 # (kg m^-3) density
GfundStorVals[m] = 2.947e10 # (Pa) storage shear modulus
GVals[m] = 2.947e10 # (Pa) shear modulus; effectively only has real part here
cVal = np.sqrt(GfundStorVals[m]/rhoVals[m]) # (m s^-1) speed of sound; ~3.34E3
LVal = cVal / (2*fundFreq) # (m) thickness of QCM; ~0.334E-3

# setting the values in the appropriate cells of the parameter arrays to what has been given as a guess for the top and bottom of the film; top here
rhoVals[0] = parameters['rhoTop'].value
GfundStorVals[0] = parameters['GfundStorTop'].value*1e6
GfundLossVals[0] = parameters['GfundLossTop'].value*1e6

# bottom layer values are entered here
if m > 1:
    rhoVals[m-1] = parameters['rhoBot'].value
    GfundStorVals[m-1] = parameters['GfundStorBot'].value*1e6
    GfundLossVals[m-1] = parameters['GfundLossBot'].value*1e6

# setting symbols to something readable in case you choose to display some of the math; you can use, for example, display(kVals[1])
for i in range(0, m+1):
    
    indexName = str(i)
    if i == m:
        indexName = 'm'
    
    # providing symbols at each layer
    rho[i] = sp.symbols('ρ_{%s}' % indexName)
    G[i] = sp.symbols('\widetilde{G}_{%s}' % indexName)
    GfundStor[i] = sp.symbols('G\'_{%s, n=1}' % indexName)
    GfundLoss[i] = sp.symbols('G\'\'_{%s, n=1}' % indexName)
    
    z[i] = sp.symbols('z_{%s}' % indexName)
    k[i] = sp.symbols('\widetilde{k}_{%s}' % indexName)
    imped[i] = sp.symbols('\widetilde{Z}_{%s}' % indexName)

    # setting the mathematical relations of symbols which produce the values of the wavenumber and impedance
    kVals[i] = omega*sp.sqrt(rho[i]/G[i])
    impedVals[i] = sp.sqrt(rho[i]*G[i])

# fills out zVals array; uses standard units (e.g. meters not nanometers)
def makeLayers(hTop, hBot):
    global botSideAmount # setting this to global so this function is able to reassign it if needed
    
    if m == 1:
        hBot = 0
    
    # top of QCM, bottom of film
    zVals[m] = hTop + hBot
    
    if m == 2:
        zVals[1] = hTop       
    
    # for the gradient models associated with m > 2, this is where layers are sliced up around the polymer-polymer interface according to botSideAmount for optimal resolution
    if m > 2:
        # imposing botSideAmount equals hTop when hTop varies is useful for preserving expected behavior
        if parameters['hTop'].vary:
            botSideAmount = hTop*1e6 # note the hTop entered here is in meters, so we set botSideAmount in micrometers to match the hTop parameter's units in the broader scope outside this function

        # there are m-1 film layers contained in the focused region; the thickness of the focused region is split up equally among these layers
        gradDivWidth = (hTop + (botSideAmount*1e-6)) / (m-1)
        # z[0] already equals 0 from the initialization; each subsequent layer of the focused region is one gradDivWidth deeper into the film
        for i in range(1, m):
            zVals[i] = zVals[i-1] + gradDivWidth

# this is where the position of the middle of a layer is assigned a value depending on the functional form of the gradient; we use a tanh form; uses standard units (e.g. meters not nanometers)
# this is like taking a median of the layer's region; a mean may be more physical and could be an interesting change to try, though the difference in the result is unlikely to be significant when m is reasonably high
def gradientModel(hTop, topVal, botVal, w, mu, pos):
    # tanh curve; subtracting hTop from z centers the curve by default on the compositional interface
    
    # a check that modulus is currently the property being handled and not density
    if not topVal == parameters['rhoTop'].value:
               
        # this handles cases where the code passes just the storage or loss component rather than the full complex value
        if isinstance(topVal, float):
            
            # conversion to a log scale for implementation in the hyperbolic tangent function
            topVal = np.log10(topVal)
            botVal = np.log10(botVal)
            
            # function outputs are converted back to a linear scale before being returned for processing
            return np.power(10, ((botVal + topVal)/2) + (((botVal - topVal)/2) * np.tanh(2*(pos - hTop - mu) / w)))
            
        # full complex modulus
        else:
          
            # conversion to a log scale for implementation in the hyperbolic tangent function
            topValRe = np.log10(float(sp.re(topVal)))
            topValIm = np.log10(float(sp.im(topVal)))
            botValRe = np.log10(float(sp.re(botVal)))
            botValIm = np.log10(float(sp.im(botVal)))
                       
            # function outputs are converted back to a linear scale before being returned for processing
            return np.power(10, ((botValRe + topValRe)/2) + (((botValRe - topValRe)/2) * np.tanh(2*(pos - hTop - mu) / w))) + 1j*np.power(10, ((botValIm + topValIm)/2) + (((botValIm - topValIm)/2) * np.tanh(2*(pos - hTop - mu) / w)))
                                         
    # density is simply handled on a linear scale, and is a sharp step change in our work anyway
    else:
        return ((botVal + topVal)/2) + (((botVal - topVal)/2) * np.tanh(2*(pos - hTop - mu) / w))

# fills out a parameter array using the gradient model when m > 2; uses standard units (e.g. meters not nanometers)
def gradFill(hTop, w, mu, valArray):   

    # parameter value at the top of the film
    topVal = valArray[0]
    # parameter value at the bottom of the film
    botVal = valArray[m-1]
          
    # note that this can overwrite valArray[0] in cases of strong gradients and/or low m
    # this is expected behavior, as you would not want to impose an unphysical abrupt localized deviation from the model curve in the softer top layer which is more sensitive to precise values
    for i in range(0, m-1):
        # position at center of layer
        centerPos = (zVals[i] + zVals[i+1]) / 2
        # using gradientModel to assign a value to the paramter according to the position at the center of the layer
        valArray[i] = gradientModel(hTop, topVal, botVal, w, mu, centerPos)
   
# these three funcions are used to scale axes in the graphs 
# get the maximum of a value array excluding the QCM's value
def filmMax(valArray):
    maxVal = valArray[0]
    for i in range(1, m):
        if valArray[i] > maxVal:
            maxVal = valArray[i]
    return maxVal

# get the minimum of a value array excluding the QCM's value
def filmMin(valArray):
    minVal = valArray[0]
    for i in range(1, m):
        if valArray[i] < minVal:
            minVal = valArray[i]
    return minVal

# get the span of a value array excluding the QCM's value
def span(valArray):
    # when m is 1, just returns a simple alternative because there is no meaningful span
    if m == 1:
        if not valArray[0] == 0:
            return valArray[0]
        else:
            return 1
        
    return filmMax(valArray) - filmMin(valArray)

layerNames = []
for i in range(m + 1):
    if i == m:
        layerNames.append('layer '+str(m)+' (QCM)')
    else:
        layerNames.append('layer '+str(i))

# function to display information about the model
def displayModel(hTop = parameters['hTop'].value, rhoTop = parameters['rhoTop'].value, GfundStorTop = parameters['GfundStorTop'].value, GfundLossTop = parameters['GfundLossTop'].value,
                 betaStorTop = parameters['betaStorTop'].value, betaLossTop = parameters['betaLossTop'].value,
                 hBot = parameters['hBot'].value, rhoBot = parameters['rhoBot'].value, GfundStorBot = parameters['GfundStorBot'].value, GfundLossBot = parameters['GfundLossBot'].value,
                 betaStorBot = parameters['betaStorBot'].value, betaLossBot = parameters['betaLossBot'].value,
                 wG = parameters['wG'].value, muG = parameters['muG'].value, wRho = parameters['wRho'].value, muRho = parameters['muRho'].value):
    
    # we run makeLayers, as well as gradFill on density and the storage and loss shear modulus for the purpose of displaying the input model
    makeLayers(hTop*1e-6, hBot*1e-6)
    if m > 2:
        gradFill(hTop*1e-6, wRho*1e-9, muRho*1e-9, rhoVals)
        gradFill(hTop*1e-6, wG*1e-9, muG*1e-9, GfundStorVals)
        gradFill(hTop*1e-6, wG*1e-9, muG*1e-9, GfundLossVals)
      
    # used to display the values at each layer in a table
    
    parameterMatrix = []
    parameterMatrix.append(np.multiply(zVals, 1e9))
    parameterMatrix.append(rhoVals)
    parameterMatrix.append(np.multiply(GfundStorVals, 1e-6))
    parameterMatrix.append(np.multiply(GfundLossVals, 1e-6))
    
    print(pd.DataFrame(np.transpose(parameterMatrix), columns=['z (nm)', 'ρ (kg m^-3)', 'G_f \' (MPa)', 'G_f \'\' (MPa)'], index = layerNames))
    
    if m > 2:
        print('\nG width: ' + str(wG) + ' nm')
        print('G asymmetry: ' + str(muG) + ' nm')
        print('ρ width: ' + str(wRho) + ' nm')
        print('ρ asymmetry: ' + str(muRho) + ' nm')
            
    # these arrays will contain information about the layers in nanometers for the purpose of plotting
    midPoints = []
    layerWidths = []
    
    for i in range(0, m):
        midPoints.append(((zVals[i] + zVals[i+1]) / 2)*1e9)
        layerWidths.append((zVals[i+1] - zVals[i])*1e9)
        
    # QCM layers, though the plots will not display them anyway
    midPoints.append((zVals[m] + LVal/2)*1e9)
    layerWidths.append(LVal)
    
    # a set of many positions in nanometers for plotting the full model curve on top of the corresponding layers
    zModelVals = np.linspace(0, (hTop + hBot)*1e3, 1000)
    
    # the figure that will contain the layer model plots
    fig = plt.figure(figsize = (8, 6))
    
    # interface position as data points for drawing interface line
    interfaceX = [hTop*1e3, hTop*1e3]
        
    # plotting density layer model
    ax1 = fig.add_subplot(221)
    
    # vertical bounds for interface line
    interfaceY = [filmMin(rhoVals), filmMax(rhoVals)]
    
    # model curve
    modelVals = []
    for i in range(len(zModelVals)):
        modelVals.append(gradientModel(hTop*1e-6, rhoTop, rhoBot, wRho*1e-9, muRho*1e-9, zModelVals[i]*1e-9))
    if m > 2:
        ax1.plot(zModelVals, modelVals, color = 'blue', label = 'model equation')
        
    # layers
    ax1.errorbar(midPoints, rhoVals, xerr = np.divide(layerWidths, 2), fmt = 'none', ecolor = 'blue', elinewidth = 3, capsize = 6, capthick = 4, label = 'model layers')
    # interface line
    if m > 1:
        ax1.plot(interfaceX, interfaceY, linestyle = 'dashed', color = 'black', label = 'interface')
        
    ax1.set(xlabel = 'Distance from Sample Top (nm)', ylabel = 'ρ (kg/(m^3)', xlim = (0, (zVals[m])*1e9), ylim = (filmMin(rhoVals) - (0.1 * span(rhoVals)), filmMax(rhoVals) + (0.1 * span(rhoVals))))
    
    # plotting the storage modulus layer model
    ax2 = fig.add_subplot(222)
    
    # vertical bounds for interface line
    interfaceY = [filmMin(GfundStorVals)*1e-6, filmMax(GfundStorVals)*1e-6]
    
    # model curve
    modelVals = []
    for i in range(len(zModelVals)):
        modelVals.append(gradientModel(hTop*1e-6, GfundStorTop*1e6, GfundStorBot*1e6, wG*1e-9, muG*1e-9, zModelVals[i]*1e-9))
    if m > 2:
        ax2.plot(zModelVals, np.multiply(modelVals, 1e-6), color = 'blue')
        
    # layers
    ax2.errorbar(midPoints, np.multiply(GfundStorVals, 1e-6), xerr = np.divide(layerWidths, 2), fmt = 'none', ecolor = 'blue', elinewidth = 3, capsize = 6, capthick = 4)
    # interface line
    if m > 1:
        ax2.plot(interfaceX, interfaceY, linestyle = 'dashed', color = 'black')
        
    ax2.set(xlabel = 'Distance from Sample Top (nm)', ylabel = 'G_f \' (MPa)', xlim = (0, (zVals[m])*1e9),
            ylim = ((filmMin(GfundStorVals) - (0.1 * span(GfundStorVals))) * 1e-6, (filmMax(GfundStorVals) + (0.1 * span(GfundStorVals))) * 1e-6))
    
    # plotting the loss modulus layer model
    ax3 = fig.add_subplot(223)
    
    # vertical bounds for interface line
    interfaceY = [filmMin(GfundLossVals)*1e-6, filmMax(GfundLossVals)*1e-6]
    
    # model curve
    modelVals = []
    for i in range(len(zModelVals)):
        modelVals.append(gradientModel(hTop*1e-6, GfundLossTop*1e6, GfundLossBot*1e6, wG*1e-9, muG*1e-9, zModelVals[i]*1e-9))
    if m > 2:
        ax3.plot(zModelVals, np.multiply(modelVals, 1e-6), color = 'blue')
    
    # layers
    ax3.errorbar(midPoints, np.multiply(GfundLossVals, 1e-6), xerr = np.divide(layerWidths, 2), fmt = 'none', ecolor = 'blue', elinewidth = 3, capsize = 6, capthick = 4)
    # interface line
    if m > 1:
        ax3.plot(interfaceX, interfaceY, linestyle = 'dashed', color = 'black')
        
    ax3.set(xlabel = 'Distance from Sample Top (nm)', ylabel = 'G_f \'\' (MPa)', xlim = (0, (zVals[m])*1e9),
            ylim = ((filmMin(GfundLossVals) - (0.1 * span(GfundLossVals))) * 1e-6, (filmMax(GfundLossVals) + (0.1 * span(GfundLossVals))) * 1e-6))
    
    # plot formatting and displaying
    handles, labels = ax1.get_legend_handles_labels()
    plt.suptitle('Layer Models Using the Input Parameter Values')
    plt.figlegend(handles, labels, loc = 'lower right', fontsize = 15)
    plt.tight_layout()
    plt.show();

# the core of the math takes place here; ultimately a Newton root finding algorithm outputs a list of complex resonance shifts
def calcResonance(hTop = parameters['hTop'].value, rhoTop = parameters['rhoTop'].value, GfundStorTop = parameters['GfundStorTop'].value, GfundLossTop = parameters['GfundLossTop'].value,
                  betaStorTop = parameters['betaStorTop'].value, betaLossTop = parameters['betaLossTop'].value,
                  hBot = parameters['hBot'].value, rhoBot = parameters['rhoBot'].value, GfundStorBot = parameters['GfundStorBot'].value, GfundLossBot = parameters['GfundLossBot'].value,
                  betaStorBot = parameters['betaStorBot'].value, betaLossBot = parameters['betaLossBot'].value,
                  wG = parameters['wG'].value, muG = parameters['muG'].value, wRho = parameters['wRho'].value, muRho = parameters['muRho'].value):
    
    resShifts = []   
    
    # resetting top and bottom layer density and modulus to the function inputs for case of fitting (gradFill needs these assigned correctly ahead); changes them to standard units first
    
    # top layer
    rhoVals[0] = rhoTop
    GfundStorVals[0] = GfundStorTop*1e6
    GfundLossVals[0] = GfundLossTop*1e6
    
    # bottom layer
    if m > 1:
        rhoVals[m-1] = rhoBot
        GfundStorVals[m-1] = GfundStorBot*1e6
        GfundLossVals[m-1] = GfundLossBot*1e6
    
    # remake the layers in case of iterative calls to this function due to fitting
    makeLayers(hTop*1e-6, hBot*1e-6)
    # assigning density values to remaining layers
    if m >  2:
        gradFill(hTop*1e-6, wRho*1e-9, muRho*1e-9, rhoVals)
        
    # initialize with the identity matrix so that when you first multiply prodMatrix by another matrix, it just gives that same matrix back
    prodMatrix = sp.Matrix([[1, 0], [0, 1]])

    # this for-loop executes the product series from the associated written work; iterates from 1 (inclusive) to m+1 (non-inclusive), carrying out all of the necessary matrix inversions and multiplications
    # it comes from continuity of stress across the internal interfaces between the modeled layers, as well as continuity of displacement imposed by the reasonable assumption of no slip between the modeled layers
    for i in range(1, m+1):
        
        # matrix associated with the modeled layer on the top side of z[i]
        topMatrix = sp.Matrix([[sp.exp(-sp.I*k[i-1]*z[i]), sp.exp(sp.I*k[i-1]*z[i])], [-imped[i-1]*sp.exp(-sp.I*k[i-1]*z[i]), imped[i-1]*sp.exp(sp.I*k[i-1]*z[i])]])
        # the symbols are substituted for values where it's currently appropriate; note the wavenumber and impedance values are functions of other symbols, which will need to be subbed for numerical values
        topMatrix = topMatrix.subs(z[i], zVals[i]).subs(k[i-1], kVals[i-1]).subs(imped[i-1], impedVals[i-1]).subs(rho[i-1], rhoVals[i-1])
        
        # matrix associated with the modeled layer on the bottom side of z[i]
        botMatrix = sp.Matrix([[sp.exp(-sp.I*k[i]*z[i]), sp.exp(sp.I*k[i]*z[i])], [-imped[i]*sp.exp(-sp.I*k[i]*z[i]), imped[i]*sp.exp(sp.I*k[i]*z[i])]])
        # similar substitutions to topMatrix
        botMatrix = botMatrix.subs(z[i], zVals[i]).subs(k[i], kVals[i]).subs(imped[i], impedVals[i]).subs(rho[i], rhoVals[i])
        
        # each iteration multiplies in another set of the bottom matrix and inverse of the top matrix
        prodMatrix = prodMatrix * (topMatrix.inv()*botMatrix)

    # the expression that relates terms in the needed way so that a root-finding algorithm can obtain Δf and ΔΓ from it; see the associated written work for details on how we arrived at this
    # it should equal zero when there is correctly continuity of stress at the bottom air interface
    expr = -sp.exp(-2*sp.I*k[m]*(z[m]+LVal)) + (prodMatrix.row(1)[0] - prodMatrix.row(0)[0])/(prodMatrix.row(0)[1] - prodMatrix.row(1)[1])

    # substitutions are made in the above expression
    expr = expr.subs(k[m], kVals[m]).subs(rho[m], rhoVals[m]).subs(z[m], zVals[m])

    # the root-finding algorithm seems pretty robust with regard to the sensitivity to this guess for the complex resonance shift; -10000 Hz has pretty much always worked
    resShiftGuess = -10000  
    
    # harmonic-specific substitutions and the root-finding are performed here
    for i in range(0, numHarmonics):
        
        subbed = expr
        
        GVals = [0 for i in range(m+1)]
        
        # top and bottom film layer complex modulus is calculated
        GVals[0] = complex(GfundStorTop*1e6*np.power(harmonics[i], betaStorTop) + 1j*(GfundLossTop*1e6*np.power(harmonics[i], betaLossTop)))
        if m > 1:
            GVals[m-1] = complex(GfundStorBot*1e6*np.power(harmonics[i], betaStorBot) + 1j*(GfundLossBot*1e6*np.power(harmonics[i], betaLossBot)))
        
        # complex modulus values are assigned across the system according to the gradient model
        if m > 2:
            gradFill(hTop*1e-6, wG*1e-9, muG*1e-9, GVals)
        
        # QCM value
        GVals[m] = complex(GfundStorVals[m])
        
        # for clarity/analysis, provides the 5 MHz storage and loss components of the shear modulus at each layer
        for j in range(0, m+1):
            if i == 0:
                GfundStorVals[j] = float(sp.re(GVals[j]))
                GfundLossVals[j] = float(sp.im(GVals[j]))
                
            # substituting the complex shear modulus values in for the symbols
            subbed = subbed.subs(G[j], GVals[j])
        
        # substituting the harmonic and expanding the form of the complex angular frequency
        subbed = subbed.subs(omega, omegaVal).subs(n, harmonics[i])
        # substitutes the provided initial values if you have chosen to analyze changes from an initial state
        if analyzeChange:
            subbed = subbed.subs(freqShiftInit, freqShiftInitList[i]).subs(dissipShiftInit, dissipShiftInitList[i])
        
        # transforms the math to a form that the root-finding algorithm can understand
        lambdified  = sp.lambdify((resShift), subbed, modules = ['mpmath'])
        
        # Newton root-finding algorithm finds the resonance shifts needed to make the expression equal 0
        result = mp.findroot(lambdified, (resShiftGuess), solver = 'newton')
        resShifts.append(complex(result))
    
    return resShifts

# a function that takes a set of parameters and returns a list of all of the residuals associated with each Δf and ΔΓ; this is called repeatedly by the fitting algorithm as it varies parameter values
def residuals(paramGuesses):
    
    paramGuessVals = list(paramGuesses.valuesdict().values())
          
    # resonance shifts are modeled
    modelVals = calcResonance(*paramGuessVals)
    
    # residuals are calculated and added to a list; the fitting algorithm expects the residuals to already be normalized by uncertainty if you have known uncertainties to include, hence why this is where it is done
    residualList = []
    for i in range(0, numHarmonics):
        residualList.append(float((sp.re(modelVals[i]) - freqShiftList[i]) / freqShiftUncList[i]))
        residualList.append(float((sp.im(modelVals[i]) - dissipShiftList[i]) / dissipShiftUncList[i]))
    
    return residualList

# calculates χ^2 error given modeled values for Δf and ΔΓ
def chiSquared(modelVals):

    totalError = 0
    for i in range(0, numHarmonics):
            error = ((sp.re(modelVals[i]) - freqShiftList[i])**2 / freqShiftUncList[i]**2) + ((sp.im(modelVals[i]) - dissipShiftList[i])**2 / dissipShiftUncList[i]**2)
            totalError = totalError + error

    return float(totalError)

# a list of the input parameter value guesses
paramVals = list(parameters.valuesdict().values())
paramKeys = list(parameters.valuesdict().keys())

harmonicNames = []
for i in range(numHarmonics):
    harmonicNames.append('n = ' + str(harmonics[i]))

# only displays and analyzes input guess if fitting is not set to be performed; this is just to save runtime and avoid screen clutter
if not performFitting:
    
    # running the displayModel function with all the default guess values set
    print('Layer attributes according to the guessed parameter values: \n')
    displayModel()
    
    print('\nRunning root-finding algorithm...')
    
    # calculate the resonance shifts given the input parameter values
    resonances = calcResonance(*paramVals)
    
    # initialize data matrix which will hold the solved Δf and ΔΓ values
    modelResult = [[0 for j in range(numHarmonics)] for i in range(2)]
    for i in range(0, numHarmonics):
        modelResult[0][i] = float(sp.re(resonances[i]))
        modelResult[1][i] = float(sp.im(resonances[i]))
    
    print('\nResult of parameter guess:\n')
        
    if analyzeChange:
        print(pd.DataFrame(np.transpose(modelResult), columns=['Δf - Δf_i (Hz)', 'ΔΓ - ΔΓ_i (Hz)'], index = harmonicNames))
    else:
        print(pd.DataFrame(np.transpose(modelResult), columns=['Δf (Hz)', 'ΔΓ (Hz)'], index = harmonicNames))
    
    # here we plot the modeled resonance shifts on top of the experimental values
    
    # experimental Δf data
    plt.errorbar(harmonics, freqShiftList, yerr = freqShiftUncList, fmt = 'o', ms = 9, mfc = 'white', markeredgewidth = 2.5, mec = 'black', ecolor = 'black', elinewidth = 2.5, capsize = 10, capthick = 2.5)
    # this is just here so that the marker in the legend doesn't have error bars    
    plt.scatter(harmonics, freqShiftList, label='Data', marker = 'o', s = 90, linewidths = 2.5, facecolors = 'white', edgecolors = 'black')
    
    # model result
    plt.plot(harmonics, modelResult[0], color = 'blue', alpha = 0.4, zorder = 10, linewidth = 2)
    plt.scatter(harmonics, modelResult[0], label='Model', marker = 'o', s = 30, color = 'blue', zorder = 10)
    
    #  formatting and displaying the Δf graph
    plt.legend(prop={'size': 10})
    plt.xlabel('Harmonic')
    if analyzeChange:
        plt.ylabel('Δf - Δf_i (Hz)')
    else:
        plt.ylabel('Δf (Hz)')
    plt.show();
    
    # experimental ΔΓ data
    plt.errorbar(harmonics, dissipShiftList, yerr = dissipShiftUncList, fmt = 'o', ms = 9, mfc = 'white', markeredgewidth = 2.5, mec = 'black', ecolor = 'black', elinewidth = 2.5, capsize = 10, capthick = 2.5)
    # this is just here so that the marker in the legend doesn't have error bars    
    plt.scatter(harmonics, dissipShiftList, label='Data', marker = 'o', s = 90, linewidths = 2.5, facecolors = 'white', edgecolors = 'black')
    
    # model result
    plt.plot(harmonics, modelResult[1], color = 'blue', alpha = 0.4, zorder = 10, linewidth = 2)
    plt.scatter(harmonics, modelResult[1], label='Model', marker = 'o', s = 30, color = 'blue', zorder = 10)
    
    # formatting and displaying the ΔΓ graph
    plt.legend(prop={'size': 10})
    plt.xlabel('Harmonic')
    if analyzeChange:
        plt.ylabel('ΔΓ - ΔΓ_i (Hz)')
    else:
        plt.ylabel('ΔΓ (Hz)')
        
    plt.show();
    
    # the χ^2 error for these Δf and ΔΓ values
    print('\nχ^2 error = %.4f' % chiSquared(resonances))

# if performFitting
else:
    
    print('Running fitting algorithm...\n')
    print('Parameters being varied:')

    for i in range(len(paramKeys)):
        if parameters[paramKeys[i]].vary:
            print(paramKeys[i])
    print()

    # this is where the fitting algorithm is run; it minimizes the sum of the residuals squared, which is χ^2 if you have your residuals normalized by uncertainty
    fitResult = minimize(residuals, parameters)
    # the parameter values obtained from fitting
    fitResultVals = list(fitResult.params.valuesdict().values())
   
    print('\nFit values:')
    
    # printing the result of the varied parameters
    for i in range(len(paramKeys)):
        if parameters[paramKeys[i]].vary:
            print(paramKeys[i] + ': ' + '%.4f' % fitResultVals[i])
    print()
            
    # calculating Δf and ΔΓ for the fit parameter values
    resonances = calcResonance(*fitResultVals)
  
    # displaying the best fit model
    print('Layer attributes according to the best fit parameter values: \n')
    displayModel(*fitResultVals)
    
    # initializing the data matrix which will hold the fit Δf and ΔΓ values
    fitShifts = [[0 for j in range(numHarmonics)] for i in range(2)]
    for i in range(0, numHarmonics):
        fitShifts[0][i] = float(sp.re(resonances[i]))
        fitShifts[1][i] = float(sp.im(resonances[i]))
    
    # displaying the results as a table
    print('\nΔf and ΔΓ values resulting from the fit:\n')
    if analyzeChange:
        print(pd.DataFrame(np.transpose(fitShifts), columns=['Δf - Δf_i (Hz)', 'ΔΓ - ΔΓ_i (Hz)'], index = harmonicNames))
    else:
        print(pd.DataFrame(np.transpose(fitShifts), columns=['Δf (Hz)', 'ΔΓ (Hz)'], index = harmonicNames))
       
    # plot of resonance shifts
    plt.rcParams['figure.figsize'] = (7, 5)

    # plotting experimental Δf data
    plt.errorbar(harmonics, freqShiftList, yerr = freqShiftUncList, fmt = 'o', ms = 9, mfc = 'white', markeredgewidth = 2.5, mec = 'black', ecolor = 'black', elinewidth = 2.5, capsize = 10, capthick = 2.5)
    # this is just here so that the marker in the legend doesn't have error bars    
    plt.scatter(harmonics, freqShiftList, label='Data', marker = 'o', s = 90, linewidths = 2.5, facecolors = 'white', edgecolors = 'black')

    # fit result
    plt.plot(harmonics, fitShifts[0], color = 'blue', alpha = 0.4, zorder = 10, linewidth = 2)
    plt.scatter(harmonics, fitShifts[0], label='Model', marker = 'o', s = 30, color = 'blue', zorder = 10)

    # formatting and displaying graph
    plt.legend(prop={'size': 10})
    plt.xlabel('Harmonic')
    if analyzeChange:
        plt.ylabel('Δf - Δf_i (Hz)')
    else:
        plt.ylabel('Δf (Hz)')
    plt.show();

    # plotting experimental ΔΓ data
    plt.errorbar(harmonics, dissipShiftList, yerr = dissipShiftUncList, fmt = 'o', ms = 9, mfc = 'white', markeredgewidth = 2.5, mec = 'black', ecolor = 'black', elinewidth = 2.5, capsize = 10, capthick = 2.5)
    # this is just here so that the marker in the legend doesn't have error bars    
    plt.scatter(harmonics, dissipShiftList, label='Data', marker = 'o', s = 90, linewidths = 2.5, facecolors = 'white', edgecolors = 'black')

    # fit result
    plt.plot(harmonics, fitShifts[1], color = 'blue', alpha = 0.4, zorder = 10, linewidth = 2)
    plt.scatter(harmonics, fitShifts[1], label='Model', marker = 'o', s = 30, color = 'blue', zorder = 10)

    # formatting and displaying graph
    plt.legend(prop={'size': 10})
    plt.xlabel('Harmonic')
    if analyzeChange:
        plt.ylabel('ΔΓ - ΔΓ_i (Hz)')
    else:
        plt.ylabel('ΔΓ (Hz)')
    plt.show();
    
    # the χ^2 error of these fit Δf and ΔΓ values
    print('\nχ^2 error = %.4f' % chiSquared(resonances))
