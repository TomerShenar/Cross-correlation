# 1D Cross-correlation following Zucker+ 2003 MNRAS, 342, 1291
# Input: path to folder containing normalized spectra in fits format
# Output: RVs_CCF.txt: FILENAME || JD || RV || error || S2N
# Output: co-added spectrum in frame-of-reference of star.
# Method described in: 
# Zucker & Mazeh 1994, ApJ, 420, 806
# Zucker+ 2003 MNRAS, 342, 1291
# Shenar+ 2019, A&A, 627A, 151
# Developer: Tomer Shenar, T.Shenar@uva.nl
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import sys
from scipy.interpolate import interp1d
import astropy.io.fits as fits
from astropy.io import ascii
import pandas as pd

clight = 2.9979E5
eps = 1E-10

# ***************************************
# ********* USER INPUT ******************
# ***************************************

# Determines type of interpolation when reading spectra: should not be changed
intr_kind = 'cubic'

# What kind of observations to read?  If FITS, then will try to read *fits file in given directory (below) 
# If TXT, then will look for "phases.txt" file which contains all the dates...

#WhatToRead = "FITS"
WhatToRead = "TXT"

# Determines whether diagnostic plots are shown (first / all)
PlotFirst = True
PlotAll = False

# Determines the range in which a parable is fit to the CCF function.
Fit_Range_in_fraction = 0.95

# Oversample wavelength grid (helpful for low-res spectra). Default: OverSampleFac = 1 (=No), else: set value (recommended ~5-10)
OverSampleFac = 1

CrossCorRangeA = [[4000., 4500.]]

# Minimum and maximum RVs to be searched
CrossVeloMin = -400.
CrossVeloMax = 400.


# Where to calculate S2N (not important for RV measurement)
S2Nrange = [4450., 4455.]


# Path to observations: Directory in which the spectra are found:
PathToObservations = './obs/'
# PathToObservations = './'

# Path to output: Directory in which output should be written:
PathToOutput = "./"

# Mask path
# MaskPath = (PathToObservations +
#             '00956419_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf_norm.fits')
MaskPath = PathToObservations + 'Atemp.txt'
#MaskPath = PathToOutput + 'coadded.txt'

# ***************************************
# ********* END USER INPUT **************
# ***************************************


# ***************************************
# ********* Functions *******************
# ***************************************

# Returns CCF of two functions f1(x) and f2(x) [f1 = observation, f2 = mask]
def CCF(f1, f2):
    f1 = f1
    f2 = f2
    return np.sum(f1 * f2) / np.std(f1) / np.std(f2) / N


# Returns RV and error following Zucker+ 2003
def crosscorreal(Observation, Mask):
    global PlotFirst
    CCFarr = np.array([CCF(np.copy(Observation),
                           (np.roll(Mask, s))[CrossCorInds]) for s in sRange])
    IndMax = np.argmax(CCFarr)
    vmax = veloRange[IndMax]
    CCFMAX1 = CCFarr[IndMax]
    LeftEdgeArr = np.abs(Fit_Range_in_fraction * CCFMAX1 - CCFarr[:IndMax])
    RightEdgeArr = np.abs(Fit_Range_in_fraction * CCFMAX1 - CCFarr[IndMax+1:])

    if len(LeftEdgeArr) == 0 or len(RightEdgeArr) == 0:
        print("Can't find local maximum in CCF\n")
        return np.array([np.nan, np.nan])

    IndFit1 = np.argmin(np.abs(Fit_Range_in_fraction * CCFMAX1 -
                               CCFarr[:IndMax]))
    IndFit2 = np.argmin(np.abs(Fit_Range_in_fraction * CCFMAX1 -
                               CCFarr[IndMax+1:])) + IndMax + 1
    a, b, c = np.polyfit(veloRange[IndFit1:IndFit2+1],
                         CCFarr[IndFit1:IndFit2+1], 2)
    vmax = -b/(2*a)
    CCFAtMax = min(1-1E-20, c - b**2/4./a)
    FineVeloGrid = np.arange(veloRange[IndFit1], veloRange[IndFit2], .1)
    parable = (a*FineVeloGrid**2 + b*FineVeloGrid + c)

    if PlotFirst is True or PlotAll is True:
        # plot the ccf
        fig1, ax1 = plt.subplots()
        ax1.plot(veloRange, CCFarr, color='C0')
        ax1.plot(FineVeloGrid, parable, color='C1', linewidth=1.5)
        ax1.set_xlabel('Radial velocity [km/s]')
        ax1.set_ylabel('Normalized CCF')
        try:
            cutname = observation.split(PathToObservations)[1].split('.fits')[0]
        except:
            cutname = PathToObservations
        if not os.path.isdir('CCFParabolas'):
            os.mkdir('CCFParabolas')
        #fig1.savefig('CCFParabolas/CCF_parabola_' + cutname + '.pdf')
        plt.show()
        # plot the spectrum and the template
        fig2, ax2 = plt.subplots()
        ax2.plot(wavegridlog[CrossCorInds], Observation, color='k',
                 label='observation', alpha=0.8)
        ax2.plot(wavegridlog, (Mask - np.mean(Mask)), color='orchid',
                 label='template, unshifted', alpha=0.9)
        ax2.plot((wavegridlog*(1+vmax/clight)), (Mask - np.mean(Mask)),
                 color='turquoise', label='shifted', alpha=0.9)
        ax2.set_xlabel(r'Wavelength [$\AA$]')
        ax2.set_ylabel('Normalized flux')
        ax2.legend(loc='best')
        plt.show()
        PlotFirst = False

    if CCFAtMax > 1:
        print("Failed to cross-correlate: template probably sucks!")
        print("Check cross-correlation function + parable fit.")
        return np.nan, np.nan

    CFFdvdvAtMax = 2*a
    return np.array([vmax, np.sqrt(-1./(NRes * CFFdvdvAtMax *
                                        CCFAtMax / (1 - CCFAtMax**2)))])


# Reads different types of files, returns 2D array of waves vs fluxes
def read_file(infile, col0=0, col1=1):
    ext = str(infile.split('.')[-1])

    # Check type of input file (fits or ascii) to read in data correctly
    if (ext == 'fits') or (ext ==  'fit'):
        wave, flux = read_fits(infile)

    elif (ext == 'gz'):
        wave, flux = read_tlusty(infile)

    elif (ext == 'dat' or ext == 'ascii' or ext == 'txt' or ext == 'nspec'):
        wave, flux = read_ascii(infile)

    elif (ext == 'tfits'):
        wave, flux = read_uvespop(infile)

    elif (ext == 'hfits'):
        wave, flux = read_hermes_normalized(infile)

    else:
        wave, flux = read_ascii(infile, col0, col1)
    return np.array([wave, flux]).T

# Reads ASCII file
def read_ascii(infile, col0=0, col1=1):
    # any type of ascii file (typically I call them .dat)
    # assumes that first column is wave and second column is flux
    print("%s: Input file is an ascii file." % infile)
    #wave, flux = np.loadtxt(infile, usecols=(0, 1), unpack=True)
    AsciiFile = (pd.read_csv(infile, header=None, delim_whitespace=True, comment='#')).values
    #bla = np.loadtxt(infile)
    ##print bla
    wave = AsciiFile[:,col0]
    flux = AsciiFile[:,col1]
    return wave, flux


#Reads fits file
def read_fits(infile):
    print("%s: Input file is a fits file." % infile)

    header = fits.getheader(infile)

    if 'HIERARCH SPECTRUM EXTRACTION' in header:
        wave, flux = read_psfSpec(infile)

    elif 'INSTRUME' in header:
        ins = header['INSTRUME']
        if (ins == 'MUSE'):
            wave, flux = read_pampelMUSE(infile)

        elif (ins == 'HERMES'):
            wave, flux = read_HERMES(infile)

        elif (ins == 'FEROS'):
            wave, flux = read_FEROS(infile)
        elif (ins == 'XSHOOTER'):
            wave, flux = read_XSHOOTER(infile)

        elif (ins == 'UVES'):
            wave, flux = read_UVES(infile)
        elif (ins == 'GIRAFFE' and 'nLR' in infile):
            wave, flux = read_GIRAFFE(infile)     
        elif (ins == 'GIRAFFE'):
            wave, flux = read_GIRAFFE(infile)               
            #wave, flux = read_GIRAFFE2(infile)               
        elif (ins == 'ESPCOUDE'):
            wave, flux = read_NLA(infile)   
        elif (ins == 'COS'):
            wave, flux = read_COS(infile)               
        elif (ins == 'STIS'):
            wave, flux = read_STIS(infile)              
        else:
            print('File type unkown, trying HERMES')
            wave, flux = read_HERMES(infile)                        
    else:
        wave, flux = read_HERMES(infile)

    return wave, flux


# Reads HERMES files
def read_HERMES(infile):
    # print("%s: Input file is a HERMES file." % infile)
    header = fits.getheader(infile)
    if "norm" in str(infile):
        if 'CRVAL1' in header:
            wave, flux = read_spec(infile)

    # for files with standard wavelegth array
    if ((header['CTYPE1'] == 'WAVELENGTH') or (header['CTYPE1'] == 'AWAV')):
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = crval + np.arange(0, naxis1) * cdelt

    # for files that are given in logarithmic wl array
    if (header['CTYPE1'] == 'log(wavelength)'):
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = np.exp(crval + np.arange(0, naxis1)*cdelt)

    else:
        print("Could not read in HERMES fits file - unknown file type.")
        sys.exit()

    return wave, flux

# Reads FEROS files
def read_FEROS(infile):
    print("%s: Input file is a FEROS file." % infile)

    if "norm" in str(infile):
        header = fits.getheader(infile)
        if 'CRVAL1' in header:
            wave, flux = read_spec(infile)
        else:
            wlmin = header['WAVELMIN']
            wlmax = header['WAVELMAX']
            spec_bin = header['SPEC_BIN']

            wave = np.arange(wlmin, wlmax+spec_bin, spec_bin) * 10  # AA
            flux = fits.getdata(infile)

    else:

        header = fits.getheader(infile)
        if 'CRVAL1' in header:
            print('here')
            flux = fits.getdata(infile)
            crval = header['CRVAL1']
            crpix = header['CRPIX1']
            cdelt = header['CDELT1']

            wave = crval - (cdelt * crpix -
                            cdelt) + np.arange(flux.shape[0]) * cdelt

        else:
            data = fits.getdata(infile)
            wave = data.field(0)[0]
            flux = data.field(1)[0]
    return wave, flux

# Reads fits files in form of standard spectrum
def read_spec(infile):
    header = fits.getheader(infile)
    data = fits.getdata(infile)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    wave = wl0 - (delt * pix - delt) + np.arange(data.shape[0]) * delt

    flux = data

    return wave, flux

#Retrieved key value from header
def get_key_from_header(infile, key):
    header = fits.getheader(infile)
    info = header[key]
    return info


# ***************************************
# ********* End of Functions ************
# ***************************************

for CCFplot in glob.glob('*CCF_parable_*'):
    os.remove(CCFplot)

# Reading all normalized spectra in specified folder

# File with obsname + HJD

if WhatToRead == "TXT":
    SpectraDatesPath = PathToObservations + 'phases.txt'
    Observations = ascii.read(SpectraDatesPath)['obsname']
    Dates = np.array(ascii.read(SpectraDatesPath)['HJD'])
    ObsArray = enumerate(Observations)
else:
    Observations = []
    Dates = []
    ObsArray = enumerate(glob.glob(PathToObservations +  '*rect.fits'))

# Read all observations:
for i, observation in ObsArray:
    obspath = str(observation)
    sys.stdout.write("trying to read observation " + str(i + 1) + " in "
                     + obspath + " ...")
    spectrum =read_file(obspath)
    if WhatToRead=='FITS':
        Observations.append(observation)
        date = get_key_from_header(obspath, 'MJD-OBS')
        Dates.append(date)
    S2N = (spectrum[:, 0] > S2Nrange[0]) * (spectrum[:, 0] < S2Nrange[1])
    print('Done')
    if 'Spectra' in locals():
        Spectra.append(spectrum)
    if 'S2Ns' in locals():
        S2Ns.append(1/(np.std(spectrum[:, 1][S2N])))
    else:
        S2Ns = [1/(np.std(spectrum[:, 1][S2N]))]
        Spectra = [spectrum]

print("Preparing wavelength and velocity grids...")

LambdaRangeUser = CrossCorRangeA * \
    np.array([1. - 1.1*CrossVeloMax/clight, 1 - 1.1*CrossVeloMin/clight])
LamRangeB = LambdaRangeUser[0, 0]
LamRangeR = LambdaRangeUser[-1, 1]
Dlam = Spectra[0][1, 0] - Spectra[0][0, 0]
Resolution = Spectra[0][1, 0]/Dlam

#For N in error formula (NRes):

vbin = clight / Resolution 

Nwaves = int(np.log(LamRangeR/LamRangeB) / np.log(1. + vbin/clight))
wavegridlog = LamRangeB*(1. + vbin/clight)**np.arange(Nwaves)


IntIs = np.array([np.argmin(np.abs(wavegridlog - CrossCorRangeA[i][0]))
                  for i in np.arange(len(CrossCorRangeA))])
IntFs = np.array([np.argmin(np.abs(wavegridlog - CrossCorRangeA[i][1]))
                  for i in np.arange(len(CrossCorRangeA))])

Ns = IntFs - IntIs
NRes = np.sum(Ns)


if OverSampleFac == 1:
# N = number of grid points
    N = NRes
else:
#If oversample: redefine grids...
    print("Oversampling grids by factor ", OverSampleFac)
    vbin = clight / Resolution / OverSampleFac

    Nwaves = int(np.log(LamRangeR/LamRangeB) / np.log(1. + vbin/clight))
    wavegridlog = LamRangeB*(1. + vbin/clight)**np.arange(Nwaves)


    IntIs = np.array([np.argmin(np.abs(wavegridlog - CrossCorRangeA[i][0]))
                    for i in np.arange(len(CrossCorRangeA))])
    IntFs = np.array([np.argmin(np.abs(wavegridlog - CrossCorRangeA[i][1]))
                    for i in np.arange(len(CrossCorRangeA))])

    Ns = IntFs - IntIs
    N = np.sum(Ns)


sRange = np.arange(int(CrossVeloMin/vbin), int(CrossVeloMax/vbin)+1, 1)
veloRange = vbin*sRange
CrossCorInds = np.concatenate(([np.arange(IntIs[i], IntFs[i])
                                for i in np.arange(len(IntFs))]))

Observations, Dates = np.array(Observations), np.array(Dates)
sort = np.argsort(Dates)


print("Reading masks...")

MaskAll = read_file(MaskPath)
Mask = interp1d(MaskAll[:,0], np.nan_to_num(MaskAll[:,1]), bounds_error=False,
                fill_value=1., kind=intr_kind)(wavegridlog)

print("Measuring RVs...")

Vs = []
sigs = []
Vrels = []
# Perform CCF for each observation
for i, spectrum in enumerate(Spectra):
    fluxes = interp1d(spectrum[:, 0], np.nan_to_num(
        spectrum[:, 1]), bounds_error=False, fill_value=1.,
        kind='cubic')(wavegridlog[CrossCorInds])
    CCFeval = crosscorreal(np.copy(fluxes - np.mean(fluxes)),
                           np.copy(Mask - np.mean(Mask)))
    sigs.append(CCFeval[1])
    Vs.append(CCFeval[0])
    print('RV for ' + str(Observations[i]) + ': ' + str(CCFeval[0]) + ' +/- ' +
          str(CCFeval[1]))


# Save RVs
if 'phases' in locals():
    np.savetxt(PathToOutput + '/RVs_CCF.txt',
               np.c_[Observations[sort], Dates[sort], phases[sort],
                     (np.array(Vs))[sort], np.array(sigs)[sort],
                     np.array(S2Ns)[sort]],
               header='PATH  MJD phase v1 sig S2N', fmt="%s")
else:
    np.savetxt(PathToOutput + '/RVs_CCF.txt',
               np.c_[Observations[sort], Dates[sort],
                     (np.array(Vs))[sort], np.array(sigs)[sort],
                     np.array(S2Ns)[sort]],
               header='PATH  MJD v1 sig S2N', fmt="%s")

# Make co-added spectrum
wavegrid = Spectra[0][:, 0]
S2NsqrdSum = np.sum(np.array(S2Ns)**2)
for i, spectrum in enumerate(np.array(Spectra)[sort]):
    v = (np.array(Vs)[sort])[i]
    weightS2N = (np.array(S2Ns)[sort])[i]**2/S2NsqrdSum
    shiftspec = interp1d(spectrum[:, 0]*(1. - v/clight),
                         np.nan_to_num(spectrum[:, 1]),
                         bounds_error=False, fill_value=1.,
                         kind='cubic')(wavegrid)
    if 'CoaddSpec' in locals():
        CoaddSpec += weightS2N*shiftspec
    else:
        CoaddSpec = weightS2N*shiftspec
np.savetxt('coadded.txt', np.c_[wavegrid, CoaddSpec])

# plot measured RVs
fig3, ax3 = plt.subplots()
ax3.errorbar(Dates[sort], (np.array(Vs))[sort],  yerr=np.array(sigs)[sort],
             fmt='o', color='midnightblue')
ax3.set_xlabel('MJD [d]')
ax3.set_ylabel(r'$\Delta$RV [km/s]')
plt.show()
