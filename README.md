##########################

1D Cross-correlation following Zucker+ 2003 MNRAS, 342, 1291

Input: path to folder containing normalized spectra in fits format

Input: a template (or mask) used to cross-correlate (e.g., observation, model...)

Output: RVs_CCF.txt: FILENAME || JD || RV || error || S2N

Output: co-added spectrum in frame-of-reference of star.

Method described in: 

Zucker & Mazeh 1994, ApJ, 420, 806

Zucker+ 2003 MNRAS, 342, 1291

Shenar+ 2019, A&A, 627A, 151

Developer: Tomer Shenar, T.Shenar@uva.nl

##########################

Important note:


The script can either read fits files or txt files (set by the user under "WhatToRead" variable)

If fits: all fits files in "PathToObservations" are read;  dates are read from header

If txt: program expects file called "phases.txt" in the directory in which the observations are stored ("PathToObservations") with the following format:

OBSNAME1 MJD1

OBSNAME2 MJD2

etc...


##########################

User-adjustable variables:

##########################

intr_kind (str; interpolation type: linear, cubic, etc... default cubic)

WhatToRead (str; FITS or TXT)

PlotFirst (bool; Plot the CCF + match between template and spectrum only for first spectrum... recommended TRUE)

PlotAll  (bool; Plot the above for all spectra.... recommended TRUE for careful check, false for faster run)

Fit_Range_in_fraction (float; threshold fraction of the CCF for which a parabola is fit [play with it and see impact])

CrossCorRangeA (array of tuples; range borders at which to compute the CCF [can be multiple ranges, e.g. lines]. User can introduce lines as variables if they please)

CrossVeloMin (float; minimum RV in km/s to be scanned for)

CrossVeloMax (float; maximum RV in km/s to be scanned for)

S2Nrange (tuple; where to compute S2N... only important for final stacking of spectra)

PathToObservations (str; relative/absolute path to where observations are stored)

PathToOutput (str; where output should be stored)

MaskPath (str; where template [or mask] is stored).

########

Example:

########

Download the available CCF.py script, the make_spectra_SB2.py script, and the two TLUSTY model atmospheres provided here (source: http://tlusty.oca.eu/)
The make_spectra_SB2.py script will create mock observations using these models, based on the parameters specified by the user. For a flux ratio of Q = 0.0, these will be SB1 observations. For the first run, no need to change anything.

After installing all necessary packages (via pip install)...

1. Run make_spectra_SB2.py 

- 10 mock observations will be created with S2N = 40, R = 20000, sampling of 2 pixels per resolution element for a binary comprising the two TLUSTY models with P = 473d, T0=32.38 [d], e=0, omega=90deg, K1= 100 km/s (and K2=50km/s, but for Q=0 this is not relevant).

- Additionally, the file "phases.txt" will be created, needed to run CCF.py

- Finally, the files Atemp.txt and Btemp.txt will be created. These are the input spectra used to create the mock data after convolution with the instrument (i.e. convolution with Gaussian and re-sampling).

2. Run CCF.py

The script should compute the CCFs and their maxima (by fitting parabolas) using the template and all observations.

The template used by default in this example is "Atempt.txt", i.e., the template used to compute the data. Feel free to experiment with other templates (e.g., one of the observations).

As an output, the user will find:

- RVs_CCF.txt : a file containing the RVs and the errors

- coadded.txt : the coadded spectrum, formed by shifting-and-adding all spectra in the frame-of-reference of the template.

Check how the measured RVs compare with the intrinsic RVs (coded in the file names).

Play with the parameters! Good luck!

