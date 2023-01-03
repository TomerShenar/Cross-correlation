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


The script can either read fits files or txt files (set by the user under "WhatToRead" variable)
If fits: all fits files in "PathToObservations" are read;  dates are read from header
If txt: program expects file called "phases.txt" in the directory in which the observations are stored ("PathToObservations") with the following format:
OBSNAME1 MJD1
OBSNAME2 MJD2
etc...

User-adjustable variables:
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
