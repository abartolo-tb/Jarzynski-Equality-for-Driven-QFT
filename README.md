# Jarzynski-Equality-for-Driven-QFT
Numerical simulation code for the paper arXiv:1710.00829

The Mathematica "Generator-NtoN*.m" files are scripts meant to be run remotely. These scripts calculate the work distribution functions for the various scattering processes. Note, even with parallel processing, these scripts can take several days to run. These scripts were tailored for running on an 8-core machine and will require modification before being run on a machine with a different core count.

The Mathematica "NtoN*-workDist.m" files are consolidated versions of the work distribution functions calculated by the respective "Generator-NtoN*.m" file.

The Mathematica "Work-Analysis.nb" file is meant to be used interactively to investigate and plot the work distribution functions. This notebook has hard-coded numerical values for the work distribution functions taken from the .m files. As such, the "Generator-NtoN*.m" files are only needed if one desires to recalculate the work distribution functions at different points or at higher accuracy.
