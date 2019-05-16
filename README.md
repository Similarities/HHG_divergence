# HHG_divergence
spectral divergence for High Harmonic spectrum, multi-lineouts with nonlinear dependency (e.g ~reciprocal), FWHM fit tool

High harmonic radiation is for instance emitted when realitvistic laser pulses interact with a solid surface. 
It has to be understood as an relativistic phenomena: relativistic electron oscillations reflect the incident 
laser pulse and introduce a relativistic Dopplershift to the reflected light. Since the oscillation is coupled
to the laser pulse period, an atto second pulse train is emitted. In detail the mechanism is a bit more complex :)
The HHG spectrum consits of Harmonics of the fundamental wavelength (lamdaL/N), recorded with a spectrometer
that is detected as well defined pattern of lines, that more or less have a 1/x spacing. 
By divergence we understand the spatial collmination of this emitted light, that usually changes with the HHG
order N, leading to a spatial narrowing of the HHG lines on the detector. 

This programm is a simple tool to determine the HHG divergence via lineouts over each harmonic number N. 
In detail, each lineout takes a integration over a number of lines into account (ROI_y) and determines 
the FWHM via stepfunction. Grating equation has to be inserted into the code. 
is not yet compatible to Python 3.x, since some input parameters (e.g. ROI_y have to be converted to int, 
which python 3 interprets as float.) . Reads 16 bit tiff files.

Example pictures show its output.
