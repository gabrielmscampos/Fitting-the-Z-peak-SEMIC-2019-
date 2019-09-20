# Fitting the Z peak SEMIC 2019
In this repository is available all the code used in the study and the references used.

Pre-filtered Run2011A data can be found at https://github.com/cms-opendata-education/zboson-exercise or at the data folder of this repository. All fitting tools were used from LMFIT (https://lmfit.github.io/lmfit-py/) library.

All mathematical functions and probability density functions were hardcoded in the script and the convoluted integral is implemented by numerical integration.

In a nutshell the best mathematical model found to fit the Z⁰ resonance peak was given by the equation

<img src="http://latex.codecogs.com/gif.latex?%5Cpsi%28x%3B%5CGamma%2CM%2C%20%5Calpha%2C%20n%2C%20%5Csigma%2CA%29%20%3D%20Be%5E%7B-%5Cfrac%7Bx%7D%7B%5Ctau%7D%7D&plus;%20%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7DBW%28u%3B%5CGamma%2CM%2CA%29%5Ccdot%20CB%28x-u%2C%20%5Calpha%2C%20n%2C%20%5Csigma%2C%20M%2C%20A%29du." />

This model produced the fallowing fit along the selection criteria in the pT and eta: 

![alt text](http://url/to/img.png)
