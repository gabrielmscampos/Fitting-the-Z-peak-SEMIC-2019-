# Fitting the Z peak SEMIC 2019
In this repository is available the code and references used in this study.

Pre-filtered Run2011A data can be found at https://github.com/cms-opendata-education/zboson-exercise or at the data folder of this repository. All fitting tools were used from LMFIT (https://lmfit.github.io/lmfit-py/) library.

All mathematical functions and probability density functions and the convoluted integral are hardcoded in the script.

In a nutshell the best mathematical model that fit the Z⁰ resonance peak was given by the equation

<p align="center">
  <img src="http://latex.codecogs.com/gif.latex?%5Cpsi%28x%3B%5CGamma%2CM%2C%20%5Calpha%2C%20n%2C%20%5Csigma%2CA%29%20%3D%20Be%5E%7B-%5Cfrac%7Bx%7D%7B%5Ctau%7D%7D&plus;%20%5Cint_%7B-%5Cinfty%7D%5E%7B%5Cinfty%7DBW%28u%3B%5CGamma%2CM%2CA%29%5Ccdot%20CB%28x-u%2C%20%5Calpha%2C%20n%2C%20%5Csigma%2C%20M%2C%20A%29du." />
</p>

This model produced the following fit along with the selection criteria in transverse momentum (pT) and pseudorapidity (eta): 

<p align="center">
  <img src="https://github.com/gabrielmscampos/Fitting-the-Z-peak-SEMIC-2019-/blob/master/plots/model4:convoluted_breitwigner_crystalball_exponential.png" width="350" height="350" />
</p>

References:

[1] MOREIRA, M. A. The Standard Model of Particle Physics. Rev. Bras. Ensino Fís. 2009, vol. 31, n.1, pp.1306.1-1306.11.
[2] CMS Colaboration, The CMS Experiment at the Cern LHC, JINST 3  (2008) S08004.
[3] CMS Open Data, Z Boson exercise. (2018). Disponível em: https://github.com/cms-opendata-education/zboson-exercise.
[4] CAMPOS, G. M. S. Fitting the Z peak with CMS Open Data. (2019). DIsponível em: https://gist.github.com/gabrielmscampos/a528c942a039575000b4daae9da248c3.
[5] WIKIPEDIA. Crystal Ball.  Disponível em: en.wikipedia.org/wiki/Crystall_Ball_function.
[6] M. Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018).
[7]  Albert M. Sirunyan et al. Measurements of the inclusive W and Z production cross sections in pp collisions at sqrt(s) = 7 TeV with the CMS experiment.  JHEP10 (2011) 132.
