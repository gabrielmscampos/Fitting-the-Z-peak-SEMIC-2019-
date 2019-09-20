#########################
## MODULES
#########################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import scipy.integrate as integrate
import scipy.special as special

print("\nLoading modules: Ok!")

#########################
## SETUPING ENV
#########################

constant = {"z_mass"         : 91.1876,                   # z boson mass [Gev/c²]
            "error_z_mass"   : 0.0021,                    # z boson mass error [Gev/c²]
            "hbar"           : 6.62607015e-34/(2*np.pi),  # reduced Planck constant [J.s]
            "Az"             : 0.3978,                    # acceptance for the Z -> mu+mu- selection
            "error_Az"       : 0.0005,                    # acceptance error for the Z -> mu+mu- selection
            "eff_Zmumu"      : 0.87,                      # corrected selection efficiency for the Z -> mu+mu-
            "error_eff_Zmumu": 0.011,                     # corrected selection efficiency error for the Z -> mu+mu-
            "cs_Zmumu"       : 968,                       # cross section for inclusive Z -> mu+mu- production [pb]
            "error_cs_Zmumu" : 8}                         # cross section error for inclusive Z -> mu+mu- production [pb]

dpi = 300 # matplotlib savefig argument
style = dict(size = 11, color = "black")  # Styling text in matplotlib

print("Setuping constants: Ok!")

############################################
## Importing Run2011A di-muon's channel data
############################################

'''
Pre-filtered Run2011A data can be found at:
	https://github.com/cms-opendata-education/zboson-exercise
'''

filepath = "/home/gabriel/Documents/UERJ/iniciação-científica/z-tutorial/DoubleMuRun2011A.csv"
run2011A  = pd.read_csv(filepath, sep = ",")

print("Importing data: Ok!")

##############################
## Defining functions and pdfs
##############################

class Vector2:
    '''A class for computing 2-dimensional vectors operations '''
    
    def __init__(self, x, y):
        self.x = np.array(x)
        self.y = np.array(y)
        
    def __str__(self):
        return "({}, {})".format(str(self.x), str(self.y))
    
    def __add__(self, v):
        return Vector2(self.x + v.x, self.y + v.y)
    
    def __sub__(self, v):
        return Vector2(self.x - v.x, self.y - v.y)
    
    def __mul__(self, n):
        return Vector2(n*self.x, n*self.y)
    
    def __truediv__(self, n):
        return Vector2(self.x/n, self.y/n)
    
    def magnitude(self):
        return np.power(np.power(self.x, 2) + np.power(self.y, 2), 1/2)
    
    def dot_product(self, v):
        return self.x*v.x + self.y*v.y
	
def inv_mass(pt1, pt2, eta1, eta2, phi1, phi2):
    '''Compute the invariant mass'''
    
    # pt1: transverse momentum for the first particle
    # pt2: transverse momentum for the second particle
    # eta1: pseudorapidity for the first particle
    # eta2: pseudorapidity for the second particle
    # phi1: azimuthal angle for the second particle
    # phi2: azimuthal angle for the second particle
    return ( 2*pt1*pt2*(np.cosh(eta1 - eta2) - np.cos(phi1 - phi2)) )**(1/2)

def linear(E, a, b):
	'''Evaluate the Linear function'''
	return a*E + b

def breitwigner(E, M, gamma, A):
	'''Evaluate the Relativistic Breit-Wigner function'''
	return A*((2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))))/((E**2-M**2)**2+M**2*gamma**2)

def pseudo_breitgauss(E, gamma, M, A, alfa, sigma):
	'''Evaluate a pseudo convoluted relativistic Breit-Wigner -- Standard Gaussian Distribution PDF'''
	return alfa*A*( (2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))) )/((E**2-M**2)**2+M**2*gamma**2) + (1-alfa)*A*np.exp(-(E-M)**2/(2*sigma**2))/(sigma*(2*np.pi)**(1/2))

def gaussian(E, M, sigma, A):
	'''Evaluate the Standard Gaussian PDF'''
	return A*np.exp(-(E-M)**2/(2*sigma**2))/(sigma*(2*np.pi)**(1/2))

def breitgauss(E, M, gamma, sigma, A):
	'''Evaluate the convoluted relativistic Breit-Wigner -- Gaussian PDF'''

	np = 100    ## Number of convolution steps
	sc = 5.0    ## Convolution extends to +-sc Gaussian sigmas

	Elow = E - sc*sigma   ## Minimum range of convolution integral
	Eupp = E + sc*sigma   ## Maximum range of convolution integral

	step = (Eupp - Elow)/np   ## Step

	## Convolution integral by discretized sum
	summ = 0
	for i in range(1, 50):
		xx = Elow + (i-0.5)*step
		fbw = breitwigner(xx, M, gamma, A)
		summ = summ + fbw*gaussian(E, xx, sigma, A)

		xx = Eupp - (i-0.5)*step
		fbw = breitwigner(xx, M, gamma, A)
		summ = summ + fbw*gaussian(E, xx, sigma, A)

	return (step*summ)/sigma

def exponential(E, tau, B):
	'''Evaluate the Exponential function'''
	return B*np.exp(-E/tau)

def crystalball_function(x, alpha, n, sigma, mean):
	'''Evaluate the Crystal Ball function'''

	if sigma < 0.:
	    return 0.

	a = np.abs(alpha)
	z = (x-mean)/sigma

	if alpha < 0:
	    z = -z

	A = np.power(n/a, n)*np.exp(-0.5*np.power(a, 2))
	B = (n/a) - a

	condition = np.where(True, z > -alpha, z <= -alpha)
	evaluate = np.where(condition, np.exp(-0.5*np.power(z, 2)), A*np.power(B - z, -n))

	return evaluate

def crystalball_pdf(x, alpha, n, sigma, mean, A):
    '''Evaluation of the PDF (is defined only for n > 1)'''
    
    if sigma < 0.:
        return 0.
    
    if n <= 1:
        return np.NaN
    
    aa = np.abs(alpha)
    C  = (n/aa)*(1./(n-1.))*np.exp(-alpha*alpha/2.)
    D  = np.sqrt(np.pi/2.)*(1 + special.erf(aa/np.sqrt(2.)))
    N  = 1./(sigma*(C + D))
    
    return A*N*crystalball_function(x, alpha, n, sigma, mean)

def breitball(E, M, gamma, sigma, A, alpha, n):
	'''Evaluation of the Convoluted Relativistic Breit-Wigner with a Crystal Ball (is defined only for n > 1)'''
    
	np = 100    ## Number of convolution steps
	sc = 5.0    ## Convolution extends to +-sc Gaussian sigmas

	Elow = E - sc*sigma   ## Minimum range of convolution integral
	Eupp = E + sc*sigma   ## Maximum range of convolution integral

	step = (Eupp - Elow)/np   ## Step

	## Convolution integral by discretized sum
	summ = 0
	for i in range(1, 50):
	    xx = Elow + (i-0.5)*step
	    fbw = breitwigner(xx, M, gamma, A)
	    summ = summ + fbw*crystalball_pdf(x = E, alpha = alpha, n = n, sigma = sigma, mean = xx, A = A)
	    
	    xx = Eupp - (i-0.5)*step
	    fbw = breitwigner(xx, M, gamma, A)
	    summ = summ + fbw*crystalball_pdf(x = E, alpha = alpha, n = n, sigma = sigma, mean = xx, A = A)
	    
	return (step*summ)/sigma

print("Defining functions: Ok!")

####################################
## Computing the transverse momentum
####################################

PtVector = Vector2(run2011A['px1'] + run2011A['px2'], run2011A['py1'] + run2011A['py2'])
P1Vector = Vector2(run2011A['px1'], run2011A['py1'])
P2Vector = Vector2(run2011A['px2'], run2011A['py2'])

momentumA = pd.DataFrame( {'pt' : Vector2.magnitude(PtVector),
                           'pt1': Vector2.magnitude(P1Vector),
                           'pt2': Vector2.magnitude(P2Vector) } )

print("Computing pT: Ok!")

####################################
## Reconstructing di-muon's channel
####################################

spectrumA = pd.DataFrame( { 'dimu': inv_mass(momentumA['pt1'], momentumA['pt2'], run2011A['eta1'], run2011A['eta2'], run2011A['phi1'], run2011A['phi2']) } )

nbins    = 500
weights  = np.array([nbins/np.log(10)/inv_mass for inv_mass in spectrumA.dimu.values])
log_mass = np.log10(spectrumA.dimu.values)

plt.hist(log_mass, bins = nbins, range = (-0.5, 2.5), weights = weights, lw = 0, color = "darkgrey")
plt.yscale("log")
plt.xlim(-0.5, 2.5)
plt.xlabel("Log( Massa [GeV/c²] )", fontsize = '12')
plt.ylabel("Eventos", fontsize = '12')
plt.title(r"Espectro de massa do canal de decaimento do sistema $\mu+\mu-$")
plt.savefig("img/dimuon_mass_spectrum.png", dpi = dpi)
plt.show()

print("Reconstructing di-muon's channel: Ok!")

####################################
## Selecting Z Boson Resonance Peak
####################################

z_boson  = pd.DataFrame( { 'mass' : [ mass for mass in spectrumA['dimu'] if ((mass >= 60) and (mass <= 120)) ] } )

nbins               = 200
y_zboson, bin_edges = np.histogram(z_boson.mass.values, bins = nbins, range = (60,120))
x_zboson            = (bin_edges[:-1] + bin_edges[1:])/2.

plt.hist(z_boson.mass.values, bins = nbins, color = "blue", label = "Entries: {0:.0f} \nMean: {1:.2f} \nStd Dev: {2:.2f}".format(len(z_boson.mass.values), np.mean(z_boson.mass.values), np.std(z_boson.mass.values)))
plt.xlim(60, 120)
plt.ylim(0)
plt.xlabel("Massa [GeV/c²]", fontsize = '12')
plt.ylabel("Eventos", fontsize = '12')
plt.legend(loc = "upper right", handlelength = 0, handletextpad = 0)
plt.title(r"Pico de massa do bóson $Z \rightarrow \mu+\mu-$")
plt.savefig("img/z_boson_peak_noCuts.png", dpi = dpi)
plt.show()

print("Selecting Z Boson Resonance Peak: Ok!")

###############################
## Model: Breit-Wigner + Linear
###############################

pars = [2.3, 91.18, 6300, 1, 2, 50]

bwl      = Model(breitwigner) + Model(linear)
bwl_pars = bwl.make_params(gamma = pars[0], M = pars[1], A = pars[2], a = pars[4], b = pars[5])
bwl      = bwl.fit(y_zboson, params = bwl_pars, E = x_zboson, weights = np.sqrt(1.0/y_zboson), nan_policy = "omit")

M_bwl     = bwl.params['M'].value
gamma_bwl = bwl.params['gamma'].value
a_bwl     = bwl.params['a'].value
b_bwl     = bwl.params['b'].value

stderrM_bwl      = bwl.params['M'].stderr
stderr_gamma_bwl = bwl.params["gamma"].stderr
Nchi2_bwl        = bwl.redchi

entries = "Entries: {0:.0f}".format(len(z_boson.mass.values))
mean = "Mean: {0:.2f}".format(np.mean(z_boson.mass.values))
std = "Std Dev: {0:.2f}".format(np.std(z_boson.mass.values))
mass = "M = ({0:.2f} \u00B1 {1:.2f}) GeV/c²".format(M_bwl, stderrM_bwl)
decaywidth = "\u0393 = ({0:.2f} \u00B1 {1:.2f}) GeV".format(gamma_bwl, stderr_gamma_bwl)
chi2 = "χ²/ndf = {0:.2f}".format(Nchi2_bwl)

plt.text(62, 1600, entries, **style)
plt.text(62, 1530, mean, **style)
plt.text(62, 1460, std, **style)
plt.text(62, 1360, mass, **style)
plt.text(62, 1290, decaywidth, **style)
plt.text(62, 1220, chi2, **style)

plt.scatter(x_zboson, y_zboson, marker = ".", color = "black", label = "Data")
plt.hist(z_boson.mass.values, bins = nbins, color = "skyblue")
plt.plot(x_zboson, bwl.best_fit, color = "red", label = "PDF")
plt.plot(x_zboson, a_bwl*x_zboson + b_bwl, color = "green", label = "BKG")
plt.ylim(0)
plt.xlim(x_zboson.min(), x_zboson.max())
plt.legend()
plt.xlabel("Massa [GeV/c²]")
plt.ylabel("Eventos")
#plt.suptitle(r"Ajuste do pico de massa do bóson $Z \rightarrow \mu+\mu-$")
plt.title("Fundo: Função Linear\nSinal: Distribuição Breit-Wigner")
plt.savefig("img/model1:bretiwigner_linear.png", dpi = dpi)
plt.show()

print("Regression Model 1: Ok!")

#########################################################################
## Model: Pseudo Convoluted Relativistic Breit-Wigner - Gaussian + Linear
#########################################################################

pars = [3.5, 91.18, 6300, -1, 100, 0.7, 1]

psConvBgLin      = Model(pseudo_breitgauss) + Model(linear) 
psConvBgLin_pars = psConvBgLin.make_params(gamma = pars[0], M = pars[1], A = pars[2], a = pars[3], b = pars[4], alfa = pars[5], sigma = pars[6])
psConvBgLin      = psConvBgLin.fit(y_zboson, params = psConvBgLin_pars, E = x_zboson, weights = np.sqrt(1.0/y_zboson))

gamma_ps = psConvBgLin.params['gamma'].value
M_ps     = psConvBgLin.params['M'].value
a_ps     = psConvBgLin.params['a'].value
b_ps     = psConvBgLin.params['b'].value

stderrM_ps      = psConvBgLin.params['M'].stderr
stderr_gamma_ps = psConvBgLin.params["gamma"].stderr
Nchi2_ps        = psConvBgLin.redchi

entries = "Entries: {0:.0f}".format(len(z_boson.mass.values))
mean = "Mean: {0:.2f}".format(np.mean(z_boson.mass.values))
std = "Std Dev: {0:.2f}".format(np.std(z_boson.mass.values))
mass = "M = ({0:.2f} \u00B1 {1:.2f}) GeV/c²".format(M_ps, stderrM_ps)
decaywidth = "\u0393 = ({0:.2f} \u00B1 {1:.2f}) GeV".format(gamma_ps, stderr_gamma_ps)
chi2 = "χ²/ndf = {0:.2f}".format(Nchi2_ps)

plt.text(62, 1500, entries, **style)
plt.text(62, 1430, mean, **style)
plt.text(62, 1360, std, **style)
plt.text(62, 1260, mass, **style)
plt.text(62, 1190, decaywidth, **style)
plt.text(62, 1120, chi2, **style)

plt.scatter(x_zboson, y_zboson, marker = ".", color = "black", label = "Data")
plt.hist(z_boson.mass.values, bins = nbins, color = "lightgreen")
plt.plot(x_zboson, psConvBgLin.best_fit, color = "blue", label = "PDF")
plt.plot(x_zboson, a_ps*x_zboson + b_ps, color = "red", label = "BKG")
plt.ylim(0)
plt.xlim(x_zboson.min(), x_zboson.max())
plt.xlabel("Massa [GeV/c²]")
plt.ylabel("Eventos")
plt.legend()
#plt.suptitle(r"Ajuste do pico de massa do bóson $Z \rightarrow \mu+\mu-$")
plt.title("Fundo: Função Linear\nSinal: Pseudo convolução entre as distribuições Breit-Wigner e de Gauss")
plt.savefig("img/model2:pseudoConvoluted_breitwigner_gaussian_linear.png", dpi = dpi)
plt.show()

print("Regression Model 2: Ok!")

##################################################################
## Model: Convoluted Relativistic Breit-Wigner - Gaussian + Linear
##################################################################

pars = [2.3, 91.18, 6300, 1, 2, 50]

convBgLin      = Model(breitgauss) + Model(linear)
convBgLin_pars = convBgLin.make_params(gamma = pars[0], M = pars[1], A = pars[2], sigma = pars[3], a = pars[4], b = pars[5])
convBgLin      = convBgLin.fit(y_zboson, params = convBgLin_pars, E = x_zboson, weights = np.sqrt(1.0/y_zboson))

M_hw     = convBgLin.params['M'].value
gamma_hw = convBgLin.params['gamma'].value
A_hw     = convBgLin.params['A'].value
a_hw     = convBgLin.params['a'].value
b_hw     = convBgLin.params['b'].value
sigma_hw = convBgLin.params['sigma'].value

stderrM_hw      = convBgLin.params['M'].stderr
stderr_gamma_hw = convBgLin.params["gamma"].stderr
Nchi2_hw        = convBgLin.redchi

entries = "Entries: {0:.0f}".format(len(z_boson.mass.values))
mean = "Mean: {0:.2f}".format(np.mean(z_boson.mass.values))
std = "Std Dev: {0:.2f}".format(np.std(z_boson.mass.values))
mass = "M = ({0:.2f} \u00B1 {1:.2f}) GeV/c²".format(M_hw, stderrM_hw)
decaywidth = "\u0393 = ({0:.2f} \u00B1 {1:.2f}) GeV".format(gamma_hw, stderr_gamma_hw)
chi2 = "χ²/ndf = {0:.2f}".format(Nchi2_hw)

plt.text(62, 1500, entries, **style)
plt.text(62, 1430, mean, **style)
plt.text(62, 1360, std, **style)
plt.text(62, 1260, mass, **style)
plt.text(62, 1190, decaywidth, **style)
plt.text(62, 1120, chi2, **style)

plt.scatter(x_zboson, y_zboson, marker = ".", color = "black", label = "Data")
plt.hist(z_boson.mass.values, bins = nbins, color = "pink")
plt.plot(x_zboson, convBgLin.best_fit, color = "red", label = "PDF")
plt.plot(x_zboson, a_hw*x_zboson + b_hw, color = "green", label = "BKG")
plt.ylim(0)
plt.xlim(x_zboson.min(), x_zboson.max())
plt.legend()
plt.xlabel("Massa [GeV/c²]")
plt.ylabel("Eventos")
#plt.suptitle(r"Ajuste do pico de massa do bóson $Z \rightarrow \mu+\mu-$")
plt.title("Fundo: Função Linear\nSinal: Convolução entre as distribuições Breit-Wigner e de Gauss")
plt.savefig("img/model3:convoluted_breitwigner_gaussian_linear.png", dpi = dpi)
plt.show()

print("Regression Model 3: Ok!")

#############################################
## Selecting Z Boson Resonance Peak with Cuts
#############################################

### Selection criteria, pT > 25 and |eta| < 2.1
dimuon_pt_cut = momentumA[momentumA["pt"] > 25]
muon1_eta_cut = run2011A[np.abs(run2011A["eta1"]) < 2.1]
muon2_eta_cut = run2011A[np.abs(run2011A["eta2"]) < 2.1]

### Computing the spectrum and selecting the Z boson resonance peak
spectrum = pd.DataFrame( { 'dimu': inv_mass(dimuon_pt_cut['pt1'], dimuon_pt_cut['pt2'], muon1_eta_cut["eta1"], muon2_eta_cut["eta2"], run2011A['phi1'], run2011A['phi2']) } )
spectrum = spectrum[np.isnan(spectrum["dimu"]) == False]
z_boson  = pd.DataFrame( { 'mass': np.array([ mass for mass in spectrum["dimu"] if ((mass >= 60) and (mass <= 120)) ]) } )

### Plotting the resonance peak
nbins               = 200
y_zboson, bin_edges = np.histogram(z_boson.mass.values, bins = nbins, range = (60,120))
x_zboson            = (bin_edges[:-1] + bin_edges[1:])/2.

#plt.scatter(x_zboson, y_zboson, marker = '+', color = 'black')
plt.hist(z_boson.mass.values, bins = nbins, color = "red", label = "Entries: {0:.0f} \nMean: {1:.2f} \nStd Dev: {2:.2f}".format(len(z_boson.mass.values), np.mean(z_boson.mass.values), np.std(z_boson.mass.values)))
plt.xlim(60, 120)
plt.ylim(0)
plt.xlabel("Massa [GeV/c²]", fontsize = '12')
plt.ylabel("Eventos", fontsize = '12')
plt.legend(loc = "upper right", handlelength = 0, handletextpad = 0)
plt.title(r"Pico de massa do bóson $Z \rightarrow \mu+\mu-$")
plt.savefig("img/z_boson_peak_withCuts.png", dpi = dpi)
plt.show()

print("Selecting Z Boson Resonance Peak with Cuts: Ok!")

#############################################################################
## Model: Convoluted Relativistic Breit-Wigner and Crystal Ball + Exponential
#############################################################################

bbexp = Model(breitball) + Model(exponential)

pars = [3.5, 91.18, 106, 1, 0.9, 1.5, 5, 100]

pars_bbexp = bbexp.make_params(gamma = pars[0], M = pars[1], A = pars[2], sigma = pars[3], alpha = pars[4], n = pars[5], tau = pars[6], B = pars[7])
pars_bbexp["n"].set(pars[5], min = 1.01)
pars_bbexp["alpha"].set(pars[4], min = 0)
pars_bbexp["sigma"].set(pars[3], min = 0)
pars_bbexp["tau"].set(pars[6], min = 0)
res_bbexp = bbexp.fit(y_zboson, params = pars_bbexp, E = x_zboson, weights = np.sqrt(1.0/y_zboson), nan_policy = "omit")

gamma_bbexp = res_bbexp.params['gamma'].value
M_bbexp     = res_bbexp.params['M'].value
B_bbexp     = res_bbexp.params["B"].value
tau_bbexp   = res_bbexp.params['tau'].value

stderrM_bbexp = res_bbexp.params['M'].stderr
stderr_gamma_bbexp = res_bbexp.params["gamma"].stderr
Nchi2_bbexp  = res_bbexp.redchi

entries = "Entries: {0:.0f}".format(len(z_boson.mass.values))
mean = "Mean: {0:.2f}".format(np.mean(z_boson.mass.values))
std = "Std Dev: {0:.2f}".format(np.std(z_boson.mass.values))
mass = "M = ({0:.2f} \u00B1 {1:.2f}) GeV/c²".format(M_bbexp, stderrM_bbexp)
decaywidth = "\u0393 = ({0:.2f} \u00B1 {1:.2f}) GeV".format(gamma_bbexp, stderr_gamma_bbexp)
chi2 = "χ²/ndf = {0:.2f}".format(Nchi2_bbexp)
pt = r"$p_{T}$ > 25 GeV"
eta = r"$|\eta|$ < 2.1"

plt.text(62, 300, pt, **style)
plt.text(62, 285, eta, **style)
plt.text(62, 265, entries, **style)
plt.text(62, 250, mean, **style)
plt.text(62, 235, std, **style)
plt.text(62, 215, mass, **style)
plt.text(62, 200, decaywidth, **style)
plt.text(62, 185, chi2, **style)

plt.scatter(x_zboson, y_zboson, marker = ".", color = "black", label = "Data")
plt.hist(z_boson.mass.values, bins = nbins, color = "tomato")
plt.plot(x_zboson, res_bbexp.best_fit, color = "green", label = "PDF")
plt.plot(x_zboson, B_bbexp*np.exp(-x_zboson/tau_bbexp), color = "blue", label = "BKG")
plt.ylim(0)
plt.xlim(x_zboson.min(), x_zboson.max())
plt.xlabel("Massa [GeV/c²]")
plt.ylabel("Eventos")
plt.legend()
#plt.suptitle(r"Ajuste do pico de massa do bóson $Z \rightarrow \mu+\mu-$")
plt.title("Fundo: Função Exponencial\nSinal: Convolução entre as distribuições Breit-Wigner e Crystal Ball")
plt.savefig("img/model4:convoluted_breitwigner_crystalball_exponential.png", dpi = dpi)
plt.show()

print("Regression Model 4: Ok!")

####################################
## Mass measurement compability
####################################

Z_mass     = constant["z_mass"]
err_Z_mass = constant["error_z_mass"]

discrep_mass = np.abs(M_bbexp - Z_mass)
sigma_mass   = np.sqrt(np.power(stderrM_bbexp, 2) + np.power(err_Z_mass, 2))
C_mass       = discrep_mass/sigma_mass
E_mass       = stderrM_bbexp*100/np.abs(M_bbexp)

print("Compability calculation: Ok!")

####################################
## Number of events
####################################

pdf = res_bbexp.best_fit
exp = B_bbexp*np.exp(-x_zboson/tau_bbexp)

n_signal     = integrate.simps(pdf, x_zboson)
n_background = integrate.simps(exp, x_zboson)
n_Zmumu      = n_signal - n_background

err_n_signal     = np.sqrt(n_signal)
err_n_background = np.sqrt(n_background)
err_n_Zmumu      = np.sqrt(n_Zmumu)

print("Number of events calculation: Ok!")

###################################
## Reporting
###################################

print("\nReporting...")

print("\nModel 1:")
print("Z mass = {0:.2f} \u00B1 {1:.2f} \n\u0393 (Decay Width) = {2:.2f} \u00B1 {3:.2f} \nχ² = {4:.2f}".format(M_bwl, stderrM_bwl, gamma_bwl, stderr_gamma_bwl, Nchi2_bwl))

print("\nModel 2:")
print("Z mass = {0:.2f} \u00B1 {1:.2f} \n\u0393 (Decay Width) = {2:.2f} \u00B1 {3:.2f} \nχ² = {4:.2f}".format(M_ps, stderrM_ps, gamma_ps, stderr_gamma_ps, Nchi2_ps))

print("\nModel 3:")
print("Z mass = {0:.2f} \u00B1 {1:.2f} \n\u0393 (Decay Width) = {2:.2f} \u00B1 {3:.2f} \nχ² = {4:.2f}".format(M_hw, stderrM_hw, gamma_hw, stderr_gamma_hw, Nchi2_hw))

print("\nModel 4:")
print("Z mass = {0:.2f} \u00B1 {1:.2f} \n\u0393 (Decay Width) = {2:.2f} \u00B1 {3:.2f} \nχ² = {4:.2f}".format(M_bbexp, stderrM_bbexp, gamma_bbexp, stderr_gamma_bbexp, Nchi2_bbexp))

print("\nrelative error = {0:.2f}%".format(E_mass))
print("discrepancy = {0:.2f} sigma".format(C_mass))

print("\nn_signal = {0:.0f} \u00B1 {1:.0f}".format(n_signal, err_n_signal))
print("n_background = {0:.0f} \u00B1 {1:.0f}".format(n_background, err_n_background))
print("n_Zmumu = {0:.0f} \u00B1 {1:.0f}\n".format(n_Zmumu, err_n_Zmumu))