import numpy as np
import pycbc
import lal
from pycbc.fft import fft
from pycbc.filter import matchedfilter

class LISAantennaPattern:
    
    '''
    psi - polarisation angle
    beta - latitude in SSB frame.
    lamda - longitude in SSB frame.
    t - time
    '''
    
    Tperiod = 1.  #year
    omega = 2*np.pi/Tperiod
    
    def __init__(self, psi, beta, lamda, t):
        
        self.psi = psi
        self.beta = beta
        self.lamda = lamda
        self.t = t
        
    def Dplus(self):
        
        const1 = (np.sqrt(3)/64)
        
        D0 = -36*np.cos(self.beta)**2*np.sin(2*LISAantennaPattern.omega*self.t)
        D1 = (3 - np.cos(2*self.beta))*np.cos(2*self.lamda)*np.sin(4*LISAantennaPattern.omega*self.t)
        D2 = np.sin(2*self.lamda)*(np.cos(4*LISAantennaPattern.omega*self.t) - 9)
        D3 = -4*np.sqrt(3)*np.sin(2*self.beta)*(np.sin(3*LISAantennaPattern.omega*self.t - self.lamda) - 3*np.sin(LISAantennaPattern.omega*self.t + self.lamda))
        
        D_plus = const1*(D0 + D1 + D2 + D3)
        
        return D_plus
    
    def Dcross(self):
        
        const2 = 1/16.
        
        D0 = np.sqrt(3)*np.sin(self.beta)*(9*np.cos(2*self.lamda) - np.cos(4*LISAantennaPattern.omega*self.t - 2*self.lamda))
        D1 = -6*np.cos(self.beta)*(np.cos(3*LISAantennaPattern.omega*self.t - self.lamda) + 3*np.cos(LISAantennaPattern.omega*self.t + self.lamda))

        D_cross = const2*(D0 + D1)

        return D_cross

    def Fplus(self):
        
        F_plus = 0.5*(np.cos(2*self.psi)*self.Dplus() + np.sin(2*self.psi)*self.Dcross())
    
        return F_plus
    
    def Fcross(self):
        
        F_cross = 0.5*(-np.sin(2*self.psi)*self.Dplus() + np.cos(2*self.psi)*self.Dcross())

        return F_cross
    
def calculate_snr(signal, psd, f_lower, f_upper):
    
    '''
    This function calculates the integral using summation.
    
    inputs
    --------
    vec: signal frequency series.
    psd: Noise power spectral density frequency series.
    f_lower: lower limit of integral.
    f_upper: upper limit of integral.

    output
    --------
    snr: float.
    '''
    
    # if signal is a time series, then convert it into frequency series.
    # signal = matchedfilter.make_frequency_series(signal)
    delta_f = psd.delta_f
    
    if isinstance(signal, pycbc.types.TimeSeries):
        signal = signal.to_frequencyseries(delta_f = delta_f)
    
    # resize signal and psd to same length.
    if len(signal) < len(psd):
        
        signal.resize(psd.sample_rate/psd.delta_f // 2 + 1)
        
    elif len(psd) < len(signal):
        
        psd.resize(signal.sample_rate/signal.delta_f // 2 + 1)
        
    N = (len(signal) - 1) * 2
    kmin, kmax = matchedfilter.get_cutoff_indices(f_lower, f_upper, delta_f, N)
    
    Sum = (np.conjugate(signal)[kmin : kmax] * signal[kmin : kmax] / psd[kmin : kmax]).sum()
    result = np.sqrt(4 * Sum.real * delta_f)
    
    return result

# def calculate_snrUsingTrapz(signal, psd, f_lower, f_upper):
    
#     '''
#     This function calculates the integral using numpy.trapz.
    
#     inputs
#     --------
#     vec: signal frequency series.
#     psd: Noise power spectral density frequency series.
#     f_lower: lower limit of integral.
#     f_upper: upper limit of integral.

#     output
#     --------
#     signal power: float.
#     '''
    
#     # if signal is a time series, then convert it into frequency series.
#     # signal = matchedfilter.make_frequency_series(signal)
#     delta_f = psd.delta_f
    
#     if isinstance(signal, pycbc.types.TimeSeries):
#         signal = signal.to_frequencyseries(delta_f = delta_f)
    
#     # resize signal and psd to same length.
#     if len(signal) < len(psd):
        
#         signal.resize(psd.sample_rate/psd.delta_f // 2 + 1)
        
#     elif len(psd) < len(signal):
        
#         psd.resize(signal.sample_rate/signal.delta_f // 2 + 1)
        
#     N = (len(signal) - 1) * 2
#     kmin, kmax = matchedfilter.get_cutoff_indices(f_lower, f_upper, delta_f, N)
    
#     # Sum = (np.conjugate(signal)[kmin : kmax] * signal[kmin : kmax] / psd[kmin : kmax]).sum()
#     integrand = np.conjugate(signal)[kmin : kmax] * signal[kmin : kmax] / psd[kmin : kmax]
    
#     integral = np.trapz(integrand, dx = delta_f, axis = 0)
#     result = np.sqrt(4 * integral.real)
    
#     return result