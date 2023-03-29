import numpy as np
import pycbc
import h5py
from pycbc.waveform import get_fd_waveform
import lal
from pycbc import psd
from pycbc.fft import fft
from pycbc.filter import matchedfilter
from pycbc.psd.read import from_txt
from utils import *
from tqdm import tqdm

print('Reading psd file...')
PSD = np.loadtxt('LISA_psd.txt')
print('Done!')

# Generate a GW signal from Super-massive BHB system.

m1, m2 = 1e6, 1e6
distance = 6500 #Mpc
duration = 1000000
delta_f = 1/duration
fLow, fHigh = 1e-4, 0.1

print('Generating waveform in source frame.')
hpTilde, hcTilde = get_fd_waveform(approximant='IMRPhenomD', mass1=m1, mass2=m2, delta_f=delta_f, f_lower=fLow, f_final = fHigh, distance = distance)

print('Done!')

LISA_psd = from_txt('LISA_psd.txt', length = len(hpTilde),\
                  delta_f = delta_f, low_freq_cutoff = fLow, is_asd_file=False)

print('Reading D_plus and D_cross (detector tensors).')
with h5py.File('detectorTensor.hdf', 'r') as f:
    D_plus = f['/psi=piOver6/D_plus'][()]
    D_cross = f['/psi=piOver6/D_cross'][()]
    
psi = np.pi/6

# calculating antenna pattern factors from the detector tensor.

F_plus = 0.5*(np.cos(2*psi)*D_plus + np.sin(2*psi)*D_cross)
F_cross = 0.5*(-np.sin(2*psi)*D_plus + np.cos(2*psi)*D_cross)

F_plus = F_plus[0:-1: 30]    # As we are going to calculate SNR for few instants of time that are separated by one month.

print('Calculating SNR.')

SNR = np.empty(F_plus.shape)

for t in tqdm(range(F_plus.shape[0])):
    for b in range(F_plus.shape[1]):
        for l in range(F_plus.shape[2]):
            
            # Projected waveform
            hTilde = F_plus[t, b, l] * hpTilde + F_cross[t, b, l] * hcTilde
    

            SNR[t, b, l] = calculate_snr(hTilde, psd = LISA_psd, f_lower = fLow, f_upper = fHigh)
        
with h5py.File('LISA_SNR.hdf', 'w') as f:
    f.create_dataset('d=6500Mpc_snr', data = SNR, compression = 'gzip')