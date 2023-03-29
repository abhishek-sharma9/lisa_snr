# lisa_snr
This repository contains files that are used to calculate signal-to-noise ratio of supermassive BHB system in LISA.

- `Fplus.jpg`, `Fcross.jpg` and `Fmag.jpg` are plots for the LISA antenna pattern as a function of all sky positions (in SSB frame) at different instants of time assuming circular orbits of LISA around the Sun. I have used equations 47, 53 (with e = 0), 54 and 55 of [Cornish and Rubbo](https://arxiv.org/pdf/gr-qc/0209011.pdf) with the following change of coordinates from `$\theta, \; \phi, \;$` and `$\psi$` to `$\beta, \; \lambda$` and `$\psi$` (as mentioned on page 16 of [Synthetic LISA](https://arxiv.org/pdf/gr-qc/0407102.pdf)). 

`$$\beta = \pi/2 - \theta, \; \lambda = \phi, \; \psi = -\psi$$`

- `snr.jpg` shows the variation of SNR for all sky positions at different times. SNR is calculated using,
`$$\mathrm{SNR} = 4 \Re  \int_{f_l}^{f_u} \frac{\Tilde{h}^*(f)h(f)}{S_n(f)} \,df$$`

- I have used frequency domain approximant `IMRPhenomD`. The values of various parameters used to generate waveform are as follows:

`$$m_1 = m_2 = 10^6 M_\odot, \; d_L = 6500 \; Mpc, \; df = 10^{-6}, \; f_l = 10^{-4} Hz, \; f_u = 0.1 Hz$$`