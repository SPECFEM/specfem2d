#!/usr/bin/env python

# written by Paul Cristini for CNRS, LMA, Marseille, France in January 2012

from pylab import *
from numpy import *
import numpy.fft as f

######################################################
# USER PARAMETERS
######################################################
# source dominant frequency
f0_dominant_frequency = 20.0

# start time
t0_start_time = 1.0     # 2.8/f0_dominant_frequency

# seismogram length
T_length = 4.0

# receiver position
receiver_X = 3000.

######################################################

# Ricker source time function
def RickerF(t,t0,f0):
    a = pi*pi*f0*f0
    x = -2.0 * a * (1.0 - 2.0 * a *(t-t0)*(t-t0)) * exp(-a * (t-t0)*(t-t0))
    return x

def compute_solution():
    # frequency & start time
    f0 = f0_dominant_frequency
    t0 = t0_start_time

    # Size of the fourier transform
    N = 32*32*8
    T = T_length  # Length of the time signal
    t = linspace(0.,T,N)
    dt = t[1]-t[0]

    # source time function
    Rick = RickerF(t,t0,f0)
    # in frequency-domain
    Sf = f.fft(Rick)
    freq = f.fftfreq(N,d=dt)

    # plots figure
    fig = figure()
    ax = plot(freq[0:N/2],Sf[0:N/2],freq[0:N/2],abs(Sf[0:N/2]))

    #
    c0 = 1500
    k0 = 2*pi*f0/c0
    a0 = 0.000002

    print 'Q =',1/(c0*a0)

    # Position of the receiver
    PosX = receiver_X

    # Modification of the spectrum to include linear attenuation
    attf = k0*freq/f0 - 4 * a0 * freq * log((freq+0.0000001)/f0)
    Sf_att = Sf * exp(-1j * attf * PosX) * exp(-a0 * 2 * pi * freq * PosX)

    # Spectrum of the non attenuated signal
    Sf_0 = Sf * exp(-1j * 2 * pi * freq/c0 * PosX)

    #
    Sf_att_nan = nan_to_num(Sf_att)

    # plots figure
    plot(freq[0:N/2],Sf_att_nan[0:N/2],freq[0:N/2],abs(Sf_att_nan[0:N/2]))
    xlim(xmax=100)
    title('Spectra with and without attenuation')

    # signal in time-domain
    Sig_attn = f.ifft(Sf_att_nan)

    # plots figure
    figure(11)
    plot(t,Sig_attn)
    plot(t,f.ifft(Sf_0))
    grid()
    title('Snapshots with and without attenuation')
    show()


if __name__=='__main__':
    compute_solution()
