
# by Philip Knaute, June 2012

import numpy as np
import scipy.interpolate as interp
import scipy.integrate as integr
import scipy.signal as sgnl
import matplotlib.pyplot as plt
import ownModules.proc.filter as filt 

def calAdjCCTTFromTrace(nt,dt,tStartIn,tEndIn,dataIn, synthIn):
    """ calculate the cross correlation traveltime adjoint sources for one seismogram 
    IN:
        nt          : number of timesteps in each seismogram 
        dt          : timestep of seismograms 
        tStartIn      : float starting time for trace
        tEndIn        : float end time for trace
    OUT:
        fBar        : array containing the adjoint seismogram for the trace
        t           : ndarray containing the time steps 
    """
    isCalculateWeights = False
    if isCalculateWeights: 
        dSeism = np.zeros(nt)
        weight = 0

    # -- time vector
    t = np.ogrid[0:(nt-1)*dt:nt*1j]
    # -- the norm 
    norm = 0
        
    # -- numpy arrays initialisation 
    velSynth = np.zeros(nt)  
    accSynth = np.zeros(nt)  
    timeWind = np.zeros(nt)  
    fBar = np.zeros(nt)  

    # -- calculate time time-window
    tStart = tStartIn
    tEnd = tEndIn
    # -- the starting and ending sample numbers
    iStart = int(np.floor(tStart/dt))
    iEnd = int(np.ceil(tEnd/dt))
    # -- sample length of the window
    iWind = iEnd - iStart
    #print iStart,iEnd,iWind
    timeWind[iStart:iEnd]=sgnl.hann(iWind)

    # -- calculate the adjoint
    synth = synthIn
    interpTrc = interp.InterpolatedUnivariateSpline(t,synth)
    velSynth = interpTrc(t,1)
    accSynth = interpTrc(t,2) 

    integrArgument = timeWind*synth*accSynth
    # -- calculating the norm
    norm = integr.simps(integrArgument,dx=dt,axis=-1,even='last')
    
    # -- divide every trace (row in matrices) by their norm (row in vector norm)
    fBar = timeWind*velSynth / norm
    if  isCalculateWeights: 
        # -- read in the data seismograms
        data = dataIn
        # -- calculate the difference between data and synthetics (amplitude) per trace
        dSeism = data - synth                     
        # -- calculate the weight per trace
        integrArgument = timeWind*velSynth*dSeism
        weight = integr.simps(integrArgument,dx=dt,axis=-1,even='last')
        print "weight", weight/norm
        # -- multiply weight with every adj trace
        fBar = fBar*weight   
        print weight
    return [fBar,t]
