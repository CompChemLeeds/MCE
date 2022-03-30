# Hamiltonian Parameters File (Parameters are in atomic units)

# Spin Boson Parameters (exp cutoff, DL, UBO, or DL base for LHC)

SB={
    #Delta
    'SBDelta':'1.0d0',
    #Epsilon
    'SBEps':'0.0d0',
    #Beta
    'SBBeta':'5.0d0',

    #Exponential cutoff    
    'SBEwc':'2.5d0',
    #Kondo Paramter
    'SBEkondo':'0.09d0',
    #wmax exponential cutoff
    'SBEwmaxfact':'5.0d0',

    #wc drude lorentz cutoff
    'SBDwc':'0.05d0',
    #lambda value drude lorentz
    'SBDlambda':'2.5d0',
    #wmax for drude lorentz cutoff
    'SBDwmaxfact':'20.0d0',

    #Upper limit of norm
    'SBupnorm':'0.9999d0',
    #lower limit of norm
    'SBdownnorm':'0.9998d0',

    # 's' exponent value
    'SB_s_factor':'1.0d0',

    # w0 for underdamped brownian oscillator
    'SBBw0':'6.0d0',
    #Lambda underdamped brownian oscillator
    'SBBlambda':'0.20d0',
    
    'SBBGamma':'0.66666666666d0',
    
    'SBBwmaxfact':'2.0d0'
    
}
# Spin Boson LHC Parameters
SBLHC={
    'SBlambda_1':0.0125,
    'SBlambda_2':0.0083333333,
    'SBlambda_3':0.00625,
    'SBlambda_4':0.005,
    'SBlambda_5':0.004166666667,
    'SBlambda_6':0.00357143,
    'SBlambda_7':0.003125,
    'SBGamma_1':0.25,
    'SBGamma_2':0.25,
    'SBGamma_3':0.25,
    'SBGamma_4':0.25,
    'SBGamma_5':0.25,
    'SBGamma_6':0.25,
    'SBGamma_7':0.25,
    'SBw0_1':10.0,
    'SBw0_2':15.0,
    'SBw0_3':20.0,
    'SBw0_4':25.0,
    'SBw0_5':30.0,
    'SBw0_6':35.0,
    'SBw0_7':40.0
}
# Harmonic Potential Parameters
HP={
    'HPw':'1.0d0',
    'HPupnorm':'0.999d0',
    'HPdownnorm':'0.998d0'
}
# Free Particle Parameters
FP={
    'FPmass':'1.0d0',
    'FPupnorm': '0.999d0',
    'FPdownnorm': '0.998d0'
}
# Morse Potential Parameters
MP={
    'MPw':'1.0d0',
    'MPmass':'1.0d0',
    'MPDissEn':'10.25d0',
    'MPWellParam':'0.2209d0',
    'MPupnorm':'0.999d0',
    'MPdownnorm':'0.998d0'
}
# Inverted Gaussian Parameters
IV={
    'IVm': 1.0,
    'IVw': 0.05,
    'IVlambda': 0.5,
    'IVIntensity': 0.1,
    'IVupnorm': 0.999999,
    'IVdownnorm': 0.999998
}
# Coulomb Potential Parameters
CP={
    'CPm':1.0,
    'CPfrequency': 0.05,
    'CPRc': 0.0,
    'CPIntensity': 0.1,
    'CPupnorm':0.999,
    'CPdownnorm': 0.998
}
# Henon-Heiles Parameters
HH={
    'HHCoupling': 0.111803,
    'HHupnorm': 0.999,
    'HHdownnorm': 0.998
}
#Energy limit parameters
EL={
    'ECheck':'No',
    'Ntries':10,
    'Ebfmax':'5000000.0d0',
    'Ebfmin':'-5000000.0d0'
}

