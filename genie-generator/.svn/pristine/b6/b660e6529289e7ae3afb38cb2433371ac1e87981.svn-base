from math import exp, pow, log, pi, cos

def RecalculateOneWeight(wdict):
    """Recalculate OneWeight from the values stored in the in the frame
    """
    try:
        TotalCrosssection=  wdict["TotalCrosssection"]
        InteractionColumnDepth = wdict["InteractionColumnDepth"]
        TotalColumnDepth = wdict["TotalColumnDepth"]
        tgtpdg = wdict["TargetPDG"]
        PrimaryNeutrinoEnergy = wdict["PrimaryNeutrinoEnergy"]
        PowerLawIndex = wdict["PowerLawIndex"]
        GenVolRadius = wdict["InjectionSurfaceR"]
        Emin = pow(10,wdict["MinEnergyLog"])
        Emax = pow(10,wdict["MaxEnergyLog"])
        ZenithMin  = wdict["MinZenith"]  
        ZenithMax  = wdict["MaxZenith"]  
        AzimuthMin = wdict["MinAzimuth"] 
        AzimuthMax = wdict["MaxAzimuth"] 
    except KeyError:
        print "Not all necessary variables found in weight dictionary"
        return {}

    # The values are taken from GENIE through PDG library
    if tgtpdg == 1000080160:
        tgtmass=2.65444e-23
    elif tgtpdg == 1000010010:
        tgtmass=1.67267e-24
    else: 
        print "Mass for target ", tgtpdg," not defined, please update oneweightcalc.py"
        return {}
    
    exponential_factor = (TotalCrosssection * 1.0e-31) * (InteractionColumnDepth / tgtmass);
    probability        = (TotalCrosssection * 1.0e-31) * (TotalColumnDepth / tgtmass) * exp(-exponential_factor);

    energyFactor = pow( PrimaryNeutrinoEnergy , -PowerLawIndex );

    if PowerLawIndex == 1:
        # if E^-1 then integral over Emin and Emax is
        energyIntegral = log(Emax/Emin)
    else:
        # if not E^-1 then integral over Emin and Emax is
        energyIntegral = (pow(Emax, (1.-PowerLawIndex)) -
                          pow(Emin, (1.-PowerLawIndex))) / (1.-PowerLawIndex)

    areaNorm = GenVolRadius*GenVolRadius*pi*1e4
    solidAngle = (cos(ZenithMin)-cos(ZenithMax)) *(AzimuthMax-AzimuthMin)
        
    OneWeight = (probability / energyFactor) * (energyIntegral * areaNorm * solidAngle);

    return {"OneWeight":OneWeight, "TotalInteractionProbabilityWeight": probability}

def UpdateOneWeightFrame(frame):
    """Update the OneWeight value in the frame
    """
    weightdict=frame["I3MCWeightDict"]
    UpdateOneWeightDict(weightdict)
    del frame["I3MCWeightDict"]
    frame["I3MCWeightDict"]=weightdict

def UpdateOneWeightDict(wdict):
    """Update the OneWeight value in the weight dict
    """
    recalc=RecalculateOneWeight(wdict)
    wdict["TotalInteractionProbabilityWeight"] = recalc["TotalInteractionProbabilityWeight"]
    wdict["OneWeight"] = recalc["OneWeight"]
