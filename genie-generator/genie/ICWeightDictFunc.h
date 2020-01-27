/**
 *
 *
 *
 */

#ifndef ICWEIGHTDICTFUNC_H_INCLUDED
#define ICWEIGHTDICTFUNC_H_INCLUDED

// #include <sstream>
// #include <vector>
// #include <map>
// #include <string>

//ROOT
// #include <TSystem.h>
// #include <TVector3.h>

//GENIE
// #include "Conventions/Units.h"       // added
#include "EVGCore/EventRecord.h"
// #include "GHEP/GHepParticle.h"       //added
#include "EVGDrivers/GMCJDriver.h"

#include "genie/ICResultDict.h"

namespace IC_Helper
{
  /**
   *
   *
   *
   */

  // weight dict contents  
  double  _InjectionSurfaceR;                  // gOptCylinderRadius [m]
  double  _PowerLawIndex;                      // gOptPowerLawIndex
  double  _MinEnergyLog;                       // log10( gOptNuEnergyMin/GeV )
  double  _MaxEnergyLog;                       // log10( gOptNuEnergyMax/GeV )
  double  _MinZenith;                          // gOptZenithMin [rad]
  double  _MaxZenith;                          // gOptZenithMax [rad]
  double  _MinAzimuth;                         // gOptAzimuthMin [rad]
  double  _MaxAzimuth;                         // gOptAzimuthMax [rad]

  int32_t  _TargetPDGCode;                     // evt.Particle(1)->Pdg()
  double   _GENIEWeight;                       // evt.Weight()
  double   _GlobalProbabilityScale;            // mcj_driver->GlobProbScale()
  double   _GeneratorVolume;                   // _InjectionSurfaceR*_InjectionSurfaceR * pi * _TotalDetectionLength [pi*R^2*L]
  double   _Crosssection;                      // evt.XSec() / (1e-38*units::cm2)
  double   _InteractionProbabilityWeight;      // _GENIEWeight * evt.Probability()
  double   _TotalInteractionProbabilityWeight; // _GENIEWeight * _GlobalProbabilityScale

  double   _TotalDetectionLength;              // gOptCylinderLength [m]
  double   _LengthInVolume;                    // _TotalDetectionLength/2. * (u_min + 1.); u_min = int_pos * cyl_end / (half_l*half_l); half_l = _TotalDetectionLength/2.
  double   _EnergyLost;                        // 0.0 --why??
  double   _InteractionType;                   // CC->1, NC->2, else->0; is_CC = evt.Summary()->ProcInfo().IsWeakCC(); is_NC = evt.Summary()->ProcInfo().IsWeakNC();
  
  double   _PrimaryNeutrinoEnergy;             // k1.Energy(); const TLorentzVector & k1 = *(neutrino->P4()); GHepParticle * neutrino = evt.Probe();
  double   _OneWeight;                         // (_TotalInteractionProbability / energyFactor) * (energyIntegral * areaNorm * solidAngle); see ConvertToMCTree.cxx
  double   _NEvents;                           // static_cast<double>( gOptNevents ) []
  
  void     ClearMCWeight(weightmap_t weightdict_);
  void     CalculateMCWeight(genie::EventRecord * evt, GMCJDriver * mcj_driver, weightmap_t weightdict_);
  void     FillMCWeight(weightmap_t weightdict_);
  // double   CalculateColumnDepth(double d, weightmap_t weightdict_);
  void     WriteMCWeightToFile(FILE* f,int ievt, EventRecord* evt, weightmap_t wdm);
};
  
#endif //ICWEIGHTDICTFUNC_H_INCLUDED
