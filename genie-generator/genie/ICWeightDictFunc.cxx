/**
 *
 *
 *
 */

#include "ICWeightDictFunc.h"

#include <sstream>
#include <vector>
#include <map>
#include <string>

//ROOT
// #include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>

//GENIE
// #include "Conventions/Units.h"       // added
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"       //added
#include "EVGDrivers/GMCJDriver.h"
#include "genie/ICResultDict.h"

//____________________________________________________________________________
void ClearMCWeight(weightmap_t weightdict_){
  weightdict_.clear();
}
//____________________________________________________________________________
void FillMCWeight(weightmap_t weightdict_){

  ClearMCWeight();

  // add checks for if weights dict values were changed from initial dummy values?

  // fill weight dict with predefined variables
  weightdict_["InjectionSurfaceR"]                 = _InjectionSurfaceR;
  weightdict_["PowerLawIndex"]                     = _PowerLawIndex;
  weightdict_["MinEnergyLog"]                      = _MinEnergyLog;
  weightdict_["MaxEnergyLog"]                      = _MaxEnergyLog;
  weightdict_["MinZenith"]                         = _MinZenith;
  weightdict_["MaxZenith"]                         = _MaxZenith;
  weightdict_["MinAzimuth"]                        = _MinAzimuth;
  weightdict_["MaxAzimuth"]                        = _MaxAzimuth;

  weightdict_["TargetPDGCode"]                     = _TargetPDGCode;
  weightdict_["GENIEWeight"]                       = _GENIEWeight;
  weightdict_["GlobalProbabilityScale"]            = _GlobalProbabilityScale;
  weightdict_["GeneratorVolume"]                   = _GeneratorVolume;
  weightdict_["Crosssection"]                      = _Crosssection;
  weightdict_["InteractionProbabilityWeight"]      = _InteractionProbabilityWeight;
  weightdict_["TotalInteractionProbabilityWeight"] = _TotalInteractionProbabilityWeight;

  weightdict_["TotalDetectionLength"]              = _TotalDetectionLength;
  weightdict_["LengthInVolume"]                    = _LengthInVolume;
  weightdict_["EnergyLost"]                        = _EnergyLost;
  weightdict_["InteractionType"]                   = _InteractionType;

  weightdict_["PrimaryNeutrinoEnergy"]             = _PrimaryNeutrinoEnergy;
  weightdict_["OneWeight"]                         = _OneWeight;
  weightdict_["NEvents"]                           = _NEvents;
}
//____________________________________________________________________________
void CalculateMCWeight(EventRecord * evt, GMCJDriver * mcj_driver, weightmap_t weightdict_){

  // values taken from input
  _NEvents              = static_cast<double>( gOptNevents );
  _PowerLawIndex        = gOptPowerLawIndex;
  _MinEnergyLog         = TMath::Log10(gOptNuEnergyMin);
  _MaxEnergyLog         = TMath::Log10(gOptNuEnergyMax);
  _MinZenith            = gOptZenithMin;
  _MaxZenith            = gOptZenithMax;
  _MinAzimuth           = gOptAzimuthMin;
  _MaxAzimuth           = gOptAzimuthMax;
  
  // geometry
  _InjectionSurfaceR    = gOptCylinderRadius;
  _TotalDetectionLength = gOptCylinderLength;

  _GeneratorVolume      = gOptCylinderRadius*gOptCylinderRadius * gOptCylinderLength; // V = pi*R^2*L

  // LengthInVolume -- projection of [beginning of the cylinder, vertex position] line onto neutrino direction (== cylinder axis)
  const double half_l = gOptCylinderLength/2.;  // half length of cylinder
  
  GHepParticle * neutrino = evt->Probe();
  const TLorentzVector & k1 = *(neutrino->P4());  // v 4-p (k1)                                                                                                                                              
  TLorentzVector * vtx = evt->Vertex();            // vertex in detector coord system

  TVector3 cyl_end(k1.Px(), k1.Py(), k1.Pz());    // direction of neutrino
  cyl_end.SetMag(half_l);                         // multiply by half length cylinder gives end point
  TVector3 cyl_beg = -1.*cyl_end;                 // the opposite point is the beginning
  TVector3 int_pos(vtx->X(), vtx->Y(), vtx->Z()); // interaction position

  // calculate the relative position along the cylinder axis of the perpendicular projection of the                                                                                                    
  // vertex onto this axis                                                                                                                                                                             
  const double u_min = int_pos * cyl_end / (half_l*half_l);
  if (fabs(u_min)>1.){
    LOG("gevgen", pFATAL) << "Vertex point outside generation cylinder - should not happen! u_min=" << u_min;
    exit(3);
  }

  _LengthInVolume = half_l * (u_min+1.); // [m] 

  // primary nu energy
  _PrimaryNeutrinoEnergy = k1.Energy();
  
  // this value is hardcoded to be 0.0 in genie-icetray; see https://wiki.icecube.wisc.edu/index.php/Neutrino_Generator/I3MCWeightDict#EnergyLost; -- do we need to keep this?
  _EnergyLost = 0.0;

  // target PDG
  GHepParticle * target = evt->Particle(1);
  _TargetPDGCode        = target->Pdg();

  // interaction type: CC=1, NC=2, else=0;
  const bool is_CC = evt->Summary()->ProcInfo().IsWeakCC();
  const bool is_NC = evt->Summary()->ProcInfo().IsWeakNC();
  if ((is_CC) && (is_NC)){
    LOG("gevgen", pFATAL) << "Internal error, event is CC *and* NC!";
    exit(3);
  }
  _InteractionType = 0;
  if (is_CC) {
    _InteractionType=1;
  } else if (is_NC) {
    _InteractionType=2;
  }

  // Total crosssection
  _Crosssection = evt->XSec() / (1e-27*units::cm2); // [mb]

  // Weights
  _GENIEWeight = evt->Weight();
  _GlobalProbabilityScale = mcj_driver->GlobProbScale();
  _InteractionProbabilityWeight = _GENIEWeight * evt->Probability();
  _TotalInteractionProbabilityWeight = _GENIEWeight * _GlobalProbabilityScale;

  // OneWeight
  const double areaNorm = TMath::Pi() * gOptCylinderRadius*gOptCylinderRadius * 1e4; // [cm^2]
  const double solidAngle = (TMath::Cos(_MinZenith)-TMath::Cos(_MaxZenith))*(_MaxAzimuth-_MinAzimuth);
  
  double energyIntegral=0;
  if (_PowerLawIndex == 1.) {
    // if E^-1 then integral over Emin and Emax is                                                                                                                                                   
    energyIntegral = TMath::Log(gOptNuEnergyMax/gOptNuEnergyMin);
  } else {
    // if not E^-1 then integral over Emin and Emax is                                                                                                                                               
    energyIntegral = (TMath::Power(gOptNuEnergyMax, (1.-_PowerLawIndex)) -
		      TMath::Power(gOptNuEnergyMin, (1.-_PowerLawIndex))) / (1.-_PowerLawIndex);
  }
    
  const double energyFactor = TMath::Power( _PrimaryNeutrinoEnergy, -_PowerLawIndex );

  _OneWeight = (_TotalInteractionProbabilityWeight / energyFactor) * (energyIntegral * areaNorm * solidAngle);

}
//____________________________________________________________________________
// double CalculateColumnDepth(double length, weightmap_t weightdict_)
// {
//   double coldepth = ICE_DENSITY * length; // is in g/cm^3 * m, multiply with 1e6 to get g/m^2                        
//   return coldepth*1e6;
// }
//____________________________________________________________________________
void WriteMCWeightToFile(FILE* ofile,int ievt, EventRecord* evt, weightmap_t  wdm)
{
  // first line in the event is the event number, the number of elements in the weightdict,                          
  // the # entries in the event tree, the pdg of incoming neutrino and the energy                                    
  // the last three are used for consistency checks with the hepevt file                                             
  unsigned int len_wdm = wdm.size();

  fprintf(ofile,"%d %u %d %d %e \n",ievt,len_wdm,evt->GetEntriesFast(),evt->Probe()->Pdg(),evt->Probe()->E());

  // now print all the members of the weightdict map                                                                 
  weightmap_t::iterator it;
  for (it=wdm.begin(); it != wdm.end(); it++){
    fprintf(ofile,"    %s %.9e \n",(*it).first.c_str(),(*it).second);
  }
}
//____________________________________________________________________________
