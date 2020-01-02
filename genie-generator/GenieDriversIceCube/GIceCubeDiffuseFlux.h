//____________________________________________________________________________
/*!

\class   genie::flux::GIceCubeDiffuseFlux

\brief   Generates a flux uniformely between given polar and
         azimuthal angles, and configured by providing a TH1D histogram.
         Events are generated in a cylinder with radius R and length L centered 
         at the detector depth D (normally the depth of the center of the 
         detector) with axis parallel to the the neutrino direction. In practice 
         events are generated in a cube just large enough
         to hold cylinders of fixed sizes in all possible directions. The length
         of the side of the cube is given by Lc = 2*sqrt(R^2+L^2/4).
         Multiple neutrino species can be generated (you will need to supply
         an energy spectrum for each).
         Based on GCylindTH1Flux.


\author  Mark Dierckxsens <Mark.Dierckxsens \at ulb.ac.be>
         ULB

\created March 4, 2010

*/
//____________________________________________________________________________

#ifndef _G_ICECUBE_DIFFUSE_FLUX_H_
#define _G_ICECUBE_DIFFUSE_FLUX_H_

#include <string>
#include <vector>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>

#include "EVGDrivers/GFluxI.h"

class TH1D;
class TF1;
class TVector3;

using std::string;
using std::vector;

namespace genie {
namespace flux  {

class GIceCubeDiffuseFlux: public GFluxI {

public :
  GIceCubeDiffuseFlux();
 ~GIceCubeDiffuseFlux();

  // methods specific to this flux object
  void SetZenithMin        (double zm);
  void SetZenithMax        (double zm);
  void SetZenithMinDegrees (double zm);
  void SetZenithMaxDegrees (double zm);
  void SetAzimuthMin       (double am);
  void SetAzimuthMax       (double am);
  void SetAzimuthMinDegrees(double am);
  void SetAzimuthMaxDegrees(double am);
  void SetGenerationRadius  (double ir);
  void SetGenerationLength  (double il);
  void SetDetectorDepth    (double dd);
  void AddEnergySpectrum   (int nu_pdgc, TH1D * spectrum);
  

  // methods implementing the GENIE GFluxI interface
  const PDGCodeList &    FluxParticles (void) { return *fPdgCList; }
  double                 MaxEnergy     (void) { return  fMaxEv;    }
  double                 ZenithMin     (void) { return  fZenithMin; }
  double                 ZenithMax     (void) { return  fZenithMax; }
  double                 AzimuthMin    (void) { return  fAzimuthMin; }
  double                 AzimuthMax    (void) { return  fAzimuthMax; }
  double                 GenerationRadius(void) { return  fGenRadius; }
  double                 GenerationLength(void) { return  fGenLength; }
  double                 DetectorDepth  (void) { return  fDetDepth; }
  bool                   GenerateNext  (void);
  int                    PdgCode       (void) { return  fgPdgC;    }
  double                 Weight        (void) { return  1.0;       }
  const TLorentzVector & Momentum      (void) { return  fgP4;      }
  const TLorentzVector & Position      (void) { return  fgX4;      }
  bool                   End           (void) { return  false;     }

private:

  // private methods
  void   Initialize        (void);
  void   CleanUp           (void);
  void   ResetSelection    (void);
  /// Calculate the radius of the injection sphere/ half length of cube
  double InjectionRadius(double gL,double gR); 
  bool   InsideGenCylinder(const TVector3& pos,const TVector3& mom) const;
  void   AddAllFluxes      (void);
  int    SelectNeutrino    (double Ev);
  double GeneratePhi       (void) const;
  double GenerateAzimuth   (void) const;
  double GenerateZenith    (void) const;
  double GenerateR         (void) const;
  double GenerateZ         (void) const;
  void   GenerateXY(double& x, double& y) const;
  TVector3   GenerateXYZ() const;
  // private data members
  double         fMaxEv;       ///< maximum energy
  double         fZenithMin;   ///< min zenith angle (>=0) in rad, default 0
  double         fZenithMax;   ///< max zenith angle (<= pi) in rad, default pi
  double         fAzimuthMin;  ///< min azimuthal angle (>=0) in rad, default 0 
  double         fAzimuthMax;  ///< max azimuthal angle (<=2*pi) in rad, default 2*pi
  double         fDetDepth;    ///< depth of detector - center of cylinder
  double         fGenRadius;   ///< radius of cylinder - generation volume
  double         fGenLength;   ///< length of cylinder - generation volume
  double         fInjRadius;   ///< radius of sphere - injection volume 
  PDGCodeList *  fPdgCList;    ///< list of neutrino pdg-codes
  int            fgPdgC;       ///< running generated nu pdg-code
  TLorentzVector fgP4;         ///< running generated nu 4-momentum
  TLorentzVector fgX4;         ///< running generated nu 4-position
  vector<TH1D *> fSpectrum;    ///< flux = f(Ev), 1/neutrino species
  TH1D *         fTotSpectrum; ///< combined flux = f(Ev)

};

//___________________________________________________________________________
inline void GIceCubeDiffuseFlux::SetZenithMinDegrees(double zm)
{
  this->SetZenithMin(zm*TMath::DegToRad());
}
//___________________________________________________________________________
inline void GIceCubeDiffuseFlux::SetZenithMaxDegrees(double zm)
{
  this->SetZenithMax(zm*TMath::DegToRad());
}
//___________________________________________________________________________
inline void GIceCubeDiffuseFlux::SetAzimuthMinDegrees(double am)
{
  this->SetAzimuthMin(am*TMath::DegToRad());
}
//___________________________________________________________________________
inline void GIceCubeDiffuseFlux::SetAzimuthMaxDegrees(double am)
{
  this->SetAzimuthMax(am*TMath::DegToRad());
}

} // flux namespace
} // genie namespace


#endif // _G_ICECUBE_DIFFUSE_FLUX_H_
