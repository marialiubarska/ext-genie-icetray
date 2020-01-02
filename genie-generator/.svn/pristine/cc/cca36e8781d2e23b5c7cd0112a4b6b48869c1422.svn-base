//____________________________________________________________________________
/*
 Author: Mark Dierckxsens - ULB

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cassert>
#include <algorithm>
#include <map>

#include <TH1D.h>
#include <TH1I.h>
#include <TF1.h>
#include <TFile.h>

#include "Conventions/Constants.h"
#include "GIceCubeDiffuseFlux.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::flux;

// change this to true for producing checking histograms
bool makeHistoCheck=false;
TH1I* hcubecount=0;
TH1I* hspherecount=0;

//____________________________________________________________________________
GIceCubeDiffuseFlux::GIceCubeDiffuseFlux()
{
  this->Initialize();
}
//___________________________________________________________________________
GIceCubeDiffuseFlux::~GIceCubeDiffuseFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
bool GIceCubeDiffuseFlux::GenerateNext(void)
{
  //-- Reset previously generated neutrino code / 4-p / 4-x
  this->ResetSelection();

  //-- Generate an energy from the 'combined' spectrum histogram
  //   and compute the momentum vector
  double Ev = (double) fTotSpectrum->GetRandom();

  TVector3 numom(1,1,1);  // momentum along the neutrino direction
  numom.SetTheta(TMath::Pi()-this->GenerateZenith());  // polar angle = pi - zenith
  numom.SetPhi(this->GenerateAzimuth());  // azimuth angle
  numom.SetMag(Ev);         // with |p|=Ev

  fgP4.SetPxPyPzE(numom.Px(), numom.Py(), numom.Pz(), Ev);

  //-- Select a neutrino species from the flux fractions at the
  //   selected energy
  fgPdgC = (*fPdgCList)[this->SelectNeutrino(Ev)];

  //-- Compute neutrino 4-x

  //double phi = this->GeneratePhi();  // rndm angle [0,2pi]
  //double R  = this->GenerateR();   // rndm R [0,fInjRadius]

  // This generates events inside a cylinder with fixed size and direction
  //double x,y=0.;
  //this->GenerateXY(x,y);
  //double z = this->GenerateZ();
  
  // Generate event inside a sphere with radius fInjRadius
  // then check if it is inside cylinder with direction parallel to the 
  // neutrino centered at the detector center and with radius fGenRadius 
  // and length fGenLength
  TVector3 nupos(0,0,0);
  bool incyl=false;
  int ncount=0;
  while (!incyl) {
    ncount++;
    nupos=this->GenerateXYZ();
    incyl=this->InsideGenCylinder(nupos,numom);
    if (ncount>100){
      LOG("IceCubeDiffuseFlux", pFATAL) << "Stuck in infinite loop ";
      break;
    }
  }
  if (hspherecount) hspherecount->Fill(ncount);

  fgX4.SetXYZT(nupos.X(),nupos.Y(),nupos.Z(),0.);
 
  LOG("IceCubeDiffuseFlux", pINFO) << "Generated neutrino pdg-code: " << fgPdgC;
  LOG("IceCubeDiffuseFlux", pINFO)
        << "Generated neutrino p4: " << utils::print::P4AsShortString(&fgP4);
  LOG("IceCubeDiffuseFlux", pINFO)
             << "Generated neutrino x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::Initialize(void)
{
  LOG("IceCubeDiffuseFlux", pNOTICE) << "Initializing GIceCubeDiffuseFlux driver";

  fMaxEv       = 0;
  fZenithMin   = 0; 
  fZenithMax   = TMath::Pi();
  fAzimuthMin  = 0; 
  fAzimuthMax  = 2*TMath::Pi();
  fDetDepth    = 1950;
  fGenRadius   = 0;
  fGenLength   = 0;
  fInjRadius   = 0;
  fPdgCList    = new PDGCodeList;
  fTotSpectrum = 0;

  this->ResetSelection();
  //this->SetRtDependence("x");
  //eg, other example: this->SetRtDependence("pow(x,2)");
  //
  if (makeHistoCheck) {
    hcubecount = new TH1I("hcubecount","Trials in cube",30,0,30);
    hspherecount = new TH1I("hspherecount","Trials in sphere",30,0,30);
  }
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::ResetSelection(void)
{
// initializing running neutrino pdg-code, 4-position, 4-momentum
  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);
}
//__________________________________________________________________________
double GIceCubeDiffuseFlux::InjectionRadius(double gl,double gr)
{
  // calculate the half length of the side of the cude of that just encompasses 
  // a cylinder with length gl and radius rl in all possible directions
  double rs = gl*gl/4 + gr*gr;
  return TMath::Sqrt(rs);
}

//___________________________________________________________________________
void GIceCubeDiffuseFlux::CleanUp(void)
{
  LOG("IceCubeDiffuseFlux", pNOTICE) << "Cleaning up...";

  if (fPdgCList   ) delete fPdgCList;
  if (fTotSpectrum) delete fTotSpectrum;

  unsigned int nspectra = fSpectrum.size();
  for(unsigned int i = 0; i < nspectra; i++) {
     TH1D * spectrum = fSpectrum[i];
     delete spectrum;
     spectrum = 0;
  }

  if (makeHistoCheck){
    // write out histograms
    TFile f("./histocounts.root","recreate");
    if (hcubecount) {
      hcubecount->Write();
      delete hcubecount;
    }
    if (hspherecount) {
      hspherecount->Write();
      delete hspherecount;
    }
    f.Close();
  }
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetZenithMin(double zm)
{
  if (zm<0 || zm>TMath::Pi()){
    LOG ("IceCubeDiffuseFlux", pWARN) << "Minimum Zenith angle outside range [0,pi]: " 
                        << zm << "\n  -> Will be set to 0 ";
    fZenithMin = 0;
    return;
  }
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting minimum zenith angle = " << zm;
  fZenithMin = zm;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetZenithMax(double zm)
{
  if (zm<0 || zm>TMath::Pi()){
    LOG ("IceCubeDiffuseFlux", pWARN) << "Maximum Zenith angle outside range [0,pi]: " 
                        << zm << "\n  -> Will be set to pi ";
    fZenithMax = TMath::Pi();
    return;
  }
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting maximum zenith angle = " << zm;
  fZenithMax = zm;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetAzimuthMin(double am)
{
  if (am<0 || am>2*TMath::Pi()){
    LOG ("IceCubeDiffuseFlux", pWARN) << "Minimum Azimuth angle outside range [0,2*pi]: " 
                        << am << "\n  -> Will be set to 0 ";
    fAzimuthMin = 0;
    return;
  }
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting minimum azimuthal angle = " << am;
  fAzimuthMin = am;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetAzimuthMax(double am)
{
  if (am<0 || am>2*TMath::Pi()){
    LOG ("IceCubeDiffuseFlux", pWARN) << "Maximum Azimuth angle outside range [0,2*pi]: " 
                        << am << "\n  -> Will be set to 2*pi ";
    fAzimuthMax = 2*TMath::Pi();
    return;
  }
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting maximum azimuthal angle = " << am;
  fAzimuthMax = am;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetGenerationRadius(double gr)
{
  if (gr<0){
    LOG ("IceCubeDiffuseFlux", pERROR) << "Only positive values for generation radius allowed";
    return;
  }
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting generation radius = " << gr;
  if (gr!=fGenRadius){       // generation radius changed, recalculate injection radius
    fInjRadius = this->InjectionRadius(fGenLength,gr);
    LOG("IceCubeDiffuseFlux", pDEBUG) << "Injection radius changed to " << fInjRadius;
  }
  fGenRadius = gr;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetGenerationLength(double gl)
{
  if (gl<0){
    LOG ("IceCubeDiffuseFlux", pERROR) << "Only positive values for generation length allowed ";
    return;
  }
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting generation length = " << gl;
  if (gl!=fGenLength){       // generation length changed, recalculate injection radius
    fInjRadius = this->InjectionRadius(gl,fGenRadius);
    LOG("IceCubeDiffuseFlux", pDEBUG) << "Injection radius changed to " << fInjRadius;
  }
  fGenLength = gl;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::SetDetectorDepth(double dd)
{
  LOG ("IceCubeDiffuseFlux", pNOTICE) << "Setting detector depth = " << dd;
  fDetDepth = dd;
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::AddEnergySpectrum(int nu_pdgc, TH1D * spectrum)
{
  LOG("IceCubeDiffuseFlux", pNOTICE) << "Adding flux spectrum for pdg = " << nu_pdgc;

  fPdgCList->push_back(nu_pdgc);

  bool accepted = (count(fPdgCList->begin(),fPdgCList->end(),nu_pdgc) == 1);
  if(!accepted) {
     LOG ("IceCubeDiffuseFlux", pWARN)
            << "The pdg-code isn't recognized and the spectrum was ignored";
  } else {
     fSpectrum.push_back(spectrum);

     int    nb  = spectrum->GetNbinsX();
     Axis_t max = spectrum->GetBinLowEdge(nb)+spectrum->GetBinWidth(nb);
     fMaxEv = TMath::Max(fMaxEv, (double)max);

     LOG("IceCubeDiffuseFlux", pNOTICE) 
          << "Updating maximum energy of flux particles to: " << fMaxEv;

     this->AddAllFluxes(); // update combined flux
  }
}
//___________________________________________________________________________
void GIceCubeDiffuseFlux::AddAllFluxes(void)
{
  LOG("IceCubeDiffuseFlux", pNOTICE) << "Computing combined flux";

  if(fTotSpectrum) delete fTotSpectrum;

  vector<TH1D *>::const_iterator spectrum_iter;

  unsigned int inu=0;
  for(spectrum_iter = fSpectrum.begin();
                       spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
     TH1D * spectrum = *spectrum_iter;

     if(inu==0) { fTotSpectrum = new TH1D(*spectrum); }
     else       { fTotSpectrum->Add(spectrum);        }
     inu++;
  }
}
//___________________________________________________________________________
int GIceCubeDiffuseFlux::SelectNeutrino(double Ev)
{
  const unsigned int n = fPdgCList->size();
  double fraction[n];

  vector<TH1D *>::const_iterator spectrum_iter;

  unsigned int inu=0;
  for(spectrum_iter = fSpectrum.begin();
                       spectrum_iter != fSpectrum.end(); ++spectrum_iter) {
     TH1D * spectrum = *spectrum_iter;
     fraction[inu++] = spectrum->GetBinContent(spectrum->FindBin(Ev));
  }

  double sum = 0;
  for(inu = 0; inu < n; inu++) {
     sum += fraction[inu];
     fraction[inu] = sum;
     LOG("IceCubeDiffuseFlux", pDEBUG) << "SUM-FRACTION(0->" << inu <<") = " << sum;
  }

  RandomGen * rnd = RandomGen::Instance();
  double R = sum * rnd->RndFlux().Rndm();

  LOG("IceCubeDiffuseFlux", pDEBUG) << "R e [0,SUM] = " << R;

  for(inu = 0; inu < n; inu++) {if ( R < fraction[inu] ) return inu;}

  LOG("IceCubeDiffuseFlux", pERROR) << "Could not select a neutrino species";
  assert(false);

  return -1;
}
//___________________________________________________________________________
double GIceCubeDiffuseFlux::GeneratePhi(void) const
{
  RandomGen * rnd = RandomGen::Instance();
  double phi = 2.*kPi * rnd->RndFlux().Rndm(); // [0,2pi]
  return phi;
}
//___________________________________________________________________________
double GIceCubeDiffuseFlux::GenerateR(void) const
{
  RandomGen * rnd = RandomGen::Instance();
  double r = fGenRadius * rnd->RndFlux().Rndm(); // [0,fGenRadius]
  return r;
}
//___________________________________________________________________________
/*
 * Generate at a random position in a circle with radius fGenRadius
 */
void GIceCubeDiffuseFlux::GenerateXY(double& x, double& y) const
{
  x=0;
  y=0;
  int ntrial=0;
  while(1){
    if (ntrial>=25){
      LOG("IceCubeDiffuseFlux",pWARN) << "No suitable X Y coordinate found after 25 trials, will return (0,0)";
      x=0; y=0;
      return;
    }
    RandomGen * rnd = RandomGen::Instance();
    x = 2*rnd->RndFlux().Rndm()-1;
    y = 2*rnd->RndFlux().Rndm()-1;
    ntrial++;
    if ( TMath::Sqrt(x*x+y*y) <= 1){
      x *= fGenRadius;
      y *= fGenRadius;
      return ;
    }
  }
}
//___________________________________________________________________________
/*
 * Generate at a random position along the generation length
 */
double GIceCubeDiffuseFlux::GenerateZ(void) const
{
  RandomGen * rnd = RandomGen::Instance();
  double z = (fGenLength)*rnd->RndFlux().Rndm(); // [0,GenLength]
  return z;
}
//___________________________________________________________________________
/*
 * Generate random events inside a sphere
 */
TVector3 GIceCubeDiffuseFlux::GenerateXYZ() const
{
  TVector3 pos(0,0,0);
  int ntrial=0;
  while(1){
    if (ntrial>=25){
      LOG("IceCubeDiffuseFlux",pWARN) << "No suitable X Y Z coordinate found after 25 trials, will return (0,0,0)";
      if (hcubecount) hcubecount->Fill(ntrial);
      pos.SetXYZ(0,0,0);
      return pos;
    }
    RandomGen * rnd = RandomGen::Instance();
    pos.SetX(2*rnd->RndFlux().Rndm()-1);
    pos.SetY(2*rnd->RndFlux().Rndm()-1);
    pos.SetZ(2*rnd->RndFlux().Rndm()-1);

    LOG("IceCubeDiffuseFlux",pDEBUG) << "Generated position (x,y,z,mag)" 
                                     << pos.X() << ", " << pos.Y() << ", " << pos.Z() 
                                     << ", " << pos.Mag2() << ", " << pos.Mag()  ;
    ntrial++;
    if ( pos.Mag2() <= 1){
      LOG("IceCubeDiffuseFlux",pDEBUG) << "Position accepted, setting mag to " << fInjRadius*pos.Mag();
      if (hcubecount) hcubecount->Fill(ntrial);
      pos.SetMag(fInjRadius*pos.Mag());
      return pos;
    }
  }  
}
//___________________________________________________________________________
bool GIceCubeDiffuseFlux::InsideGenCylinder(const TVector3& pos,const TVector3& mom) const
{
  // First determine the extreme points of the central axis of the cylinder which
  // is parallel to the neutrino direction. 
  TVector3 p1 = mom;
  p1.SetMag(fGenLength/2);

  // LOG("IceCubeDiffuseFlux",pDEBUG) << "Momentum of incoming neutrino " 
  //                                 << mom.PX() << ", " << mom.PY() << ", " << mom.PZ() << ", " << mom.P();
  // LOG("IceCubeDiffuseFlux",pDEBUG) << "Extreme point of cylinder axis " 
  //                                 << p1.X() << ", " << p1.Y() << ", " << p1.Z() << ", " << p1.Mag();
  // LOG("IceCubeDiffuseFlux",pDEBUG) << "Generation Sphere Radius: " << fInjRadius;
  
  //  
  // The vector v describing the line determined by this point and the center (thus the axis), 
  // is given by:
  //     v = p1 * u
  // with the endpoints of the cylinder (on the end planes) at u=-1 and 1
  // The distance between the test point pos and this line is:
  //     d**2 = (p1 * u - pos)**2
  // The closest distance, obtained by d(d**2)/du=0, determines u:
  //     u_min = (p1.pos)/L**2
  // Where L is the half length of the cylinder = p1**2
  // if u is inside [-1,1], then the test point pos is inside the cylinder along the 
  // axis direction (i.e. within the end planes). The only thing to check then is 
  // if the minimum distance is not larger than the radius of the cylinder. This
  // minimum distance is:
  //     d_min = sqrt(pos**2 - (pos.p1)**2/L**2)
  //           = sqrt(pos**2 - u**2 * L**2)
  //

  double half_l = fGenLength/2; // half length of cylinder
  double u_min = pos * p1;
  u_min /= (half_l*half_l);
  if (TMath::Abs(u_min)>1) 
    return false;
  double d_min = pos.Mag2() - half_l*half_l*u_min*u_min;
  if (d_min<0) {
    LOG ("IceCubeDiffuseFlux", pWARN) << "Minimum distance squared smaller than zero";
    d_min=0;
  }
  d_min = TMath::Sqrt(d_min);
  if (d_min>fGenRadius)
    return false;
  return true;
}
//___________________________________________________________________________
double GIceCubeDiffuseFlux::GenerateAzimuth(void) const
{
  RandomGen * rnd = RandomGen::Instance();
  double azim = (fAzimuthMax-fAzimuthMin) * rnd->RndFlux().Rndm(); 
  azim += fAzimuthMin;
  return azim; // [fAzimuthMin,fAzimuthMax]
}
//___________________________________________________________________________
double GIceCubeDiffuseFlux::GenerateZenith(void) const
{
  RandomGen * rnd = RandomGen::Instance();
  double cosmin = TMath::Cos(fZenithMin);
  double cosmax = TMath::Cos(fZenithMax);
    double zen = TMath::ACos(cosmin+(cosmax-cosmin)*rnd->RndFlux().Rndm()); 
  //double zen = TMath::Pi();
  return zen; // [fZentihMin,fZenithMax]
}

