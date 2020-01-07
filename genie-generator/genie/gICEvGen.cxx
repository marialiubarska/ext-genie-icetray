//____________________________________________________________________________
/*!

\program gicevgen

\brief   IceCube generator based on the 'generic' GENIE v+A event generation driver (gicevgen)

	 This event generation driver can handle event generation for neutrinos with random
         directions following a spectral flux between a given energy range on a single target
         (water by default).
         Events are generated in a cylinder with radius R and length L centered at the detector
         depth D (normally the depth of the center of the detector) with axis parallel to the 
         the neutrino direction. In practice events are generated in a sphere just large enough
         to hold cylinders of fixed sizes in all possible directions. The radius of this sphere 
         is given by Rs = sqrt(R^2+L^2/4).

         Syntax :
           gicevgen [-h] [-r runnum] -n nev -f neg_spectral_index -e emin,emax -p nupdg 
                    [-t tgtmix] [-D depth] [-R radius] [-L length] [-Z zmin,zmax] [-A amin,amax]

         Options :
           [] Denotes an optional argument.
           -h Prints-out help on using gicevgen and exits.
           -n Specifies the number of events to generate.
           -r Specifies the MC run number.
           -f flux: negative (!) of spectral index. a in x**(-1*a)
           -e Specifies the neutrino energy range: a comma separated pair of values
           -p Specifies the neutrino PDG code.
           -t Specifies the target PDG code (pdg format: 10LZZZAAAI) _or_ a target
              mix (pdg codes with corresponding weights) typed as a comma-separated 
              list of pdg codes with the corresponding weight fractions in brackets, 
              eg code1[fraction1],code2[fraction2],... 
              For example, to use a target mix of 95% O16 and 5% H type: 
              `-t 1000080160[0.95],1000010010[0.05]'.
              Default is water 1000080160[0.899],1000010010[0.111]
           -D depth of center generation volume in m, positive number! default is 1950
           -R radius of cylindrical generation volume in m, default is 1200
           -L length of cylindrical generation volume in m, default is 2000
           -Z range of zenith angles in degrees: comma separated pair of values
           -A range of azimuthal angles in degrees: comma separated pair of values
           

\author  Mark Dierckxsens - ULB

\created January 08, 2010

*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "Conventions/XmlParserStatus.h"
#include "Conventions/GBuild.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/AppInit.h" // new
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
// #include "Utils/CmdLineArgParserUtils.h" // outdated                                              
// #include "Utils/CmdLineArgParserException.h" // outdated                                          
#include "Utils/CmdLnArgParser.h" // new 

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#include "FluxDrivers/GCylindTH1Flux.h"
#include "FluxDrivers/GMonoEnergeticFlux.h"
//#include "Geo/PointGeomAnalyzer.h"
#include "GenieDriversIceCube/GIceCubeDiffuseFlux.h"
#include "GenieDriversIceCube/SimpleIceCubeGeomAnalyzer.h"
//#include "GIceCubeDiffuseFlux.h"
#endif
#endif

#include "HelperClasses/CrossSectionAccessor.h"

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;
using namespace std;

typedef map<string,double> weightmap_t;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

void            GenerateEvents();
GeomAnalyzerI * GeomDriver              (void);
GFluxI *        FluxDriver              (void);
GFluxI *        MonoEnergeticFluxDriver (void);
GFluxI *        TH1FluxDriver           (void);

void     ClearMCWeight(void);
void     FillMCWeight(EventRecord* evt);
double CalculateColumnDepth(double d);
void WriteMCWeightToFile(FILE* f,int ievt, EventRecord* evt, weightmap_t wdm);

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 0;       // n-events to generate
NtpMCFormat_t kDefOptNtpFormat = kNFGHEP; // ntuple format
Long_t        kDefOptRunNu     = 0;       // default run number

//User-specified options:
int             gOptNevents;      // n-events to generate
double          gOptNuEnergy;     // neutrino E, or min neutrino energy in spectrum
double          gOptNuEnergyRange;// energy range in input spectrum
int             gOptNuPdgCode;    // neutrino PDG code
map<int,double> gOptTgtMix;       // target mix (each with its relative weight)
Long_t          gOptRunNu;        // run number
double          gOptGenVolRadius=1200; // radius of injection cilinder 
double          gOptGenVolDepth=1950;  // Depth of the center of generation volume
double          gOptGenVolLength=2000;  // Length of generation volume
double          gOptPowerLawIndex=2;  // absolute value |x| of exponent E**(-x) 
double          gOptNuZenithMin=0; // minimum zenith angle
double          gOptNuZenithMax=180; // maximum zenith angle
double          gOptNuAzimuthMin=0; // minimum azimith angle
double          gOptNuAzimuthMax=360; // maximum azimuth angle

// weight dict map
weightmap_t weightdict_;

// densities in g/cm^3
double ICE_DENSITY=0.93;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);
  
  //-- Autoload splines (from the XML file pointed at the $GSPLOAD env. var.,
  //   if the env. var. has been set)
  // XSecSplineList * xspl = XSecSplineList::Instance();
  // xspl->AutoLoad(); //

  //-- Generate neutrino events
  //
  GenerateEvents();

  return 0;
}
//____________________________________________________________________________
void GenerateEvents(void)
{
  //-- get flux and geom drivers 
  GFluxI *        flux_driver = FluxDriver(); 
  GeomAnalyzerI * geom_driver = GeomDriver();

  //-- create the monte carlo job driver
  GMCJDriver * mcj_driver = new GMCJDriver;
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();

  //-- initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  
  string fn_prefix="genie_ic";
  ntpw.Initialize();
  ntpw.CustomizeFilenamePrefix(fn_prefix.c_str());
  
  //-- open output file to write weight_dict
  ostringstream wdic_filename;
  wdic_filename << fn_prefix  << "."
                << gOptRunNu << ".wdict.dat";

  FILE* wd_ofile;
  wd_ofile = fopen(wdic_filename.str().c_str(),"w");

  //initialize cross section accessor
  string xsec_file_root = "~/Software/genie_workspace/genie-generator/xsec_splines/GENIE_2_12_8_Water_splines.xml";
  utils::app_init::XSecTable(xsec_file_root, false);
  
  // if (gSystem->Getenv("GSPLOAD")){
  //     string xsec_file_xml = gSystem->Getenv("GSPLOAD");
  //     LOG("gicevgen", pDEBUG) << "XML cross section spline file: " << xsec_file_xml;
  //     size_t ppos = xsec_file_xml.find_last_of(".");
  //     xsec_file_root = xsec_file_xml.substr(0,ppos) + ".root";
  //     LOG("gicevgen", pDEBUG) << "ROOT cross section file name: " << xsec_file_root;
  // }

  CrossSectionAccessor::GetMe()->SetCrossSectionFile(xsec_file_root.c_str());

  //// Create an extra output file containing data necessary for i3weights
  //// modify the filename to add the run number & the ntuple format
  //ostringstream filename;
  //filename << filename_prefix  << "." 
  //         << gOptRunNu << ".MCWD.root";
  //TFile* outfile_w = new TFile(filename.str().c_str(),"RECREATE");
  //
  //TTree* wdtree = new TTree("wdtree","MC Weight Dict Tree");
  //
  //
  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gicevgen", pINFO) << " *** Generating event............ " << ievent;

     // generate a single event for neutrinos coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gicevgen", pINFO) << "Generated Event GHEP Record: " << *event;

     FillMCWeight(event);

     // add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);

     WriteMCWeightToFile(wd_ofile,ievent,event,weightdict_);

     ievent++;
     delete event;
  }

  //-- save the generated MC events
  ntpw.Save();

  //-- close weight dict output file
  fclose(wd_ofile);

  delete flux_driver;
  delete geom_driver;
  delete mcj_driver;

  //  delete outfile_w;
  // delete wdtree;
}
//____________________________________________________________________________
GeomAnalyzerI * GeomDriver(void)
{
// create a trivial point geometry with the specified target or target mix

  GeomAnalyzerI * geom_driver = new geometry::SimpleIceCubeGeomAnalyzer(gOptTgtMix);
  LOG("gicevgen", pDEBUG) << "Pointer to geo driver: " << geom_driver;

  return geom_driver;
}
//____________________________________________________________________________
GFluxI * FluxDriver(void)
{
// create & configure one of the generic flux drivers
//
  GFluxI * flux_driver = 0;

  flux_driver = TH1FluxDriver();

  return flux_driver;
}
//____________________________________________________________________________
GFluxI * TH1FluxDriver(void)
{
// 
//
  flux::GIceCubeDiffuseFlux * flux = new flux::GIceCubeDiffuseFlux;
  TH1D * spectrum = 0;

  int flux_entries = 10000000;

  double emin = gOptNuEnergy;
  double emax = gOptNuEnergy+gOptNuEnergyRange;
  //
  // ** generate the flux histogram from the input functional form
  //
  ostringstream oss;
  oss << gOptPowerLawIndex;
  string flux_formula = "x**(-1*" + oss.str() + ")";

  TF1 *  input_func = new TF1("input_func", flux_formula.c_str(), emin, emax);
  spectrum  = new TH1D("spectrum","neutrino flux", 300, emin, emax);
  spectrum->SetDirectory(0);
  spectrum->FillRandom("input_func", flux_entries);
  delete input_func;

  TFile f("./input-flux.root","recreate");
  spectrum->Write();
  f.Close();

  // set options 
  flux->SetZenithMinDegrees(gOptNuZenithMin);
  flux->SetZenithMaxDegrees(gOptNuZenithMax);
  flux->SetAzimuthMinDegrees(gOptNuAzimuthMin);
  flux->SetAzimuthMaxDegrees(gOptNuAzimuthMax);
  flux->SetGenerationRadius(gOptGenVolRadius);
  flux->SetGenerationLength(gOptGenVolLength);
  flux->SetDetectorDepth(gOptGenVolDepth);
  flux->AddEnergySpectrum(gOptNuPdgCode, spectrum);

  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  return flux_driver;
}
//............................................................................

//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gicevgen", pINFO) << "Parsing command line arguments";

  CmdLnArgParser gOptInp(argc, argv);

  // help?                                                                                            
  bool help = gOptInp.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // number of events:
  try {
    LOG("gicevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = gOptInp.ArgAsInt('n');
    // gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(...) {
    if(!gOptInp.OptionExists('n')) {
      LOG("gicevgen", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  // run number:
  try {
    LOG("gicevgen", pINFO) << "Reading MC run number";
    gOptRunNu = gOptInp.ArgAsInt('r');
    // gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(...) {
    if(!gOptInp.OptionExists('r')) {
      LOG("gicevgen", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  // flux functional form
  try {
    LOG("gicevgen", pINFO) << "Reading flux function";
    gOptPowerLawIndex = gOptInp.ArgAsDouble('f');
    // gOptPowerLawIndex = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'f');
  } catch(...) {
    if(!gOptInp.OptionExists('f')) {
      LOG("gicevgen", pFATAL) << "Unspecified (negative of) spectral index";
      PrintSyntax();
      exit(1);
    }
  }

  // neutrino energy:
  try {
    LOG("gicevgen", pINFO) << "Reading neutrino energy";
    string nue = gOptInp.ArgAsString('e');
    // string nue = genie::utils::clap::CmdLineArgAsString(argc,argv,'e');

    // is it just a range (comma separated set of values)
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);   
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0);
       gOptNuEnergy      = emin;
       gOptNuEnergyRange = emax-emin;
    } else {
      LOG("gicevgen", pFATAL) << "Invalid format energy range. Must be -e Emin,Emax - Exiting";
      PrintSyntax();
      exit(1);
    }
  } catch(...) {
    if(!gOptInp.OptionExists('e')) {
      LOG("gicevgen", pFATAL) << "Unspecified neutrino energy range - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // neutrino PDG code:
  try {
    LOG("gicevgen", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = gOptInp.ArgAsInt('p');
    // gOptNuPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(...) {
    if(!gOptInp.OptionExists('p')) {
      LOG("gicevgen", pFATAL) << "Unspecified neutrino PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // target mix (their PDG codes with their corresponding weights):
  try {
    LOG("gicevgen", pINFO) << "Reading target mix";
    string stgtmix = gOptInp.ArgAsString('t');
    // string stgtmix = genie::utils::clap::CmdLineArgAsString(argc,argv,'t');
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(stgtmix,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
    } else {
      vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
      for( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
         string tgt_with_wgt = *tgtmix_iter;
         string::size_type open_bracket  = tgt_with_wgt.find("[");
         string::size_type close_bracket = tgt_with_wgt.find("]");
         string::size_type ibeg = 0;
         string::size_type iend = open_bracket;
         string::size_type jbeg = open_bracket+1;
         string::size_type jend = close_bracket-1;
         int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend).c_str());
         double wgt = atof(tgt_with_wgt.substr(jbeg,jend).c_str());
         LOG("Main", pINFO)
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
      }//tgtmix_iter
    }//>1

  } catch(...) {
    if(!gOptInp.OptionExists('t')) {
      LOG("gicevgen", pINFO) << "Unspecified target mix - Using Default";
      gOptTgtMix.insert(map<int, double>::value_type(1000080160, 0.889));
      gOptTgtMix.insert(map<int, double>::value_type(1000010010, 0.111));
    }
  }


  // radius of generation volume:
  try {
    LOG("gicevgen", pINFO) << "Reading radius of cylindrical generation volume";
    gOptGenVolRadius = gOptInp.ArgAsDouble('R');
    // gOptGenVolRadius = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'R');
  } catch(...) {
      LOG("gicevgen", pINFO) << "Unspecified radius - Using Default";
  }

  // depth of generation volume:
  try {
    LOG("gicevgen", pINFO) << "Reading depth of deneration volume center";
    gOptGenVolDepth = gOptInp.ArgAsDouble('D');
    // gOptGenVolDepth = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'D');
  } catch(...) {
      LOG("gicevgen", pINFO) << "Unspecified depth - Using Default";
  }
  // height of generation volume
  try {
    LOG("gicevgen", pINFO) << "Reading length of cylindrical generation volume";
    gOptGenVolLength = gOptInp.ArgAsDouble('L');
    // gOptGenVolLength = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'L');
  } catch(...) {
      LOG("gicevgen", pINFO) << "Unspecified length - Using Default";
  }

  // Minimum and maximum zenith angle
  try {
    LOG("gicevgen", pINFO) << "Reading zenith angles range in degrees";
    string nuzen = gOptInp.ArgAsString('Z');
    // string nuzen = genie::utils::clap::CmdLineArgAsString(argc,argv,'Z');

    // is it just a value or a range (comma separated set of values)
    if(nuzen.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nuzen, ",");
       assert(nurange.size() == 2);   
       double zenmin = atof(nurange[0].c_str());
       double zenmax = atof(nurange[1].c_str());
       assert(zenmax>zenmin && zenmin>=0 && zenmax<=180);
       gOptNuZenithMin = zenmin;
       gOptNuZenithMax = zenmax;
    } else {
      LOG("gicevgen", pFATAL) << "Invalid format zenith range. Must be -Z ZenithMin,Zenithmax - Exiting";
      PrintSyntax();
      exit(1);
    }
  } catch(...) {
    // Assume full range
    LOG("gicevgen", pINFO) << "No option found for zenith angles, using default [0,180] ";
    gOptNuZenithMin=0;
    gOptNuZenithMax=180;
  }

  // Minimum and maximum azimuth angle
  try {
    LOG("gicevgen", pINFO) << "Reading azimuth angles range in degrees";
    string nuaz = gOptInp.ArgAsString('A');
    // string nuaz = genie::utils::clap::CmdLineArgAsString(argc,argv,'A');

    // is it just a value or a range (comma separated set of values)
    if(nuaz.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nuaz, ",");
       assert(nurange.size() == 2);   
       double azmin = atof(nurange[0].c_str());
       double azmax = atof(nurange[1].c_str());
       assert(azmax>azmin && azmin>=0 && azmax<=360);
       gOptNuAzimuthMin = azmin;
       gOptNuAzimuthMax = azmax;
    } else {
      LOG("gicevgen", pFATAL) << "Invalid format azimuth range. Must be -A AzimuthMin,AzimuthMax - Exiting";
      PrintSyntax();
      exit(1);
    }
  } catch(...) {
    // Assume full range
    LOG("gicevgen", pINFO) << "No option found for azimuth angles, using default [0,360] ";
    gOptNuAzimuthMin=0;
    gOptNuAzimuthMax=360;
  }

  // print-out the command line options
  //
  LOG("gicevgen", pNOTICE) << "MC Run Number              = " << gOptRunNu;
  LOG("gicevgen", pNOTICE) << "Number of events requested = " << gOptNevents;
  LOG("gicevgen", pNOTICE) << "Neutrino PDG code          = " << gOptNuPdgCode;
  LOG("gicevgen", pNOTICE) << "Flux =  E**(-" << gOptPowerLawIndex << ")";
  LOG("gicevgen", pNOTICE) << "Neutrino energy            = [" 
        << gOptNuEnergy << ", " << gOptNuEnergy+gOptNuEnergyRange << "]";
  LOG("gicevgen", pNOTICE) << "Radius generation volume   = " << gOptGenVolRadius;
  LOG("gicevgen", pNOTICE) << "depth generation volume    = " << gOptGenVolDepth;
  LOG("gicevgen", pNOTICE) << "length generation volume   = " << gOptGenVolLength;
  LOG("gicevgen", pNOTICE) << "Zenith angles              = [" 
                       << gOptNuZenithMin << ", " << gOptNuZenithMax << "]";
  LOG("gicevgen", pNOTICE) << "Azimuth angles              = [" 
                       << gOptNuAzimuthMin << ", " << gOptNuAzimuthMax << "]";

  map<int,double>::const_iterator iter;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
      int    tgtpdgc = iter->first;
      double wgt     = iter->second;
      LOG("gicevgen", pNOTICE) 
          << "Target mix - element = " <<  tgtpdgc << ", wgt = " << wgt;
  }

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gicevgen", pWARN)
    << "\n\n" << "Syntax:" << "\n"
    << "   gicevgen [-h] [-r runnum] -n nev -f neg_spectral_index [-s] -e emin,emax -p nupdg [-t tgtmix] [-R radius] [-D depth] [-L length] [-Z zmin,zmax] [-A amin,amax] \n";
}

//____________________________________________________________________________
void ClearMCWeight(){
  weightdict_.clear();
}

//____________________________________________________________________________
void FillMCWeight(EventRecord * evt){

  ClearMCWeight();

  weightdict_["InjectionSurfaceR"] = gOptGenVolRadius;
  weightdict_["PowerLawIndex"]     = gOptPowerLawIndex;

  weightdict_["MinEnergyLog"]      = TMath::Log10(gOptNuEnergy);
  weightdict_["MaxEnergyLog"]      = TMath::Log10(gOptNuEnergy+gOptNuEnergyRange);
  weightdict_["MinZenith"]         = gOptNuZenithMin*TMath::DegToRad();
  weightdict_["MaxZenith"]         = gOptNuZenithMax*TMath::DegToRad();       
  weightdict_["MinAzimuth"]        = gOptNuAzimuthMin*TMath::DegToRad();
  weightdict_["MaxAzimuth"]        = gOptNuAzimuthMax*TMath::DegToRad();

  // not filled
  //weightdict_["ActiveLengthBefore"]= 0;
  //weightdict_["ActiveLengthAfter"] = 0;
  //weightdict_["ActiveLengthAfter"] = 0;

  //LOG("gicevgen", pDEBUG) << "Pointer to target nulceus: " << evt->TargetNucleus();
  //LOG("gicevgen", pDEBUG) << "TargetNucleus pdg: " << evt->TargetNucleus()->Pdg();

  int nupdg=evt->Probe()->Pdg();
  int tgtpdg=1000010010; // proton by default, pointer to target nucleus not filled if it is a proton
  if (evt->TargetNucleus()) tgtpdg=evt->TargetNucleus()->Pdg();
  const double PrimaryNeutrinoEnergy=evt->Probe()->E();
  // accessor returns in units of 1e-38cm2, 
  //  converted in mb = 1e-31 m2 = 1e-27cm2= 1e11 * 1e-38 cm2  
  const double TotalCrosssection=CrossSectionAccessor::GetMe()->GetTotalCrossSection(nupdg,tgtpdg,PrimaryNeutrinoEnergy) * 1e-11; 
  
  //calculate the beginning and end position of the cylinder
  const double half_l = gOptGenVolLength/2; // half length of cylinder
  TVector3 cyl_end = evt->Probe()->P4()->Vect(); // direction of neutrino
  cyl_end.SetMag(half_l);                 // multiply by half length cylinder gives end point
  TVector3 cyl_beg = -1*cyl_end;                // the opposite point is the beginning
  TVector3 int_pos = evt->Vertex()->Vect(); // interaction position

  // calculate the relative position along the cylinder axis of the perpendicular projection of the
  // vertex onto this axis
  double u_min = int_pos * cyl_end;
  u_min /= (half_l*half_l);
  if (TMath::Abs(u_min)>1) 
    LOG("gicevgen", pERROR) << "Vertex point outside generation cylinder - should not happen! u_min = " << u_min;
  const double LengthInVolume = half_l * (u_min+1);
  LOG("gicevgen", pDEBUG) << "LengtInVolume, u_min  " << LengthInVolume << " " << u_min;

  const double InteractionColumnDepth = CalculateColumnDepth(LengthInVolume);
  const double TotalColumnDepth = CalculateColumnDepth(gOptGenVolLength);
  LOG("gicevgen", pDEBUG) << "Column depth, total interaction  " << InteractionColumnDepth 
                          << " " << TotalColumnDepth;

  double tgtmass = CrossSectionAccessor::GetMe()->GetMassInGram(tgtpdg);
  LOG("gicevgen", pDEBUG) << "Mass of target = " << tgtmass;

  const double exponential_factor = (TotalCrosssection * 1.0e-31) * (InteractionColumnDepth / tgtmass);
  const double probability        = (TotalCrosssection * 1.0e-31) * (TotalColumnDepth / tgtmass) * exp(-exponential_factor);

  LOG("gicevgen", pDEBUG) << "Total interaction probability weight = " << probability;

  double InteractionType = 0;
  string xsintstr="";
  const ProcessInfo & pci = evt->Summary()->ProcInfo();

  if (pci.IsWeakCC()) {
    InteractionType=1;
    xsintstr="tot_cc";
  }
  else if (pci.IsWeakNC()){
    InteractionType=2;
    xsintstr="tot_nc";
  }
    
  double InteractionCrosssection = 0;
  if (InteractionType>0)
    // convert units of 1e-38 cm2 to mb (=1e-31 m2)
    InteractionCrosssection = CrossSectionAccessor::GetMe()->GetCrossSection(xsintstr.c_str(),nupdg,tgtpdg,PrimaryNeutrinoEnergy) * 1e-11;

  weightdict_["GeneratorVolume"]            = gOptGenVolRadius*gOptGenVolRadius*TMath::Pi()*gOptGenVolLength;
  weightdict_["TargetPDG"]                  = tgtpdg;
  weightdict_["TotalCrosssection"]          = TotalCrosssection;
  weightdict_["TotalInteractionProbabilityWeight"] = probability;
  weightdict_["InteractionColumnDepth"]     = InteractionColumnDepth;
  weightdict_["TotalColumnDepth"]           = TotalColumnDepth;
  weightdict_["TotalDetectionLength"]       = gOptGenVolLength;
  weightdict_["LengthInVolume"]             = LengthInVolume;
  weightdict_["EnergyLost"]                 = 0.0;

  weightdict_["InteractionType"]            = InteractionType;
  weightdict_["InteractionCrosssection"]    = InteractionCrosssection;
  //

  // generation area is a circle with radius InjectionSurfaceR
  // in cm^2 !!!! Flux is in cm^-2
  double areaNorm = gOptGenVolRadius*gOptGenVolRadius*TMath::Pi()*1e4;

  // generation solid angle
  double solidAngle = (cos(gOptNuZenithMin*TMath::DegToRad())-cos(gOptNuZenithMax*TMath::DegToRad()));
  solidAngle *= (gOptNuAzimuthMax*TMath::DegToRad()-gOptNuAzimuthMin*TMath::DegToRad());

  double energyIntegral=0;
  if (gOptPowerLawIndex == 1) {
    // if E^-1 then integral over Emin and Emax is
    energyIntegral = TMath::Log(1+gOptNuEnergyRange/gOptNuEnergy);
  } else {
    // if not E^-1 then integral over Emin and Emax is
    energyIntegral = (TMath::Power(gOptNuEnergy+gOptNuEnergyRange, (1.-gOptPowerLawIndex)) -
		      TMath::Power(gOptNuEnergy, (1.-gOptPowerLawIndex))) / (1.-gOptPowerLawIndex);
  }

  //power law index is probabaly 1 indicating a spectrum generated at E^-1
  double energyFactor = TMath::Power( PrimaryNeutrinoEnergy , -gOptPowerLawIndex );

  // OneWeight
  double OneWeight = (probability / energyFactor) * (energyIntegral * areaNorm * solidAngle);

  weightdict_["PrimaryNeutrinoEnergy"] = PrimaryNeutrinoEnergy;
  weightdict_["OneWeight"] = OneWeight;
  weightdict_["NEvents"] = gOptNevents;

}

//________________________________________________________________________________
double CalculateColumnDepth(double length)
{
  double coldepth = ICE_DENSITY * length; // is in g/cm^3 * m, multiply with 1e6 to get g/m^2
  return coldepth*1e6;
}

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
