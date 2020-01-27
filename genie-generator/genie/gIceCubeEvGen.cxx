//____________________________________________________________________________
/*!

\program gicevgen

\required changes: - add IceCube specific flux and geometry drivers
                   - add output in GHEP format

\brief   To be updated: A simple 'generic' GENIE v+A event generation driver (gevgen).

	 It handles:
 	 a) event generation for a fixed init state (v+A) at fixed energy, or
         b) event generation for simple fluxes (specified either via some
            functional form, tabular text file or a ROOT histogram) and for 
            simple 'geometries' (a target mix with its corresponding weights)

         See the GENIE manual for other apps handling experiment-specific 
         event generation cases using the outputs of detailed neutrino flux 
         simulations and realistic detector geometry descriptions.

         Syntax :
           gevgen [-h] 
                  [-r run#] 
                   -n nev 
                   -e energy (or energy range) 
                   -p neutrino_pdg 
                   -t target_pdg 
		  [-f flux_description] 
                  [-o outfile_name]
		  [-w] 
                  [--seed random_number_seed] 
                  [--cross-sections xml_file]
                  [--event-generator-list list_name]
                  [--tune genie_tune]
                  [--message-thresholds xml_file]          
                  [--unphysical-event-mask mask]
                  [--event-record-print-level level]
                  [--mc-job-status-refresh-rate  rate]
                  [--cache-file root_file]

         Options :
           [] Denotes an optional argument.
           -h 
              Prints-out help on using gevgen and exits.
           -n 
              Specifies the number of events to generate.
           -r 
              Specifies the MC run number.
           -e 
              Specifies the neutrino energy.
	      If what follows the -e option is a comma separated pair of values
              it will be interpreted as an energy range for the flux specified
              via the -f option (see below).
           -p 
              Specifies the neutrino PDG code.
           -t 
              Specifies the target PDG code (pdg format: 10LZZZAAAI) _or_ a target
              mix (pdg codes with corresponding weights) typed as a comma-separated 
              list of pdg codes with the corresponding weight fractions in brackets, 
              eg code1[fraction1],code2[fraction2],... 
              For example, to use a target mix of 95% O16 and 5% H type: 
              `-t 1000080160[0.95],1000010010[0.05]'.
           -f 
              Specifies the neutrino flux spectrum.
              It can be any of:
	      -- A function:
                 eg ` -f x*x+4*exp(-x)' 
              -- A vector file:
                 The vector file should contain 2 columns corresponding to 
                 energy,flux (see $GENIE/data/flux/ for few examples). 
              -- A 1-D ROOT histogram (TH1D):
                 The general syntax is `-f /full/path/file.root,object_name'
           -o
              Specifies the name of the output file events will be saved in.   
           -w 
              Forces generation of weighted events.
              This option is relevant only if a neutrino flux is specified.
              Note that 'weighted' refers to the selection of the primary
              flux neutrino + target that were forced to interact. A weighting
              scheme for the generated kinematics of individual processes can
              still be in effect if enabled..
              ** Only use that option if you understand what it means **
           --seed
              Random number seed.
           --cross-sections
              Name (incl. full path) of an XML file with pre-computed
              cross-section values used for constructing splines.
           --event-generator-list            
              List of event generators to load in event generation drivers.
              [default: "Default"].
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
              [default: "Default"].
           --message-thresholds           
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
           --unphysical-event-mask       
              Allows users to specify a 16-bit mask to allow certain types of
              unphysical events to be written in the output file.
              [default: all unphysical events are rejected]
           --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().              
           --mc-job-status-refresh-rate   
              Allows users to customize the refresh rate of the status file.
           --cache-file                  
              Allows users to specify a cache file so that the cache can be
              re-used in subsequent MC jobs.

	***  See the User Manual for more details and examples. ***

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 05, 2004

\cpright Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
#include <fenv.h> // for `feenableexcept`
#endif

//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TF1.h>

//GENIE
#include "Conventions/XmlParserStatus.h"
#include "Conventions/GBuild.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"       // added
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"       //added
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
#include "Utils/AppInit.h"
#include "Utils/RunOpt.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/CmdLnArgParser.h"

#include "genie/ICResultDict.h"
//#include "genie/ICWeightDictFunc.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
// internal classes                                                                                 
#include "GenieDriversIceCube/GCylindPowerLawFlux.h"
#include "GenieDriversIceCube/GConstantDensityGeometryAnalyzer.h"
// #include "GenieDriversIceCube/I3GENIESystWeights.h"
#endif
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;

//typedef map<string,double> weightmap_t;

void GetCommandLineArgs (int argc, char ** argv);
void FillCommandLineArgsDict(void);
void Initialize         (void);
void PrintSyntax        (void);

// create flux and geom drivers
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
void            GenerateEventsUsingFluxOrTgtMix();
GeomAnalyzerI * GeomDriver              (void);
GFluxI *        FluxDriver              (void);
#endif

// void     ClearMCWeight(void);
void     CalculateMCWeight(EventRecord * evt, GMCJDriver * mcj_driver, weightmap_t  optdict);
void     FillMCWeight(void);
// double   CalculateColumnDepth(double d);
void     WriteMCWeightToFile(FILE* f,int ievt, EventRecord* evt, weightmap_t wdm);


//Default options (override them using the command line arguments):
Long_t        kDefOptRunNu          = 0;       // default run number
NtpMCFormat_t kDefOptNtpFormat      = kNFGHEP; // ntuple format
int           kDefOptNevents        = 0;       // n-events to generate
double        kDefOptCylinderRadius = 1200.0;  // injection radius
double        kDefOptCylinderLength = 2000.0;  // injection length

//User-specified options:
Long_t          gOptRunNu;          // run number
string          gOptOutFileName;    // Optional outfile name
string          gOptStatFileName;   // Status file name, set if gOptOutFileName was set.
int             gOptNevents;        // n-events to generate
double          gOptNuEnergyMin;    // min neutrino energy
double          gOptNuEnergyMax;    // max neutrino energy
double          gOptZenithMin;      // min zenith
double          gOptZenithMax;      // max zenith
double          gOptAzimuthMin;     // min azimuth
double          gOptAzimuthMax;     // max azimuth
double          gOptCylinderRadius; // injection radius
double          gOptCylinderLength; // injection length
map<int,double> gOptTgtMix;         // target mix (each with its relative weight)
string          gOptFlavorString;   // neutrino flavor (NuE, NuMu or NuTau)
double          gOptNuFraction;     // nu/(nu+nubar) fraction
double          gOptPowerLawIndex;  // flux power law index (gamma)
bool            gOptSingleProbScale;// force single probability scale or not
bool            gOptSystWeights;    // calculate GENIE systematic weights or not
long int        gOptRanSeed;        // random number seed
string          gOptInpXSecFile;    // cross-section splines
// options map
weightmap_t gOptDict_;

// weight dict map                                                                               
weightmap_t weightdict_;
// weight dict contents
double  _InjectionSurfaceR                  = 0.;  // gOptCylinderRadius [m]
double  _PowerLawIndex                      = 0.;  // gOptPowerLawIndex
double  _MinEnergyLog                       = 0.;  // log10( gOptNuEnergyMin/GeV )
double  _MaxEnergyLog                       = 0.;  // log10( gOptNuEnergyMax/GeV )
double  _MinZenith                          = 0.;  // gOptZenithMin [rad]
double  _MaxZenith                          = 0.;  // gOptZenithMax [rad]
double  _MinAzimuth                         = 0.;  // gOptAzimuthMin [rad]
double  _MaxAzimuth                         = 0.;  // gOptAzimuthMax [rad]

int32_t  _TargetPDGCode                     = -1.; // evt.Particle(1)->Pdg()
double   _GENIEWeight                       = -1.; // evt.Weight()
double   _GlobalProbabilityScale            = -1.; // mcj_driver->GlobProbScale()
double   _GeneratorVolume                   = -1.; // _InjectionSurfaceR*_InjectionSurfaceR * pi * _TotalDetectionLength [pi*R^2*L]
double   _Crosssection                      = -1.; // evt.XSec() / (1e-38*units::cm2)
double   _InteractionProbabilityWeight      = -1.; // _GENIEWeight * evt.Probability()
double   _TotalInteractionProbabilityWeight = -1.; // _GENIEWeight * _GlobalProbabilityScale

double   _TotalDetectionLength              = 0.;  // gOptCylinderLength [m]
double   _LengthInVolume                    = 0.;  // _TotalDetectionLength/2. * (u_min + 1.); u_min = int_pos * cyl_end / (half_l*half_l); half_l = _TotalDetectionLength/2.
double   _EnergyLost                        = 0.;  // 0.0 --why??
double   _InteractionType                   = -1.; // CC->1, NC->2, else->0; is_CC = evt.Summary()->ProcInfo().IsWeakCC(); is_NC = evt.Summary()->ProcInfo().IsWeakNC();

double   _PrimaryNeutrinoEnergy             = 0.;  // k1.Energy(); const TLorentzVector & k1 = *(neutrino->P4()); GHepParticle * neutrino = evt.Probe();
double   _OneWeight                         = 0.;  // (_TotalInteractionProbability / energyFactor) * (energyIntegral * areaNorm * solidAngle); see ConvertToMCTree.cxx
double   _NEvents                           = 0.;  // static_cast<double>( gOptNevents ) [] 

// densities in kg/m^3                                                                                               
double ICE_DENSITY=0.93e3;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);
  FillCommandLineArgsDict();
  Initialize();

  // throw on NaNs and Infs...
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  
  //
  // Generate neutrino events
  //
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
	GenerateEventsUsingFluxOrTgtMix();
#else
  LOG("gevgen", pERROR) 
    << "\n   To be able to generate neutrino events from a flux and/or a target mix" 
    << "\n   you need to add the following config options at your GENIE installation:" 
    << "\n   --enable-flux-drivers  --enable-geom-drivers \n" ;
  exit(2);
#endif
  
  return 0;
}
//____________________________________________________________________________
void Initialize()
{
  // Initialization of random number generators, cross-section table, 
  // messenger thresholds, cache file
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, false);

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
}
//____________________________________________________________________________

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
//............................................................................
void GenerateEventsUsingFluxOrTgtMix(void)
{
  // Get flux and geom drivers
  GFluxI *        flux_driver = FluxDriver();
  GeomAnalyzerI * geom_driver = GeomDriver();

  // Create the monte carlo job driver
  GMCJDriver * mcj_driver = new GMCJDriver;
  //  mcj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  //  mcj_driver->SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->UseSplines();
  if(gOptSingleProbScale) 
	mcj_driver->ForceSingleProbScale(); // do not set this if weighted events should be generated
  mcj_driver->KeepOnThrowingFluxNeutrinos(true); // always return a neutrino
  mcj_driver->Configure();

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);

  // If an output file name has been specified... use it
  if (!gOptOutFileName.empty()){
    ntpw.CustomizeFilename(gOptOutFileName);
  }
  ntpw.Initialize();

  // Create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  // If a status file name has been given... use it
  if (!gOptStatFileName.empty()){
    mcjmonitor.CustomizeFilename(gOptStatFileName);
  }

  //-- open output file to write weight_dict 
  ostringstream wdic_filename;
  wdic_filename << gOptOutFileName << ".wdict.dat";

  FILE* wd_ofile;
  wd_ofile = fopen(wdic_filename.str().c_str(),"w");

  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gevgen", pNOTICE) << " *** Generating event............ " << ievent;

     // generate a single event for neutrinos coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gevgen", pNOTICE) << "Generated Event GHEP Record: " << *event;

     // CalculateMCWeight(event, mcj_driver);
     // FillMCWeight();

     CalculateMCWeight(event, mcj_driver, gOptDict_);
     FillMCWeight();
     
     // add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);

     // ntpw.AddEventRecord(ievent, event, weightdict_);
     // mcjmonitor.Update(ievent,event, weightdict_);
     
     WriteMCWeightToFile(wd_ofile,ievent,event,weightdict_);
     
     ievent++;
     delete event;
  }

  // Save the generated MC events
  ntpw.Save();

  //-- close weight dict output file
  fclose(wd_ofile);
								  
  delete flux_driver;
  delete geom_driver;
  delete mcj_driver;;
}
//____________________________________________________________________________
GeomAnalyzerI * GeomDriver(void)
{
// create cylinder with constant density tgt nucleon mix (using internal geometry driver)

  GeomAnalyzerI * geom_driver =
    new geometry::GConstantDensityGeometryAnalyzer(gOptTgtMix,
						   gOptCylinderLength);
  return geom_driver;
}
//____________________________________________________________________________
GFluxI * FluxDriver(void)
{
// create power law flux driver (using internal flux driver)
//
  GFluxI * flux_driver = new flux::GCylindPowerLawFlux(gOptFlavorString,
						       gOptCylinderLength,
						       gOptCylinderRadius,
						       gOptPowerLawIndex,
						       gOptNuEnergyMin,
						       gOptNuEnergyMax,
						       gOptZenithMin,
						       gOptZenithMax,
						       gOptAzimuthMin,
						       gOptAzimuthMax,
						       gOptNuFraction);

  return flux_driver;
}
#endif
//____________________________________________________________________________
//____________________________________________________________________________ 
void FillCommandLineArgsDict(void)
{
  gOptDict_.clear();

  // double and int input values
  gOptDict_["RunNu"] =           static_cast<double>( gOptRunNu ); 
  gOptDict_["Nevents"] =         static_cast<double>( gOptNevents );
  gOptDict_["NuEnergyMin"] =     gOptNuEnergyMin;
  gOptDict_["NuEnergyMax"] =     gOptNuEnergyMax;
  gOptDict_["ZenithMin"] =       gOptZenithMin;
  gOptDict_["ZenithMax"] =       gOptZenithMax;
  gOptDict_["AzimuthMin"] =      gOptAzimuthMin;
  gOptDict_["AzimuthMax"] =      gOptAzimuthMax;
  gOptDict_["CylinderRadius"] =  gOptCylinderRadius;
  gOptDict_["CylinderLength"] =  gOptCylinderLength;
  gOptDict_["NuFraction"] =      gOptNuFraction;
  gOptDict_["PowerLawIndex"] =   gOptPowerLawIndex;

  //other input values
  //  gOptDict_["OutFileName"] =     gOptOutFileName;
  //  gOptDict_["StatFileName"] =    gOptStatFileName;
  //  gOptDict_["TgtMix"] =          gOptTgtMix;
  //  gOptDict_["FlavorString"] =    gOptFlavorStrin;  
  //  gOptDict_["SingleProbScale"] = gOptSingleProbScale;
  //  gOptDict_["SystWeights"] =     gOptSystWeights;
  //  gOptDict_["RanSeed"] =         gOptRanSeed;
  //  gOptDict_["InpXSecFile"] =     gOptInpXSecFile;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen", pINFO) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen", pINFO) << "Unspecified run number - Using default";
    gOptRunNu = kDefOptRunNu;
  }

  // Output file name
  if( parser.OptionExists('o') ) {
    LOG("gevgen", pINFO) << "Reading output file name";
    gOptOutFileName = parser.ArgAsString('o');

    gOptStatFileName = gOptOutFileName;
    // strip the output file format and replace with .status
    if (gOptOutFileName.find_last_of(".") != string::npos)
      gOptStatFileName = 
	gOptStatFileName.substr(0, gOptOutFileName.find_last_of("."));
    gOptStatFileName .append(".status");
  }

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = parser.ArgAsInt('n');
  } else {
    LOG("gevgen", pWARN)
       << "Unspecified number of events to generate - Using default";
    gOptNevents = kDefOptNevents;
  }
  
  // neutrino energy
  if( parser.OptionExists('e') ) {
    LOG("gevgen", pINFO) << "Reading neutrino energy";
    string lgnue = parser.ArgAsString('e');

    // is it just a value or a range (comma separated set of values)
    if(lgnue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> lgerange = utils::str::Split(lgnue, ",");
       assert(lgerange.size() == 2);   
       double lgemin = atof(lgerange[0].c_str());
       double lgemax = atof(lgerange[1].c_str());
       assert(lgemax>lgemin && lgemin>=0);
       gOptNuEnergyMin = TMath::Power(10.,lgemin);
       gOptNuEnergyMax = TMath::Power(10.,lgemax);
    } else {
       gOptNuEnergyMin = TMath::Power(10.,atof(lgnue.c_str()));
       gOptNuEnergyMax = gOptNuEnergyMin;
    }
  } else {
    LOG("gevgen", pFATAL) << "Unspecified neutrino energy - Exiting";
    PrintSyntax();
    exit(1);
  }

  // zenith range
  if( parser.OptionExists('z') ) {
    LOG("gevgen", pINFO) << "Reading zenith angle range";
    string zen = parser.ArgAsString('z');
    
    // is it just a value or a range (comma separated set of values)                                  
    if(zen.find(",") != string::npos) {
       // split the comma separated list                                                              
       vector<string> zenrange = utils::str::Split(zen, ",");
       assert(zenrange.size() == 2);
       double zenmin = atof(zenrange[0].c_str());
       double zenmax = atof(zenrange[1].c_str());
       assert(zenmax>zenmin && zenmin>=0);
       gOptZenithMin = zenmin;
       gOptZenithMax = zenmax;
    } else {
       gOptZenithMin = atof(zen.c_str());
       gOptZenithMax = gOptZenithMin;
    }
    gOptZenithMin = gOptZenithMin*TMath::DegToRad();  ///180.0*3.1415926;
    gOptZenithMax = gOptZenithMax*TMath::DegToRad();
  } else {
    LOG("gevgen", pFATAL) << "Unspecified zenith range - Exiting";
    PrintSyntax();
    exit(1);
  }
  
  // azimuth range
  if( parser.OptionExists('a') ) {
    LOG("gevgen", pINFO) << "Reading azimuth angle range";
    string az = parser.ArgAsString('a');

    // is it just a value or a range (comma separated set of values)                                                                                                                                        
    if(az.find(",") != string::npos) {
       // split the comma separated list                                                                                                                                                                 
       vector<string> azrange = utils::str::Split(az, ",");
       assert(azrange.size() == 2);
       double azmin = atof(azrange[0].c_str());
       double azmax = atof(azrange[1].c_str());
       assert(azmax>azmin && azmin>=0);
       gOptAzimuthMin = azmin;
       gOptAzimuthMax = azmax;
    } else {
       gOptAzimuthMin = atof(az.c_str());
       gOptAzimuthMax = gOptAzimuthMin;
    }
    gOptAzimuthMin = gOptAzimuthMin*TMath::DegToRad();
    gOptAzimuthMax = gOptAzimuthMax*TMath::DegToRad();
  } else {
    LOG("gevgen", pFATAL) << "Unspecified azimuth range - Exiting";
    PrintSyntax();
    exit(1);
  }
  
  // injection cylinder radius 
  if( parser.OptionExists('R') ) {
    LOG("gevgen", pINFO) << "Reading injection cylinder radius";
    gOptCylinderRadius = parser.ArgAsDouble('R');
  } else {
    LOG("gevgen", pWARN)
       << "Unspecified injection radius - Using default";
    gOptCylinderRadius = kDefOptCylinderRadius;
  }
  
  // injection cylinder length
  if( parser.OptionExists('L') ) {
    LOG("gevgen", pINFO) << "Reading injection cylinder length";
    gOptCylinderLength = parser.ArgAsDouble('L');
  } else {
    LOG("gevgen", pWARN)
       << "Unspecified injection length - Using default";
    gOptCylinderLength = kDefOptCylinderLength;
  }

  // target mix (their PDG codes with their corresponding weights)
  if( parser.OptionExists('t') ) {
    LOG("gevgen", pINFO) << "Reading target mix";
    string stgtmix = parser.ArgAsString('t');
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(stgtmix,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt*ICE_DENSITY));
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
         LOG("Main", pNOTICE)
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt*ICE_DENSITY));
      }//tgtmix_iter
    }//>1

  } else {
    LOG("gevgen", pFATAL) << "Unspecified target PDG code - Exiting";
    PrintSyntax();
    exit(1);
  }

  // neutrino flavor
  if( parser.OptionExists("nu-type") ) {
    LOG("gevgen", pINFO) << "Reading neutrino type";
    gOptFlavorString = parser.ArgAsString("nu-type");
  } else {
    LOG("gevgen", pFATAL)
       << "Unspecified neutrino type - Exiting";
    PrintSyntax();
    exit(1);
  }
  
  // nu/(nu+nubar) fraction
  if( parser.OptionExists("nu-fraction") ) {
    LOG("gevgen", pINFO) << "Reading neutrino fraction";
    gOptNuFraction = parser.ArgAsDouble("nu-fraction");
  } else {
    LOG("gevgen", pFATAL)
       << "Unspecified neutrino fraction - Exiting";
    PrintSyntax();
    exit(1);
  }
  
  // power law index for flux (gamma)
  if( parser.OptionExists("gamma") ) {
    LOG("gevgen", pINFO) << "Reading flux power law index";
    gOptPowerLawIndex = parser.ArgAsDouble("gamma");
  } else {
    LOG("gevgen", pFATAL)
       << "Unspecified flux power law index - Exiting";
    PrintSyntax();
    exit(1);
  }

  // force SingleProbScale?
  if( parser.OptionExists("force-singleprob-scale") ) {
    LOG("gevgen", pINFO) << "Forcing single probability scale";
    gOptSingleProbScale = true;
  } else {
    LOG("gevgen", pINFO)
       << "NOT forcing single probability scale";
    gOptSingleProbScale = false;
  }

  // calculate GENIE Systematic weights?
  if( parser.OptionExists("enable-syst-weights") ) {
    LOG("gevgen", pINFO) << "GENIE systematic weights will be calculated";
    gOptSystWeights = true;
  } else {
    LOG("gevgen", pINFO)
       << "GENIE systematic weights will NOT be calculated";
    gOptSystWeights = false;
  }  
  
  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen", pWARN) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  // input cross-section file
  if( parser.OptionExists("cross-sections") ) {
    LOG("gevgen", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("cross-sections");
  } else {
    LOG("gevgen", pWARN) << "Unspecified cross-section file";
    gOptInpXSecFile = "";
  }

  //
  // print-out the command line options
  //
  LOG("gevgen", pNOTICE) 
     << "\n" 
     << utils::print::PrintFramedMesg("gicecubeevgen job configuration:");

  LOG("gevgen", pNOTICE) 
     << "MC Run Number: " << gOptRunNu;

  LOG("gevgen", pNOTICE)
     << "Output file name: " << gOptOutFileName;

  LOG("gevgen", pNOTICE)
     << "Status file name: " << gOptStatFileName;

  LOG("gevgen", pNOTICE)
     << "Number of events requested: " << gOptNevents;

  LOG("gevgen", pNOTICE)                                                                        
     << "Neutrino energy range: ["
     << gOptNuEnergyMin << ", " << gOptNuEnergyMax << "]";

  LOG("gevgen", pNOTICE)
     << "Zenith range: ["
     << gOptZenithMin << ", " << gOptZenithMax << "]";

  LOG("gevgen", pNOTICE)
     << "Azimuth range: ["
     << gOptAzimuthMin << ", " << gOptAzimuthMax << "]";

  LOG("gevgen", pNOTICE)
     << "Injection radius: " << gOptCylinderRadius;

  LOG("gevgen", pNOTICE)
     << "Injection length: " << gOptCylinderLength;

  LOG("gevgen", pNOTICE)
      << "Target code (PDG) & weight fraction (in case of multiple targets): ";
  map<int,double>::const_iterator iter;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
      int    tgtpdgc = iter->first;
      double wgt     = iter->second;
      LOG("gevgen", pNOTICE)
          << " >> " <<  tgtpdgc << " (weight fraction = " << wgt << ")";
  }

  LOG("gevgen", pNOTICE)
     << "Flavor: " << gOptFlavorString;

  LOG("gevgen", pNOTICE)
     << "Nu fraction: " << gOptNuFraction;

  LOG("gevgen", pNOTICE)
     << "Gamma: " << gOptPowerLawIndex;

  LOG("gevgen", pNOTICE)
     << "Force SingleProbScale: " << gOptSingleProbScale;

  LOG("gevgen", pNOTICE)
     << "Calculate GENIE systematic weights: " << gOptSystWeights;

  if(gOptRanSeed != -1) {
     LOG("gevgen", pNOTICE) 
       << "Random number seed: " << gOptRanSeed;
  } else {
     LOG("gevgen", pNOTICE) 
       << "Random number seed was not set, using default";
  }

  if(gOptInpXSecFile.size() > 0) {
     LOG("gevgen", pNOTICE) 
       << "Using cross-section splines read from: " << gOptInpXSecFile;
  } else {
     LOG("gevgen", pNOTICE) 
       << "No input cross-section spline file";
  }

  LOG("gevgen", pNOTICE) << "\n";

  LOG("gevgen", pNOTICE) << *RunOpt::Instance();

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gevgen [-h]"
    << "\n              [-r run number]"
    << "\n              [-o outfile_name]"
    << "\n               -n number of events"
    << "\n               -e log10 of energy or energy range [log10(E/GeV)] (e.g. 1,2) " 
    << "\n               -t target mix (e.g. 1000080160[0.95],1000010010[0.05]) "
    << "\n               -z zenith angle range [deg] (e.g. 0,90) "
    << "\n               -a azimuth angle range [deg] (e.g. 0,360) "
    << "\n               -R injection cylinder radius [m] "
    << "\n               -L injection cylinder length [m] "
    << "\n               --nu-type "
    << "\n               --nu-fraction "
    << "\n               --gamma "
    << "\n              [--force-singleprob-scale] "
    << "\n              [--enable-syst-weights] "
    << "\n              [--seed random_number_seed]"
    << "\n              [--cross-sections xml_file]"
    << "\n              [--event-generator-list list_name]"
    << "\n              [--message-thresholds xml_file]"
    << "\n              [--unphysical-event-mask mask]"
    << "\n              [--event-record-print-level level]"
    << "\n              [--mc-job-status-refresh-rate  rate]"
    << "\n              [--cache-file root_file]"
    << "\n";
}
//____________________________________________________________________________
// weight dict
//____________________________________________________________________________
// void ClearMCWeight(){
//   weightdict_.clear();
// }
//____________________________________________________________________________
void FillMCWeight(void){

  weightdict_.clear();

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
void CalculateMCWeight(EventRecord * evt, GMCJDriver * mcj_driver, weightmap_t  optdict){

  // values taken from input
  _NEvents              = optdict.find("Nevents")->second;
  _PowerLawIndex        = optdict.find("PowerLawIndex")->second;

  double minEnergy = optdict.find("NuEnergyMin")->second;
  double maxEnergy = optdict.find("NuEnergyMax")->second;
  _MinEnergyLog         = TMath::Log10( minEnergy );
  _MaxEnergyLog         = TMath::Log10( maxEnergy );
  _MinZenith            = optdict.find("ZenithMin")->second;
  _MaxZenith            = optdict.find("ZenithMax")->second;
  _MinAzimuth           = optdict.find("AzimuthMin")->second;
  _MaxAzimuth           = optdict.find("AzimuthMax")->second;
  
  // geometry
  _InjectionSurfaceR    = optdict.find("CylinderRadius")->second;
  _TotalDetectionLength = optdict.find("CylinderLength")->second;

  _GeneratorVolume      = _InjectionSurfaceR*_InjectionSurfaceR * _TotalDetectionLength * TMath::Pi();
    //gOptCylinderRadius*gOptCylinderRadius * gOptCylinderLength; // V = pi*R^2*L

  // LengthInVolume -- projection of [beginning of the cylinder, vertex position] line onto neutrino direction (== cylinder axis)
  const double half_l = _TotalDetectionLength/2.;  // half length of cylinder
  
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
  const double areaNorm = TMath::Pi() * _InjectionSurfaceR*_InjectionSurfaceR * 1e4; // [cm^2]
  const double solidAngle = (TMath::Cos(_MinZenith)-TMath::Cos(_MaxZenith))*(_MaxAzimuth-_MinAzimuth);
  
  double energyIntegral=0;
  if (_PowerLawIndex == 1.) {
    // if E^-1 then integral over Emin and Emax is                                                                                                                                                   
    energyIntegral = TMath::Log(maxEnergy/minEnergy);
  } else {
    // if not E^-1 then integral over Emin and Emax is                                                                                                                                               
    energyIntegral = (TMath::Power(maxEnergy, (1.-_PowerLawIndex)) -
		      TMath::Power(minEnergy, (1.-_PowerLawIndex))) / (1.-_PowerLawIndex);
  }
    
  const double energyFactor = TMath::Power( _PrimaryNeutrinoEnergy, -_PowerLawIndex );

  _OneWeight = (_TotalInteractionProbabilityWeight / energyFactor) * (energyIntegral * areaNorm * solidAngle);

}
//____________________________________________________________________________
double CalculateColumnDepth(double length)
{
  double coldepth = ICE_DENSITY * length; // is in g/cm^3 * m, multiply with 1e6 to get g/m^2                        
  return coldepth*1e6;
}
//____________________________________________________________________________
void WriteMCWeightToFile(FILE* ofile,int ievt, EventRecord* evt, weightmap_t  wdict)
{
  // first line in the event is the event number, the number of elements in the weightdict,                          
  // the # entries in the event tree, the pdg of incoming neutrino and the energy                                    
  // the last three are used for consistency checks with the hepevt file                                             
  unsigned int len_wdict = wdict.size();

  fprintf(ofile,"%d %u %d %d %e \n",ievt,len_wdict,evt->GetEntriesFast(),evt->Probe()->Pdg(),evt->Probe()->E());

  // now print all the members of the weightdict map                                                                 
  weightmap_t::iterator it;
  for (it=wdict.begin(); it != wdict.end(); it++){
    fprintf(ofile,"    %s %.9e \n",(*it).first.c_str(),(*it).second);
  }
}
//____________________________________________________________________________
