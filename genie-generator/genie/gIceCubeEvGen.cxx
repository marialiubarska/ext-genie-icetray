//____________________________________________________________________________
/*!

\program gicecubeevgen

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
#include <TH1.h>
#include <TF1.h>

//GENIE
#include "Conventions/XmlParserStatus.h"
#include "Conventions/GBuild.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
// #include "GHEP/GHepParticle.h"       //added
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

void GetCommandLineArgs (int argc, char ** argv);
void Initialize         (void);
void PrintSyntax        (void);

// create flux and geom drivers
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
void            GenerateEventsUsingFluxOrTgtMix();
GeomAnalyzerI * GeomDriver              (void);
GFluxI *        FluxDriver              (void);
#endif

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

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);
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
  mcj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  mcj_driver->SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();
  if(gOptSingleProbScale) 
	mcj_driver->ForceSingleProbScale();

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


  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gevgen", pNOTICE) << " *** Generating event............ " << ievent;

     // generate a single event for neutrinos coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gevgen", pNOTICE) << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  // Save the generated MC events
  ntpw.Save();

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
    string nue = parser.ArgAsString('e');

    // is it just a value or a range (comma separated set of values)
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);   
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0);
       gOptNuEnergyMin = emin;
       gOptNuEnergyMax = emax;
    } else {
       gOptNuEnergyMin = atof(nue.c_str());
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
    gOptZenithMin = gOptZenithMin/180.0*3.1415926;
    gOptZenithMax = gOptZenithMax/180.0*3.1415926;
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
    gOptAzimuthMin = gOptAzimuthMin/180.0*3.1415926;
    gOptAzimuthMax = gOptAzimuthMax/180.0*3.1415926;
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
  bool using_tgtmix = false;
  if( parser.OptionExists('t') ) {
    LOG("gevgen", pINFO) << "Reading target mix";
    string stgtmix = parser.ArgAsString('t');
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(stgtmix,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
    } else {
      using_tgtmix = true;
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
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
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
    << "\n               -e energy or energy range [GeV] (e.g. 100,1000) " 
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

  // LOG("gevgen", pNOTICE)
  //   << "\n\n" << "Syntax:" << "\n"
  //   << "\n      gevgen [-h]"
  //   << "\n              [-r run#]"
  //   << "\n               -n nev"
  //   << "\n               -e energy (or energy range) "
  //   << "\n               -p neutrino_pdg" 
  //   << "\n               -t target_pdg "
  //   << "\n              [-f flux_description]"
  //   << "\n              [-o outfile_name]"
  //   << "\n              [-w]"
  //   << "\n              [--seed random_number_seed]"
  //   << "\n              [--cross-sections xml_file]"
  //   << "\n              [--event-generator-list list_name]"
  //   << "\n              [--message-thresholds xml_file]"
  //   << "\n              [--unphysical-event-mask mask]"
  //   << "\n              [--event-record-print-level level]"
  //   << "\n              [--mc-job-status-refresh-rate  rate]"
  //   << "\n              [--cache-file root_file]"
  //   << "\n";
}
//____________________________________________________________________________
