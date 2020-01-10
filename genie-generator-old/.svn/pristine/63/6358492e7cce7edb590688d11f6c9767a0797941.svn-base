//____________________________________________________________________________
/*!

\program hepevt_writer

\brief   Convert GENIE output format to hepevt file

         This simple program gets the event details from the GHEP TTree structure
         and converst it to the STDHEP format

         Syntax:
           hepevt_writer [-h] -f input [-o output] [-n nevt]

         Options
           [] Denotes an optional argument
           -h prints out help on using hepevt_writer and exits
           -f input file name
           -o output file name [default is hepevt.dat]
           -n number of events to process [defaut is all]

\author  Mark Dierckxsens - ULB

\created January 8 2010

*/
//____________________________________________________________________________

#include <string>
#include <iostream>
#include <stdio.h>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using namespace genie;
using namespace std;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

int    gOptNEvt;
string gOptInpFilename;
string gOptOutFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;


  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  // open file for writing events in hepevt format
  FILE* ofile;
  ofile = fopen(gOptOutFilename.c_str(),"r");
  if (ofile) {// file exists already, spit out warning and exit
    fclose(ofile);
    LOG("myAnalysis", pFATAL) 
      << "Output file already exists. Delete before continuing.";
    exit(1);
  } 
  else {
    ofile = fopen(gOptOutFilename.c_str(),"w");      
  }

  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    //Interaction* in = event.Summary();
    // const ProcessInfo & proc = in->ProcInfo();
    //const InitialState & inst = in->InitState();

    fprintf(ofile,"%d %d \n",i,event.GetEntriesFast());
    
    //
    // Loop over all particles

    TLorentzVector* vtx = event.Vertex();
    GHepParticle * p = 0;
    TIter event_iter(&event);
    
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      // write out all particles, momentum/energies are in GeV, 
      // get position from vertex, and put it in mm (stdhep standards)
      fprintf(ofile," %d %d %d %d %d %d %e %e %e %e %e %e %e %e %e \n",
              p->Status(),p->Pdg(),
              p->FirstMother(),p->LastMother(),
              p->FirstDaughter(),p->LastDaughter(),
              p->Px(),p->Py(),p->Pz(),p->E(),p->Mass(),
              vtx->X()*1e3,vtx->Y()*1e3,vtx->Z()*1e3,vtx->T());
    }
    
    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  // close input GHEP event file
  file.Close();

  //close output file
  fclose(ofile);
  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("myAnalysis", pINFO) << "Parsing commad line arguments";

  // help?
  bool help = genie::utils::clap::CmdLineArgAsBool(argc,argv,'h');
  if(help) {
    PrintSyntax();
    exit(0);
  }

  // get input filename
  try {
    LOG("myAnalysis", pINFO) << "Reading event sample filename";
    gOptInpFilename = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("myAnalysis", pFATAL) 
        << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // get output filename
  try {
    LOG("myAnalysis", pINFO) << "Reading output filename";
    gOptOutFilename = utils::clap::CmdLineArgAsString(argc,argv,'o');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("myAnalysis", pINFO) 
        << "Unspecified output filename - Will use hepevt.dat";
      gOptOutFilename="hepevt.dat";
    }
  }

  // number of events to analyse
  try {    
    LOG("myAnalysis", pINFO) << "Reading number of events to analyze";
    gOptNEvt = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("myAnalysis", pINFO)
        << "Unspecified number of events to analyze - Use all";
      gOptNEvt = -1;
    }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("hepevt_writer", pWARN)
    << "\n\n" << "Syntax:" << "\n"
    << "   hepect_writer [-h] -f input [-o output] [-n nevt] \n";
}
//____________________________________________________________________________



