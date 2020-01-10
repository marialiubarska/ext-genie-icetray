#include <string>
#include <sstream>

#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>

#include "CrossSectionAccessor.h"

#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"

using std::ostringstream;

using namespace genie;

CrossSectionAccessor* CrossSectionAccessor::csa_ = NULL;


//_____________________________________________________________________________
CrossSectionAccessor::CrossSectionAccessor() : xs_file_(NULL) {
}

//_____________________________________________________________________________
CrossSectionAccessor::~CrossSectionAccessor(){
  if (xs_file_) {
    xs_file_->Close();
    delete xs_file_;
  }
}

//_____________________________________________________________________________
CrossSectionAccessor* CrossSectionAccessor::GetMe(){
  if (!csa_){
    csa_ = new CrossSectionAccessor();
    LOG("CrossSectionAccessor", pDEBUG) << "CrossSectionAccessor singleton created";
  }
  return csa_;
}

//_____________________________________________________________________________
void CrossSectionAccessor::SetCrossSectionFile(const char* fn){
  if (xs_file_ && xs_file_->IsOpen()) {
    LOG("CrossSectionAccessor", pWARN) << "An input file is already open - it will be replaced";  
  }
  csa_->OpenCrossSectionFile(fn);
}

//_____________________________________________________________________________
void CrossSectionAccessor::OpenCrossSectionFile(const char* fn){
  if (xs_file_ && xs_file_->IsOpen()) {
    xs_file_->Close();
    delete xs_file_;
  }
  xs_file_ = new TFile(fn);
  if (!xs_file_->IsOpen()){
    LOG("CrossSectionAccessor", pFATAL) << "Unable to open cross section file ";
  }  
}

//_____________________________________________________________________________
double CrossSectionAccessor::GetCrossSection(const char* proc, int nupdg, int tgtpdg, double E) const 
{
  if (!xs_file_ || !xs_file_->IsOpen()){
    LOG("CrossSectionAccessor", pWARN) << "No valid cross section file opened - return 0 cross section";
    return 0;
  }  

  //-- get pdglibrary for mapping pdg codes to names   
  PDGLibrary * pdglib = PDGLibrary::Instance();

  ostringstream dptr;

  string probe_name = pdglib->Find(nupdg)->GetName();
  string tgt_name   = (tgtpdg==1000000010) ? "n" : pdglib->Find(tgtpdg)->GetName();

  dptr << probe_name << "_" << tgt_name;

  TDirectory* dir = (TDirectory*) xs_file_->Get(dptr.str().c_str());
  TGraph* graph = (TGraph*) dir->Get(proc);

  double xsec=0;
  if (graph){
    xsec=graph->Eval(E);
  }
  else {
    LOG("CrossSectionAccessor", pWARN) << "Could not find TGraph for process " << proc << " - returning 0";
  }
  return xsec;
}

//_____________________________________________________________________________

double CrossSectionAccessor::GetTotalCrossSection(int nupdg, int tgtpdg, double E) const 
{
  double xsec_cc=this->GetCrossSection("tot_cc",nupdg,tgtpdg,E);
  double xsec_nc=this->GetCrossSection("tot_nc",nupdg,tgtpdg,E);

  LOG("CrossSectionAccessor", pDEBUG) << "Total CC cross section = " << xsec_cc;
  LOG("CrossSectionAccessor", pDEBUG) << "Total NC cross section = " << xsec_nc;

  return xsec_cc+xsec_nc;
}

//_____________________________________________________________________________

double  CrossSectionAccessor::GetMassInGram(int ipdg) const
{
  PDGLibrary * pdglib = PDGLibrary::Instance();

  double mass  =  pdglib->Find(ipdg)->Mass(); // mass in GeV
  return mass*1.7827e-24; // mass in gram

}
