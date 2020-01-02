//____________________________________________________________________________
/*!

\program eventloop

\brief   Simple event loop over GENIE events making a few plots

         Histograms are saved into the file eventloop.root

         Syntax:
           eventloop -f input [-n nevt]

         Options
           [] Denotes an optional argument
           -f input file name
           -n number of events to process [defaut is all]

\author  Mark Dierckxsens - ULB

\created January 8 2010

*/
//____________________________________________________________________________

#include <string>
#include <iostream>

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
#include "Viewer/GHepPrinter.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

using std::string;
using namespace genie;
using namespace std;

void GetCommandLineArgs (int argc, char ** argv);
void BookHistos();
TH1D* GetHisto(const char* hn);
TH2D* GetHisto2D(const char* hn);
void WriteHistos();

int    gOptNEvt;
string gOptInpFilename;
TObjArray gHlist(0);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  LOG("myAnalysis", pDEBUG) << "Booking histograms";
  BookHistos();

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;


  LOG("myAnalysis", pDEBUG) << "Open file for reading...";
  TFile file(gOptInpFilename.c_str(),"READ");
  LOG("myAnalysis", pDEBUG) << "   ... opened";

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) {
    LOG("myAnalysis", pERROR) << "Could not find gtree in input file";
    return 1;
  }

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  cout << "Looping over " << nev << " events" << endl;
  //
  // Loop over all events
  //
  for(int i = 0; i < nev; i++) {

    //cout << "Event : " << i << endl;
    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    event.Print(cout);

    Interaction* in = event.Summary();
    const ProcessInfo & proc = in->ProcInfo();
    //const InitialState & inst = in->InitState();
    //const Kinematics & kine = in->Kine();

    GHepParticle* nu = event.Probe();
    GHepParticle* fsl = event.FinalStatePrimaryLepton();

    const TLorentzVector& k1 = *(nu->P4());
    const TLorentzVector& k2 = *(fsl->P4());

    const TLorentzVector q = k1 - k2;
    
    double v = q.Energy();
    double yinel = v / k1.Energy();
    if (proc.IsWeakCC()) 
      GetHisto("h_yinel_cc")->Fill(yinel);
    else if (proc.IsWeakNC()) 
      GetHisto("h_yinel_nc")->Fill(yinel);
    else
      LOG("myAnalysis", pWARN)  << "Unclassified interactions type: " \
                                << proc.InteractionTypeAsString();

    GetHisto("h_nu_e_all")->Fill(k1.Vect().Mag());
    GetHisto("h_nu_theta_all")->Fill(k1.Vect().Theta()*TMath::RadToDeg());
    GetHisto("h_nu_costheta_all")->Fill(TMath::Cos(k1.Vect().Theta()));
    GetHisto("h_nu_phi_all")->Fill(k1.Vect().Phi()*TMath::RadToDeg());

    double Xvtx = event.Vertex()->X();
    double Yvtx = event.Vertex()->Y();
    double Rvtx = TMath::Sqrt(Xvtx*Xvtx+Yvtx*Yvtx);
    GetHisto("h_nu_x_all")->Fill(Xvtx);
    GetHisto("h_nu_y_all")->Fill(Yvtx);
    GetHisto("h_nu_r_all")->Fill(Rvtx);
    GetHisto("h_nu_z_all")->Fill(event.Vertex()->Z());
    GetHisto2D("h_nu_x_y_all")->Fill(Xvtx,Yvtx);
    GetHisto2D("h_nu_r_z_all")->Fill(Rvtx,event.Vertex()->Z())
;

    // only CC events
    if (!proc.IsWeakCC()) continue;
    //if (i<10) LOG("myAnalysis", pNOTICE) << event;


    //
    // Loop over all particles in this event
    //

    GHepParticle * p = 0;
    TIter event_iter(&event);

    TLorentzVector hadsum(0,0,0,0);
    Int_t nhadr=0;
    //cout << endl;
    GHepParticle* fs_had = event.FinalStateHadronicSystem();
    TLorentzVector hadsys(0,0,0,0);
    if (fs_had) hadsys = *(fs_had->P4());
    vector<GHepParticle*> hadvec;
    while((p=dynamic_cast<GHepParticle *>(event_iter.Next())))
    {
      if ( !(p->Status() == kIStStableFinalState)) continue;
      if (TMath::Abs(p->Pdg()) >100       // hadronic particles
          || (p->Pdg()==22 &&    //gammas from hadronic particles
              TMath::Abs(event.Particle(p->FirstMother())->Pdg())>100)){     
        //cout << *p;
        nhadr++ ;
        hadsum += *(p->P4());
        hadvec.push_back(p);
      }
    }// end loop over particles	
    //cout << endl;

    Double_t p2_tot=0,p_tot=0,p_tot_back=0;    
    TMatrixDSym SpherTen(3);
    vector<GHepParticle*>::iterator itervec;
    for (itervec=hadvec.begin(); itervec != hadvec.end(); itervec++){
      // fill unnorm spher tensor
      TLorentzVector vv = *((*itervec)->P4());
      p2_tot += vv.Vect().Mag2();
      p_tot += vv.Vect().Mag();
      Double_t pangle = hadsum.Angle(vv.Vect());
      GetHisto("h_pangle")->Fill(TMath::RadToDeg()*pangle);
      GetHisto2D("h_pangle_e")->Fill(k1.Energy(),TMath::RadToDeg()*pangle);
      GetHisto2D("h_pangle_ppart")->Fill(vv.Vect().Mag(),TMath::RadToDeg()*pangle);
      GetHisto2D("h_phad_angle")->Fill(hadsum.P(),TMath::RadToDeg()*pangle);
      GetHisto2D("h_phad_angle_he")->Fill(hadsum.P(),TMath::RadToDeg()*pangle);
      if (k1.Energy()<10) 
        GetHisto2D("h_phad_angle_10gev")->Fill(hadsum.P(),TMath::RadToDeg()*pangle);
      if (TMath::Cos(pangle)<0) p_tot_back += vv.Vect().Mag();
      for (Int_t ii=0; ii<3; ++ii){
        for (Int_t jj=0; jj<3 ; ++jj){
         SpherTen(ii,jj) += vv(ii)*vv(jj);
        }
      }
    }

    if (p_tot==0) cout << "Warning 0 total momentum" << endl;
    else {
      GetHisto("h_pfrac_back")->Fill(p_tot_back/p_tot);
      GetHisto2D("h_pfrac_back_e")->Fill(k1.Energy(),p_tot_back/p_tot);
    }
    //cout << "n hadrons = : " << nhadr << " vec length = " << hadvec.size() 
    //     << " P_tot = " << p2_tot << endl;
    // norm spher tensor
    SpherTen *= (1./p2_tot);

    //for (Int_t ii=0;ii<3;++ii)
    //  cout << SpherTen(ii,0) << "\t" << SpherTen(ii,1) 
    //       << "\t" << SpherTen(ii,2) << endl;

    TMatrixDSymEigen SpherTenEigen(SpherTen);
    TVectorD SpherEigenVals = SpherTenEigen.GetEigenValues();
    TMatrixD EigenVecs = SpherTenEigen.GetEigenVectors();
    TMatrixDColumn tdc_evt_axis(EigenVecs,0);
    TVector3 evt_axis(tdc_evt_axis(0),tdc_evt_axis(1),tdc_evt_axis(2));

    GetHisto("h_angle_axis_sum")->Fill(TMath::RadToDeg()*evt_axis.Angle(hadsum.Vect()));


    //cout << "Eigenvalues: ";
    //for (Int_t k = SpherEigenVals.GetLwb(); k <= SpherEigenVals.GetUpb(); k++)
    //  cout << " " << SpherEigenVals(k);
    //cout << endl;

    Double_t sphericity = 3/2*(1-SpherEigenVals(0));
    Double_t aplanarity = 3/2*SpherEigenVals(2);

    GetHisto("h_nu_e")->Fill(k1.Energy());
    GetHisto("h_hadron_n")->Fill(nhadr);
    GetHisto2D("h_hadn_e")->Fill(k1.Energy(),nhadr);

    GetHisto("h_hsum_px")->Fill(hadsum.Px());
    GetHisto("h_hsum_py")->Fill(hadsum.Py());
    GetHisto("h_hsum_pz")->Fill(hadsum.Pz());
    GetHisto("h_hsum_ptot")->Fill(hadsum.P());
    GetHisto("h_hsum_e")->Fill(hadsum.E());

    GetHisto("h_hsys_px")->Fill(hadsys.Px());
    GetHisto("h_hsys_py")->Fill(hadsys.Py());
    GetHisto("h_hsys_pz")->Fill(hadsys.Pz());
    GetHisto("h_hsys_ptot")->Fill(hadsys.P());
    GetHisto("h_hsys_e")->Fill(hadsys.E());

    GetHisto("h_hdif_px")->Fill(hadsum.Px()-hadsys.Px());
    GetHisto("h_hdif_py")->Fill(hadsum.Py()-hadsys.Py());
    GetHisto("h_hdif_pz")->Fill(hadsum.Pz()-hadsys.Pz());
    GetHisto("h_hdif_ptot")->Fill(hadsum.P()-hadsys.P());
    Double_t edif=hadsum.E()-hadsys.E();
    GetHisto("h_hdif_e")->Fill(edif);

    //    if (TMath::Abs(edif) > 8 ){
    //  LOG("myAnalysis", pNOTICE) << "Large energy difference: " << edif;
    //  LOG("myAnalysis", pNOTICE) << "# hadrons: " << nhadr;
    //  for (itervec=hadvec.begin(); itervec != hadvec.end(); itervec++){
    //    LOG("myAnalysis", pNOTICE) << "   Pdg: " << (*itervec)->Pdg()
    //                               << "   E: " << (*itervec)->E();
    //  }
      //LOG("myAnalysis", pNOTICE) << event;
    //}

    GetHisto("h_sphericity")->Fill(sphericity);
    if (k1.Energy()<10) GetHisto("h_sphericity_10gev")->Fill(sphericity);
    GetHisto2D("h_spher_e")->Fill(k1.Energy(),sphericity);
    GetHisto2D("h_spher_phad")->Fill(hadsum.P(),sphericity);
    GetHisto2D("h_spher_nhad")->Fill(nhadr,sphericity);

    GetHisto("h_aplanarity")->Fill(aplanarity);
    GetHisto2D("h_aplan_e")->Fill(k1.Energy(),aplanarity);
    GetHisto2D("h_aplan_phad")->Fill(hadsum.P(),aplanarity);
    GetHisto2D("h_aplan_nhad")->Fill(nhadr,aplanarity);


    // clear current mc event record
    mcrec->Clear();

  }//end loop over events

  // close input GHEP event file
  file.Close();

  WriteHistos();

  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("myAnalysis", pINFO) << "Parsing commad line arguments";

  // get GENIE event sample
  try {
    LOG("myAnalysis", pINFO) << "Reading event sample filename";
    gOptInpFilename = utils::clap::CmdLineArgAsString(argc,argv,'f');
    LOG("myAnalysis", pINFO) << "   " << gOptInpFilename;
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("myAnalysis", pFATAL) 
        << "Unspecified input filename - Exiting";
      exit(1);
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
//_________________________________________________________________________________


void BookHistos()
{  
  TH1D* thd= 0;
  TH2D* thdd=0;

  thd = new TH1D("h_yinel_nc","y distribution - NC",100,0,1);
  gHlist.Add(thd);

  thd = new TH1D("h_yinel_cc","y distribution - CC",100,0,1);
  gHlist.Add(thd);

  thd = new TH1D("h_nu_e_all","Neutrino Energy - all",500,0,500);
  gHlist.Add(thd);
  
  thd = new TH1D("h_nu_theta_all","Neutrino theta - all",90,0,180);
  gHlist.Add(thd);
  
  thd = new TH1D("h_nu_costheta_all","Neutrino cos theta - all",100,-1,1);
  gHlist.Add(thd);
  
  thd = new TH1D("h_nu_phi_all","Neutrino phi - all",90,-180,180);
  gHlist.Add(thd);
  
  thd = new TH1D("h_nu_x_all","Neutrino x position - all",260,-1300,1300);
  gHlist.Add(thd);
  thd = new TH1D("h_nu_y_all","Neutrino y position - all",260,-1300,1300);
  gHlist.Add(thd);
  thd = new TH1D("h_nu_r_all","Neutrino radial position - all",260,0,1300);
  gHlist.Add(thd);
  thd = new TH1D("h_nu_z_all","Neutrino z position - all",220,-1100,1100);
  gHlist.Add(thd);

  thdd = new TH2D("h_nu_x_y_all","Neutrino y vs x",54,-1400,1400,54,-1400,1400);
  gHlist.Add(thdd);
  thdd = new TH2D("h_nu_r_z_all","Neutrino z vs R",54,0,1400,44,-1100,1100);
  gHlist.Add(thdd);

  thd = new TH1D("h_nu_e","Neutrino Energy",400,0,100);
  gHlist.Add(thd);

  thd = new TH1D("h_hadron_n","Hadron multiplicity",20,0,20);
  gHlist.Add(thd);

  thd = new TH1D("h_hsum_px","Px of hadronic sum",200,0,20);
  gHlist.Add(thd);

  thd = new TH1D("h_hsum_py","Py of hadronic sum",200,0,20);
  gHlist.Add(thd);

  thd = new TH1D("h_hsum_pz","Pz of hadronic sum",300,0,30);
  gHlist.Add(thd);

  thd = new TH1D("h_hsum_ptot","Ptot of hadronic sum",300,0,30);
  gHlist.Add(thd);

  thd = new TH1D("h_hsum_e","E of hadronic sum",400,0,40);
  gHlist.Add(thd);

  thd = new TH1D("h_hsys_px","Px of hadronic system",200,0,20);
  gHlist.Add(thd);

  thd = new TH1D("h_hsys_py","Py of hadronic sysytem",200,0,20);
  gHlist.Add(thd);

  thd = new TH1D("h_hsys_pz","Pz of hadronic sysytem",300,0,30);
  gHlist.Add(thd);

  thd = new TH1D("h_hsys_ptot","Ptot of hadronic system",300,0,30);
  gHlist.Add(thd);

  thd = new TH1D("h_hsys_e","E of hadronic system",400,0,40);
  gHlist.Add(thd);

  thd = new TH1D("h_hdif_px","Px diff of hadronic system",200,-10,10);
  gHlist.Add(thd);

  thd = new TH1D("h_hdif_py","Py diff of hadronic sysytem",200,-10,10);
  gHlist.Add(thd);

  thd = new TH1D("h_hdif_pz","Pz diff of hadronic sysytem",300,-15,15);
  gHlist.Add(thd);

  thd = new TH1D("h_hdif_ptot","Ptot diff of hadronic system",400,-20,20);
  gHlist.Add(thd);

  thd = new TH1D("h_hdif_e","E diff of hadronic system",200,-10,10);
  gHlist.Add(thd);

  thdd = new TH2D("h_hadn_e","Hadron multiplicity vs Energy",50,0,100,20,0,20);
  gHlist.Add(thdd);

  thd = new TH1D("h_sphericity","Sphericity of hadronic system",200,0,1);
  gHlist.Add(thd);
  thd = new TH1D("h_sphericity_10gev","Sphericity of had sys, Enu<10GeV",200,0,1);
  gHlist.Add(thd);
  thdd = new TH2D("h_spher_e","Sphericity vs E",50,0,100,100,0,1);
  gHlist.Add(thdd);
  thdd = new TH2D("h_spher_phad","Sphericity vs Ptot hadronic",40,0,20,100,0,1);
  gHlist.Add(thdd);
  thdd = new TH2D("h_spher_nhad","Sphericity vs N hadrons",20,0,20,100,0,1);
  gHlist.Add(thdd);
  thd = new TH1D("h_aplanarity","Aplanarity of hadronic system",200,0,0.5);
  gHlist.Add(thd);
  thdd = new TH2D("h_aplan_e","Aplanarity vs E",50,0,100,100,0,0.5);
  gHlist.Add(thdd);
  thdd = new TH2D("h_aplan_phad","Aplanarity vs Ptot hadronic",40,0,20,100,0,0.5);
  gHlist.Add(thdd);
  thdd = new TH2D("h_aplan_nhad","Aplanarity vs N hadrons",20,0,20,100,0,0.5);
  gHlist.Add(thdd);

  thd = new TH1D("h_angle_axis_sum","Angle sphericity axis - summed vector",
                 180,0,180);
  gHlist.Add(thd);

  thd = new TH1D("h_pangle","Angles of particles wrt summed system",180,0,180);
  gHlist.Add(thd);
  
  thdd = new TH2D("h_pangle_e","Angles of particles wrt summed system vs E",50,0,100,90,0,180);
  gHlist.Add(thdd);
  
  thdd = new TH2D("h_pangle_ppart","Angles of particles wrt summed system vs Ptot of particle",50,0,25,90,0,180);
  gHlist.Add(thdd);
  
  thd = new TH1D("h_pfrac_back","Fraction of momentum going backwards",100,0,1);
  gHlist.Add(thd);

  thdd = new TH2D("h_pfrac_back_e","Fraction of backwards momentum vs E",50,0,100,50,0,1);
  gHlist.Add(thdd);

  thdd = new TH2D("h_phad_angle","Angle vs momentum",50,0,10,90,0,180);
  gHlist.Add(thdd);

  thdd = new TH2D("h_phad_angle_he","Angle vs momentum",50,0,50,90,0,180);
  gHlist.Add(thdd);

  thdd = new TH2D("h_phad_angle_10gev","Angle vs momentum",50,0,10,90,0,180);
  gHlist.Add(thdd);

}


TH1D* GetHisto(const char* hn)
{
  return (TH1D*)gHlist.FindObject(hn);
}

TH2D* GetHisto2D(const char* hn)
{
  return (TH2D*)gHlist.FindObject(hn);
}

void WriteHistos()
{
  TFile f("eventloop.root","RECREATE");
  gHlist.Write();
  //f.Write();
  f.Close();
}
