void modify_graph(TGraph* tg,Double_t scale=1){

  for (Int_t i=0; i<tg->GetN();i++){
    Double_t x,y;
    tg->GetPoint(i,x,y);
    if (x==0) cout << "energy zero" << endl;
    //if (i<10 || i>990) cout << "i = " << i << " x = " << x << " y = " << y << endl; 
    y/=x;
    y*=scale;
    tg->SetPoint(i,x,y);
  }

}

void add_graph(TGraph* tg1,TGraph* tg2, TGraph* tgout){
  Int_t nstep =1000;
  tgout->Set(nstep);
  Double_t xmin=0.1;
  Double_t xmax=100;
  Double_t step=(xmax-xmin)/(nstep-1);
  for (Int_t i=0; i<nstep; ++i){
    Double_t x = xmin + i*step;
    Double_t y = tg1->Eval(x) + tg2->Eval(x);
    tgout->SetPoint(i,x,y);
    //cout << i << " " << x << " " << y << endl;
  }
    
}


void plot_xsec()
{
  TFile f("/user/mdier/scratch/icecube/genie/xsec_splines/gxspl-water-v2.5.1.root");
  
  //TDirectory* dir_nu_mu_O16 = (TDirectory)* f.Get("nu_mu_O16");
  // get total cc xsec
  TGraph* tg_tot_cc;
  f.GetObject("nu_mu_O16/tot_cc",tg_tot_cc);
  // get total nc xsec
  TGraph* tg_tot_nc;
  f.GetObject("nu_mu_O16/tot_nc",tg_tot_nc);
  // get cc-coh xsec
  TGraph* tg_coh_cc;
  f.GetObject("nu_mu_O16/coh_cc",tg_coh_cc);
  //
  TGraph* tg_dis_cc_n;
  f.GetObject("nu_mu_O16/dis_cc_n",tg_dis_cc_n);
  //
  TGraph* tg_dis_cc_p;
  f.GetObject("nu_mu_O16/dis_cc_p",tg_dis_cc_p);
  //
  TGraph* tg_qel_cc_n;
  f.GetObject("nu_mu_O16/qel_cc_n",tg_qel_cc_n);

  TGraph* tg_qel_cc_p;
  f.GetObject("nu_mu_O16/qel_cc_p_charm4222",tg_qel_cc_p);
  
  TGraph* tg_res_cc_n;
  f.GetObject("nu_mu_O16/res_cc_n",tg_res_cc_n);

  TGraph* tg_res_cc_p;
  f.GetObject("nu_mu_O16/res_cc_p",tg_res_cc_p);


  modify_graph(tg_tot_cc);
  modify_graph(tg_tot_nc);
  modify_graph(tg_coh_cc);
  modify_graph(tg_dis_cc_n);
  modify_graph(tg_dis_cc_p);
  modify_graph(tg_qel_cc_n);
  modify_graph(tg_qel_cc_p);
  modify_graph(tg_res_cc_n);
  modify_graph(tg_res_cc_p);

  TGraph* tg_dis_cc = new TGraph();
  add_graph(tg_dis_cc_n,tg_dis_cc_p,tg_dis_cc);

  TGraph* tg_qel_cc = new TGraph();
  add_graph(tg_qel_cc_n,tg_qel_cc_p,tg_qel_cc);

  TGraph* tg_res_cc = new TGraph();
  add_graph(tg_res_cc_n,tg_res_cc_p,tg_res_cc);

  TCanvas* tc = new TCanvas("tc","canvas",700,500);

  tc->SetLogy();
  gStyle->SetOptStat(0);
  

  TH1D* tdum= new TH1D("tdum","#nu_{#mu} CC Cross Sections;E [GeV];#sigma/E [10^{-38} cm^{2}/GeV]",10,0,30);
  tdum->SetMaximum(100);
  tdum->SetMinimum(0.001);

  tdum->Draw();

  tg_tot_cc->SetLineWidth(2);
  tg_tot_cc->Draw("C");

  tg_coh_cc->SetLineWidth(2);
  tg_coh_cc->SetLineColor(3);
  tg_coh_cc->Draw("C");

  tg_dis_cc->SetLineWidth(2);
  tg_dis_cc->SetLineColor(2);
  tg_dis_cc->Draw("C");

  //tg_dis_cc_n->SetLineWidth(2);
  //tg_dis_cc_n->SetLineStyle(2);
  //tg_dis_cc_n->SetLineColor(3);
  //tg_dis_cc_n->Draw("C");

  //tg_dis_cc_p->SetLineWidth(2);
  //tg_dis_cc_p->SetLineStyle(3);
  //tg_dis_cc_p->SetLineColor(3);
  //tg_dis_cc_p->Draw("C");

  tg_qel_cc->SetLineWidth(2);
  tg_qel_cc->SetLineColor(4);
  tg_qel_cc->Draw("C");

  //tg_qel_cc_n->SetLineWidth(2);
  //tg_qel_cc_n->SetLineStyle(2);
  //tg_qel_cc_n->SetLineColor(4);
  //tg_qel_cc_n->Draw("C");

  //tg_qel_cc_p->SetLineWidth(2);
  //tg_qel_cc_p->SetLineStyle(3);
  //tg_qel_cc_p->SetLineColor(4);
  //tg_qel_cc_p->Draw("C");

  tg_res_cc->SetLineWidth(2);
  tg_res_cc->SetLineColor(6);
  tg_res_cc->Draw("C");

  //  tg_res_cc_n->SetLineWidth(2);
  //tg_res_cc_n->SetLineColor(6);
  //tg_res_cc_n->Draw("C");

  //tg_res_cc_p->SetLineWidth(2);
  //tg_res_cc_p->SetLineStyle(2);
  //tg_res_cc_p->SetLineColor(6);
  //tg_res_cc_p->Draw("C");

  tc->Print("xsec_numu.png");
}

