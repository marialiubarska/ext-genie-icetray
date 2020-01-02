// Create a simple IceCube geometry and write out as root file
// M. Dierckxsens - 2010/05/31

void icecube_geo(Bool_t writeout=false){

  gSystem->Load("libGeom");

  // create the geometry manager
  TGeoManager *geom = new TGeoManager("IceCube","Simple IceCube geometry");

  // Get a few elements from the table
  TGeoElementTable *table = geom->GetElementTable();
  TGeoElement *H_elem = table->GetElement(1);
  TGeoElement *N_elem = table->GetElement(7);
  TGeoElement *O_elem = table->GetElement(8);
  TGeoElement *Si_elem = table->GetElement(14);

  // define some materials
  
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum",0,0,0);
  
  Int_t natoms=0;

  // -- Air 70% N 30% O, density = 1.29e-3 g/cm^3
  TGeoMixture *mixAir = new TGeoMixture("Air",2,1.29e-3);
  mixAir->AddElement(N_elem,0.70);  
  mixAir->AddElement(O_elem,0.30);
    
  // -- Ice H2O, density = 0.93 g/cm^3
  TGeoMixture *mixIce = new TGeoMixture("Ice",2,0.93);
  mixIce->AddElement(H_elem,natoms=2);
  mixIce->AddElement(O_elem,natoms=1);

  // -- Rock SiO2, density = 2.7g/cm^3
  TGeoMixture *mixRock = new TGeoMixture("Rock",3,2.7);
  mixRock->AddElement(O_elem,natoms=2);
  mixRock->AddElement(Si_elem,natoms=1);

  // define media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",0, matVacuum);
  TGeoMedium *Ice = new TGeoMedium("Ice",1, mixIce);
  TGeoMedium *Rock = new TGeoMedium("Rock",2, mixRock);
  TGeoMedium *Air = new TGeoMedium("Air",1, mixAir);


  // make the top container volume
  Double_t v_top_height = 2.5e5; // half height of top volume box 2.5km (in cm)
  Double_t v_top_length = 2.5e5; // half length of top volume box 2.5km (in cm)
  Double_t v_top_width  = 2.5e5; // half width of top volume box 2.5km (in cm)

  TGeoVolume *top = geom->MakeBox("TOP",Vacuum,v_top_width,v_top_length,v_top_height); 
  geom->SetTopVolume(top);
  geom->SetTopVisible(0);

  // determine the translation of the ice, air and rock boxes, wrt center
  // the driving factor here is the height of the ice
  Double_t ice_depth        = 2.820e5; // depth of the ice 
  Double_t v_ice_center_top = 1.95e5; // distance from detector center to top of ice

  Double_t v_ice_center_bot = ice_depth-v_ice_center_top; // distance from center to bottom of ice
  Double_t v_ice_height     = ice_depth/2; // half height of ice volume
  Double_t ice_center_z     = v_ice_height-v_ice_center_bot; // center position of ice volume

  Double_t v_rock_height    = (v_top_height-v_ice_center_bot)/2;  // half height of rock volume
  Double_t rock_center_z    = -1*v_ice_center_bot-v_rock_height;  // center position of rock volume

  Double_t v_air_height    = (v_top_height-v_ice_center_top)/2;  // half height of air volume
  Double_t air_center_z    = v_ice_center_top+v_air_height;  // center position of air volume

  TGeoTranslation *tr_ice = new TGeoTranslation(0.,0.,ice_center_z);
  TGeoTranslation *tr_rock = new TGeoTranslation(0.,0.,rock_center_z);
  TGeoTranslation *tr_air = new TGeoTranslation(0.,0.,air_center_z);

  //make & add the ice layer  
  TGeoVolume *ice_vol = geom->MakeBox("ice_vol",Ice,v_top_width,v_top_length,v_ice_height);
  ice_vol->SetFillColor(9);
  ice_vol->SetLineColor(9);
  top->AddNode(ice_vol,1,tr_ice);

  // make & add the rock layer
  TGeoVolume *rock_vol = geom->MakeBox("rock_vol",Rock,v_top_width,v_top_length,v_rock_height);
  rock_vol->SetFillColor(28);
  rock_vol->SetLineColor(28);
  top->AddNode(rock_vol,1,tr_rock);

  // make & add the rock layer
  TGeoVolume *air_vol = geom->MakeBox("air_vol",Air,v_top_width,v_top_length,v_air_height);
  air_vol->SetFillColor(7);
  air_vol->SetLineColor(7);
  top->AddNode(air_vol,1,tr_air);

  //top->SetVisibility(0);
  geom->CloseGeometry();

  //geom->SetVisLevel(1);
  top->Draw("ogl");

  // check for overlaps
  

  // write out geometry if requested
  if (writeout) {
    geom->Export("IceCubeGeometry.root");
  }

}
