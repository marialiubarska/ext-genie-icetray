//____________________________________________________________________________
/*

 Author: Mark Dierckxsens - ULB

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "SimpleIceCubeGeomAnalyzer.h"
#include "EVGDrivers/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
SimpleIceCubeGeomAnalyzer::SimpleIceCubeGeomAnalyzer(int pdg) :
GeomAnalyzerI()
{
  map<int,double> tgtmap;
  tgtmap.insert( map<int, double>::value_type(pdg, 1.) );

  this->Initialize(tgtmap);
}
//___________________________________________________________________________
SimpleIceCubeGeomAnalyzer::SimpleIceCubeGeomAnalyzer(
           unsigned int n, const int tgtpdgc[], const double weight[]) :
GeomAnalyzerI()
{
  map<int,double> tgtmap;
  for(unsigned int i=0; i<n; i++) 
     tgtmap.insert( map<int, double>::value_type(tgtpdgc[i], weight[i]) );

  this->Initialize(tgtmap);
}
//___________________________________________________________________________
SimpleIceCubeGeomAnalyzer::SimpleIceCubeGeomAnalyzer(const map<int,double> & tgtmap) :
GeomAnalyzerI()
{
  this->Initialize(tgtmap);
}
//___________________________________________________________________________
SimpleIceCubeGeomAnalyzer::~SimpleIceCubeGeomAnalyzer()
{
  this->CleanUp();
}
//___________________________________________________________________________
const PDGCodeList & SimpleIceCubeGeomAnalyzer::ListOfTargetNuclei(void)
{
// pdg code list contains a single code corresponding to the material passed
// at the geom analyser ctor

  return *fCurrPDGCodeList;
}
//___________________________________________________________________________
const PathLengthList & SimpleIceCubeGeomAnalyzer::ComputeMaxPathLengths(void)
{
// this is irrelevant for the 'point' geometry - return a path length of 1.
// for the only defined material

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const PathLengthList & SimpleIceCubeGeomAnalyzer::ComputePathLengths(
                  const TLorentzVector & /*x*/, const TLorentzVector & /*p*/)
{
// this is irrelevant for the 'point' geometry - return a path length of 1.
// for the only defined material

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const TVector3 & SimpleIceCubeGeomAnalyzer::GenerateVertex(
  const TLorentzVector & x, const TLorentzVector & /*p*/, int /*tgtpdg*/)
{
// Just return the position given as argument - let the flux generator
// do the job

  fCurrVertex->SetXYZ(x.X(),x.Y(),x.Z());
  return *fCurrVertex;
}
//___________________________________________________________________________
void SimpleIceCubeGeomAnalyzer::Initialize(const map<int,double> & tgtmap)
{
  fCurrVertex = new TVector3(0,0,0);

  fCurrPDGCodeList = new PDGCodeList;
  fCurrPDGCodeList->clear();

  map<int,double>::const_iterator iter;
  for(iter = tgtmap.begin(); iter != tgtmap.end(); ++iter) {
	int tgtpdgc = iter->first;
	fCurrPDGCodeList->push_back(tgtpdgc);
  }

  fCurrPathLengthList = new PathLengthList(tgtmap);

  LOG("SimpleIceCubeGeom", pNOTICE) << *fCurrPDGCodeList;
  LOG("SimpleIceCubeGeom", pNOTICE) << *fCurrPathLengthList;
  LOG("SimpleIceCubeGeom", pDEBUG) << "Done Initialization";
}
//___________________________________________________________________________
void SimpleIceCubeGeomAnalyzer::CleanUp(void)
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________

