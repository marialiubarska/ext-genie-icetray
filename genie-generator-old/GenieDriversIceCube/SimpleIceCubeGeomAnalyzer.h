//____________________________________________________________________________
/*!

\class   genie::geometry::SimpleIceCubeGeomAnalyzer

\brief   Based on the PointGeomAnalyzer, it will return the vertex postion
         given as argument to GenerateVertex.


\author  Mark Dierckxsens

\created March 15, 2010

*/
//____________________________________________________________________________

#ifndef _SIMPLE_ICECUBE_GEOMETRY_ANALYZER_H_
#define _SIMPLE_ICECUBE_GEOMETRY_ANALYZER_H_

#include <map>

#include "EVGDrivers/GeomAnalyzerI.h"

using std::map;

namespace genie    {
namespace geometry {

class SimpleIceCubeGeomAnalyzer : public GeomAnalyzerI {

public :
  SimpleIceCubeGeomAnalyzer(int tgtpdgc);
  SimpleIceCubeGeomAnalyzer(unsigned int n, const int tgt_pdg[], const double weight[]);
  SimpleIceCubeGeomAnalyzer(const map<int,double> & tgtmap /* pdg -> weight*/);
 ~SimpleIceCubeGeomAnalyzer();

  // implement the GeomAnalyzerI interface

  const PDGCodeList &    ListOfTargetNuclei    (void);
  const PathLengthList & ComputeMaxPathLengths (void);

  const PathLengthList &
           ComputePathLengths
             (const TLorentzVector & x, const TLorentzVector & p);
  const TVector3 &
           GenerateVertex
             (const TLorentzVector & x, const TLorentzVector & p, int tgtpdg);
private:

  void Initialize (const map<int,double> & tgtmap);
  void CleanUp    (void);

  TVector3 *       fCurrVertex;          ///< current generated vertex
  PathLengthList * fCurrPathLengthList;  ///< current list of path-lengths
  PDGCodeList *    fCurrPDGCodeList;     ///< current list of target nuclei
};

}      // geometry namespace
}      // genie    namespace

#endif // _SIMPLE_ICECUBE_GEOMETRY_ANALYZER_H_
