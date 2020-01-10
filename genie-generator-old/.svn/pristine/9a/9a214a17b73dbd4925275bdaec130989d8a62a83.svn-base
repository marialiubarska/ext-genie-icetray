//____________________________________________________________________________
/*!

\class   CrossSectionAccessor

\brief   Get the total cross section from the ROOT file

\author  Mark Dierckxsens

\created October 8, 2010

*/
//____________________________________________________________________________

#ifndef _CROSS_SECTION_ACCESSOR_H
#define _CROSS_SECTION_ACCESSOR_H

class TFile;

class CrossSectionAccessor{

public:

  static CrossSectionAccessor* GetMe();
  ~CrossSectionAccessor();

  void    SetCrossSectionFile(const char* fn);
  double  GetCrossSection(const char* proc, int nupdg, int tgtpdg, double E) const;
  double  GetTotalCrossSection(int nupdg, int tgtpdg, double E) const;
  double  GetMassInGram(int ipdg) const;

private:

  /// Return total NC + CC cross section in units of 1e-38 cm^2
  CrossSectionAccessor();

  void OpenCrossSectionFile(const char* fn);

  static CrossSectionAccessor* csa_;
  TFile* xs_file_;

};

#endif // _CROSS_SECTION_ACCESSOR_H
