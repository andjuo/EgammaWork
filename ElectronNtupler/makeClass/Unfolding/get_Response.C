#ifdef unfold_RMatrix_cxx
#undef unfold_RMatrix_cxx
#endif
#include "unfold_RMatrix.h"
#include <iostream>

// ------------------------------------------------------

inline void HERE(const char *msg) { std::cout << msg << std::endl; }

// ------------------------------------------------------

void get_Response() {
  TString fname="test/DY_check.root";
  TFile fin(fname);
  if (!fin.IsOpen()) { std::cout << "failed to open the file\n"; return; }
  TTree *t=(TTree*)fin.Get("ntupler/ElectronTree");
  if (!t) { std::cout << "failed to get the tree\n"; return; }
  std::cout << "Got tree" <<  std::endl;

  unfold_RMatrix unf(t);
  HERE("unf created");
  unf.Loop();
  HERE("loop done");
  fin.Close();
}

// ------------------------------------------------------
