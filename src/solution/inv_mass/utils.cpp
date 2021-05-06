// our headers
#include "utils.h"

// c++ headers
#include <iostream>
#include <numeric>

namespace Utils {

void DivideByBinWidth(TH1D *histo) {
  for (Int_t i = -1; i < histo->GetNbinsX() + 1; i++)
    histo->SetBinContent(i + 1, histo->GetBinContent(i + 1) / histo->GetBinWidth(i + 1));
}

void MultiplyByBinWidth(TH1D *histo) {
  for (Int_t i = -1; i < histo->GetNbinsX() + 1; i++)
    histo->SetBinContent(i + 1, histo->GetBinContent(i + 1) * histo->GetBinWidth(i + 1));
}

TVectorD *HistoToVector(TH1D *histo) {
  TVectorD *vector = new TVectorD(histo->GetNbinsX() + 2);

  for (unsigned int ix = 0; ix < histo->GetNbinsX() + 2; ix++) {
    vector->operator()(ix) = histo->GetBinContent(ix);
  }

  return vector;
}

TMatrixD *HistoToMatrix(TH2D *histo) {
  TMatrixD *matrix = new TMatrixD(histo->GetNbinsX() + 2, histo->GetNbinsY() + 2);

  for (unsigned int ix = 0; ix < histo->GetNbinsX() + 2; ix++) {
    for (unsigned int iy = 0; iy < histo->GetNbinsY() + 2; iy++) {
      auto globalbin = histo->GetBin(ix, iy);
      matrix->operator()(ix, iy) = histo->GetBinContent(globalbin);
    }
  }

  return matrix;
}

} // namespace Utils