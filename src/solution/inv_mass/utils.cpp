// our headers
#include "utils.h"

// c++ headers
#include <numeric>

namespace Utils {
void InitConstants() { std::iota(begin(binning), end(binning), (max_mass - min_mass) / (binning.size() - 1)); }

void DivideByBinWidth(TH1D *histo) {
  for (Int_t i = -1; i < histo->GetNbinsX() + 1; i++)
    histo->SetBinContent(i + 1, histo->GetBinContent(i + 1) / histo->GetBinWidth(i + 1));
}

void MultiplyByBinWidth(TH1D *histo) {
  for (Int_t i = -1; i < histo->GetNbinsX() + 1; i++)
    histo->SetBinContent(i + 1, histo->GetBinContent(i + 1) * histo->GetBinWidth(i + 1));
}

} // namespace Utils