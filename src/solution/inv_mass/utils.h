#ifndef utils_h_
#define utils_h_

// ROOT headers
#include <TH1D.h>

// c++ headers
#include <array>

namespace Utils {
static constexpr double min_mass = 0.0, max_mass = 10.0;
static std::array<double, 101> binning{};

static void InitConstants() {
  double delta = (max_mass - min_mass) / (binning.size() - 1);
  std::generate(begin(binning), end(binning), [n = -delta, delta]() mutable { return n += delta; });
}

void DivideByBinWidth(TH1D *histo);
void MultiplyByBinWidth(TH1D *histo);
} // namespace Utils
#endif