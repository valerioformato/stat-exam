// ROOT headers
#include <TH1D.h>

// c++ headers
#include <array>

namespace Utils {
static constexpr double min_mass = 0.0, max_mass = 10.0;
static std::array<double, 101> binning{};

void InitConstants();

void DivideByBinWidth(TH1D *histo);
void MultiplyByBinWidth(TH1D *histo);
} // namespace Utils