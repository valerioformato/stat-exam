#include "resolution_model.h"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>

namespace Utils {
double ResolutionModel::operator()(double trueMomentum) {
  auto dp_p = [this](double mom) { return 0.01 * (m_pars[2] * mom * mom + m_pars[1] * mom + m_pars[0]); };

  double sigma = trueMomentum > 250 ? dp_p(250) : dp_p(trueMomentum);

  // fmt::print("sigma = {}\n", sigma);

  return trueMomentum + gRandom->Gaus(0, sigma);
}

} // namespace Utils