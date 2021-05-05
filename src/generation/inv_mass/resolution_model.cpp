#include "resolution_model.h"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>

namespace Utils {
double ResolutionModel::operator()(double trueMomentum) {
  auto dp_p = [this](double mom) {
    constexpr double thresh = 180;
    if (mom < thresh)
      return 0.01 * (m_pars[2] * mom * mom + m_pars[1] * mom + m_pars[0]);
    else
      //              f'(thresh)                          * (x - thresh)   + f(thresh)
      return 0.01 * ((1 * m_pars[2] * thresh + m_pars[1]) * (mom - thresh) + m_pars[2] * thresh * thresh +
                     m_pars[1] * thresh + m_pars[0]);
  };

  // fmt::print("sigma = {}\n", dp_p(trueMomentum));

  return trueMomentum + gRandom->Gaus(0, dp_p(trueMomentum) * trueMomentum);
}

} // namespace Utils