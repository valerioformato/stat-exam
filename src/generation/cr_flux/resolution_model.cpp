#include "resolution_model.h"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>

namespace Utils {
ResolutionModel::ResolutionModel(TList *parList) {
  m_mean = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(0))->Clone("res_mean")));
  m_sigma1 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(1))->Clone("res_sigma1")));
  m_sigma2 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(2))->Clone("res_sigma2")));
  m_frac1 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(3))->Clone("res_frac1")));
  m_frac2 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(4))->Clone("res_frac2")));
  m_alpha = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(5))->Clone("res_alpha")));
  m_beta = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(6))->Clone("res_beta")));

  m_expoFun = std::make_unique<TF1>(
      "expoFun",
      [this](double *x, double *p) { return x[0] >= 0 ? std::exp(-1.0 * p[0] * x[0]) : std::exp(p[1] * x[0]); }, -1000,
      1000, 2);
}

double ResolutionModel::operator()(double trueRigidity) {
  LoadParameters(trueRigidity);
  // fmt::print("Parameters loaded...\n");
  return 1.0 / (1.0 / trueRigidity + GetRandom() / 1000.0);
}

void ResolutionModel::LoadParameters(double trueRigidity) {
  static unsigned int lastBin = -1;

  auto thisBin = m_mean.FindBin(trueRigidity);
  if (thisBin == lastBin) {
    return;
  }

  // fmt::print("Setting parameters\n");
  m_pars[0] = m_mean.GetBinContent(thisBin);
  m_pars[1] = m_sigma1.GetBinContent(thisBin);
  m_pars[2] = m_sigma2.GetBinContent(thisBin);
  m_pars[3] = m_frac1.GetBinContent(thisBin);
  m_pars[4] = m_frac2.GetBinContent(thisBin);
  m_pars[5] = m_alpha.GetBinContent(thisBin);
  m_pars[6] = m_beta.GetBinContent(thisBin);
  m_expoFun->SetParameter(0, m_pars[5]);
  m_expoFun->SetParameter(1, m_pars[6]);
  // fmt::print("Parameters set\n");

  lastBin = thisBin;
}

double ResolutionModel::GetRandom() {
  double frac = gRandom->Uniform();

  double result = m_pars[0];

  if (frac < m_pars[3]) {
    // first gaussian
    // fmt::print("Got first gaussian ({} < {})\n", frac, m_pars[3]);
    result += gRandom->Gaus(0.0, m_pars[1]);
  } else if (frac < m_pars[3] + m_pars[4]) {
    // second gaussian
    // fmt::print("Got second gaussian ({} < {} + {})\n", frac, m_pars[3], m_pars[4]);
    result += gRandom->Gaus(0.0, m_pars[2]);
  } else {
    // fmt::print("Got double exponential ({} > {} + {})\n", frac, m_pars[3], m_pars[4]);
    result += m_expoFun->GetRandom();
  }

  return result;
}

} // namespace Utils