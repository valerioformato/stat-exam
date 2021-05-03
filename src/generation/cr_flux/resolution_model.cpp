#include "resolution_model.h"

#include <TCanvas.h>
#include <TF1.h>

namespace Utils {
ResolutionModel::ResolutionModel(TList *parList) {
  m_mean = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(0))->Clone("res_mean")));
  m_sigma1 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(1))->Clone("res_sigma1")));
  m_sigma2 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(2))->Clone("res_sigma2")));
  m_frac1 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(3))->Clone("res_frac1")));
  m_frac2 = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(4))->Clone("res_frac2")));
  m_alpha = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(5))->Clone("res_alpha")));
  m_beta = *(static_cast<TH1D *>(static_cast<TH1D *>(parList->At(6))->Clone("res_beta")));

  auto m_fun = std::make_unique<TF1>(
      "model",
      [this](double *x, double *p) {
        return model(x[0], {p[0], p[1], p[2], p[3], p[4], p[5], p[6]});
      },
      -1000, 1000, 7);
  m_fun->SetNpx(100000);
}

double ResolutionModel::operator()(double trueRigidity) { return trueRigidity; }

void ResolutionModel::LoadParameters(double trueRigidity) {
  static unsigned int lastBin = -1;

  auto thisBin = m_mean.FindBin(trueRigidity);
  if (thisBin == lastBin) {
    return;
  }

  m_fun->SetParameter(0, m_mean.GetBinContent(thisBin));
  m_fun->SetParameter(1, m_sigma1.GetBinContent(thisBin));
  m_fun->SetParameter(2, m_sigma2.GetBinContent(thisBin));
  m_fun->SetParameter(3, m_frac1.GetBinContent(thisBin));
  m_fun->SetParameter(4, m_frac2.GetBinContent(thisBin));
  m_fun->SetParameter(5, m_alpha.GetBinContent(thisBin));
  m_fun->SetParameter(6, m_beta.GetBinContent(thisBin));

  lastBin = thisBin;
}

double ResolutionModel::model(double x, std::array<double, 7> p) {
  auto expo = [&x, &p]() { return x >= p[0] ? std::exp(-1 * p[5] * (x - p[0])) : std::exp(p[6] * (x - p[0])); };

  double result = p[3] * TMath::Gaus(x, p[0], p[1], true);
  result += p[4] * TMath::Gaus(x, p[0], p[2], true);
  result += (1 - p[3] - p[4]) * expo();

  return result;
};
} // namespace Utils