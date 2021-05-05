// our headers
#include "bayes_solver.h"
#include "utils.h"

// external headers
#include <fmt/format.h>

namespace Solvers {

void Bayes::Unfold() {

  if (!CheckBinning()) {
    return;
  }

  // first step: smooth the prior
  SmoothPrior();

  // The prior is a 'probability' so it has to be normalized. Since smoothing could change normalization, we force
  // it here...
  m_prior->Scale(
      1. / (m_prior->GetSumOfWeights() + m_prior->GetBinContent(0) + m_prior->GetBinContent(m_prior->GetNbinsX() + 1)));

  fmt::print("Starting unfolding procedure...\n");

  // build the theta matrix
  BuildBayesMatrix();

  if (!m_unfolded) {
    m_unfolded = std::make_shared<TH1D>(*m_prior);
    m_unfolded->SetName("m_unfolded");
  }
  m_unfolded->Reset();

  for (int i = -1; i < m_unfolded->GetNbinsX() + 1; i++) {
    double num = 0, den = 0, err = 0;

    // cout << "=================================================================" << endl;
    // cout << "X bin " << i << "  R = " << m_unfolded->GetBinCenter(i+1) << endl;

    for (int j = -1; j < m_bayes_matrix->GetNbinsY() + 1; j++) {
      num += m_bayes_matrix->GetBinContent(m_bayes_matrix->GetBin(i + 1, j + 1)) * m_measured->GetBinContent(j + 1);
      den += m_resolution_matrix->GetBinContent(m_resolution_matrix->GetBin(i + 1, j + 1));
      // err += pow( m_bayes_matrix->GetBinContent(m_bayes_matrix->GetBin(i+1,j+1))*m_measured->GetBinError(j+1), 1 );
      err +=
          pow(m_bayes_matrix->GetBinContent(m_bayes_matrix->GetBin(i + 1, j + 1)) * m_measured->GetBinError(j + 1), 2);

      // printf("%03i ", j);
      // cout << m_bayes_matrix->GetBinContent( m_bayes_matrix->GetBin(i+1,j+1) ) << " * " <<
      // m_measured->GetBinContent(j+1) << " +     ("
      // << m_bayes_matrix->GetBinContent( m_bayes_matrix->GetBin(i+1,j+1) )*m_measured->GetBinContent(j+1) << ")" <<
      // endl;
    }

    // cout << " = " << num << " +- " << sqrt(err) << "  -> / " << den << " = " << (den>0 ? num/den : 0) << endl;
    // cout << "=================================================================" << endl;
    // cout << endl;

    if (den) {
      m_unfolded->SetBinContent(i + 1, num / den);
      m_unfolded->SetBinError(i + 1, sqrt(err) / den);
    }
  }
}

void Bayes::BuildBayesMatrix() {
  if (!m_resolution_matrix) {
    throw std::runtime_error("Resolution matrix is not set!");
  }

  if (!m_prior) {
    throw std::runtime_error("Prior is not set!");
  }

  NormalizeResMatrix();

  if (!m_bayes_matrix) {
    m_bayes_matrix = std::make_shared<TH2D>(*m_resolution_matrix);
    m_bayes_matrix->SetName("m_bayes_matrix");
  }
  m_bayes_matrix->Reset();

  // precompute denominators in bayes theorem
  std::vector<double> prob_norms(m_bayes_matrix->GetNbinsY() + 1, 0);
  for (int j = -1; j < m_bayes_matrix->GetNbinsY() + 1; j++) {
    for (int i = -1; i < m_bayes_matrix->GetNbinsX() + 1; i++) {
      prob_norms[j + 1] +=
          m_resolution_matrix->GetBinContent(m_resolution_matrix->GetBin(i + 1, j + 1)) * m_prior->GetBinContent(i + 1);
    }
  }

  // run bayes theorem for each matrix entry
  for (int i = -1; i < m_bayes_matrix->GetNbinsX() + 1; i++) {
    for (int j = -1; j < m_bayes_matrix->GetNbinsY() + 1; j++) {

      double s =
          m_resolution_matrix->GetBinContent(m_resolution_matrix->GetBin(i + 1, j + 1)) * m_prior->GetBinContent(i + 1);
      // double ls = 0;
      // for (int k = -1; k < m_resolution_matrix->GetNbinsX() + 1; k++)
      //   ls += m_resolution_matrix->GetBinContent(m_resolution_matrix->GetBin(k + 1, j + 1)) *
      //         m_prior->GetBinContent(k + 1);

      if (prob_norms[j + 1])
        m_bayes_matrix->SetBinContent(m_bayes_matrix->GetBin(i + 1, j + 1), s / prob_norms[j + 1]);
    }
  }
}

void Bayes::NormalizeResMatrix() {
  if (IsResMatrixNormalized())
    return;

  fmt::print("Normalizing resolution matrix...");

  for (int ibinx = -1; ibinx < m_resolution_matrix->GetNbinsX() + 1; ibinx++) {
    double norm = 0;

    for (int ibiny = -1; ibiny < m_resolution_matrix->GetNbinsY() + 1; ibiny++) {
      int globalbin = m_resolution_matrix->GetBin(ibinx + 1, ibiny + 1);
      norm += m_resolution_matrix->GetBinContent(globalbin);
    }

    if (!norm)
      continue;

    for (int ibiny = -1; ibiny < m_resolution_matrix->GetNbinsY() + 1; ibiny++) {
      int globalbin = m_resolution_matrix->GetBin(ibinx + 1, ibiny + 1);
      m_resolution_matrix->SetBinContent(globalbin, m_resolution_matrix->GetBinContent(globalbin) / norm);
    }
  }

  fmt::print(" done\n");
}

bool Bayes::IsResMatrixNormalized() {
  // check if any of the columns has sum greater than 1

  double sum = 0;
  for (int i = -1; i < m_resolution_matrix->GetNbinsX() + 1; i++) {
    sum = 0;
    for (int j = -1; j < m_resolution_matrix->GetNbinsY() + 1; j++)
      sum += m_resolution_matrix->GetBinContent(m_resolution_matrix->GetBin(i + 1, j + 1));
    if (sum > 1.0000001) {
      fmt::print("Warning: Resolution matrix is not normalized.\nColumn {} has sum = {}\n", i + 1, sum);
      return false;
    }
  }

  return true;
}

void Bayes::InitFlatPrior() {
  m_prior = std::make_shared<TH1D>("m_prior", "", m_resolution_matrix->GetXaxis()->GetXbins()->GetSize() - 1,
                                   m_resolution_matrix->GetXaxis()->GetXbins()->GetArray());

  for (int ib = 0; ib < m_prior->GetNbinsX(); ib++) {
    m_prior->SetBinContent(ib + 1, 1);
  }
}

void Bayes::SmoothPrior() {
  if (!m_prior) {
    throw std::runtime_error("Prior is not set!");
  }

  if (!m_smoother) {
    m_smoother = [](TH1D *h) { h->Smooth(); };
  }

  // when dealing with non-uniform bins you should always smooth the
  // underlying p.d.f., not the histogram itself
  Utils::DivideByBinWidth(m_prior.get());

  m_smoother(m_prior.get());

  // restore the "prior" histogram
  Utils::MultiplyByBinWidth(m_prior.get());
}

bool Bayes::CheckBinning() {
  // check for binning ocnsistency between the various histograms.

  int n_meas = m_measured->GetNbinsX();
  int n_smy = m_resolution_matrix->GetNbinsY();

  if (n_meas != n_smy) {
    fmt::print(" --- Different binning between measured and smearing\n");
    fmt::print(" ---- measured: {} bins;  smearing: {} bins on Y-axis\n", n_meas, n_smy);
    return false;
  }

  const double *x_meas = m_measured->GetXaxis()->GetXbins()->GetArray();
  const double *y_sm = m_resolution_matrix->GetYaxis()->GetXbins()->GetArray();

  for (int ib = 0; ib < n_meas + 1; ib++) {
    if (fabs(x_meas[ib] - y_sm[ib]) / x_meas[ib] > 1e-4) {
      fmt::print(" --- Different binning between measured and smearing\n");
      fmt::print(" ---- x_meas[{0}] = {1}; y_sm[{0}] = {2};", ib, x_meas[ib], y_sm[ib]);
      return false;
    } else
      continue;
  }

  if (m_prior) {
    int n_prior = m_prior->GetNbinsX();
    int n_smx = m_resolution_matrix->GetNbinsX();

    if (n_prior != n_smx) {
      fmt::print(" --- Different binning between prior and smearing\n");
      return false;
    }

    const double *x_prior = m_prior->GetXaxis()->GetXbins()->GetArray();
    const double *x_sm = m_resolution_matrix->GetXaxis()->GetXbins()->GetArray();

    for (int ib = 0; ib < n_prior + 1; ib++) {
      if (fabs(x_prior[ib] - x_sm[ib]) / x_sm[ib] > 1e-4) {
        fmt::print(" --- Different binning between prior and smearing\n");
        return false;
      } else
        continue;
    }
  }

  return true;
}

} // namespace Solvers