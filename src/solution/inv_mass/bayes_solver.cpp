// our headers
#include "bayes_solver.h"
#include "utils.h"

// ROOT headers
#include <TCanvas.h>

// external headers
#include <fmt/format.h>
#include <iostream>

namespace Solvers {

void Bayes::Unfold(unsigned int max_iterations, double threshold) {

  if (!CheckBinning()) {
    return;
  }

  // default comparison function for histograms. Slightly adapted from
  // https://stats.stackexchange.com/questions/184101/comparing-two-histograms-using-chi-square-distance
  if (!m_histo_comparator) {
    m_histo_comparator = [](TH1D *h1, TH1D *h2) {
      double chi2 = 0;
      unsigned int count = 0;
      for (unsigned int i = 0; i < h1->GetNbinsX(); i++) {
        if (h1->GetBinError(i + 1) != 0 || h2->GetBinError(i + 1) != 0) {
          count++;
          chi2 += 0.5 * pow(h1->GetBinContent(i + 1) - h2->GetBinContent(i + 1), 2) /
                  (pow(h1->GetBinError(i + 1), 2) + pow(h2->GetBinError(i + 1), 2));
        }
      }

      return chi2 / (count - 1);
    };
  }

  fmt::print("Starting unfolding procedure...\n");

  auto canv = std::make_unique<TCanvas>("c", "", 0, 0, 1200 + 4, 800 + 28);
  canv->Print("debug.pdf[");

  for (unsigned int iter = 0; iter < max_iterations; iter++) {
    TH1D *old_unfolded = nullptr;

    if (iter > 0) {
      // update prior
      for (unsigned int ib = 0; ib < m_prior->GetNbinsX(); ib++) {
        m_prior->SetBinContent(ib + 1, m_unfolded->GetBinContent(ib + 1));
      }

      old_unfolded = static_cast<TH1D *>(m_unfolded->Clone("old_unfolded"));
    }

    // first step: smooth the prior
    SmoothPrior();

    // The prior is a 'probability' so it has to be normalized. Since smoothing could change normalization, we force
    // it here...
    m_prior->Scale(1. / (m_prior->GetSumOfWeights() + m_prior->GetBinContent(0) +
                         m_prior->GetBinContent(m_prior->GetNbinsX() + 1)));

    // build the theta matrix
    BuildBayesMatrix();

    if (!m_unfolded) {
      m_unfolded = std::make_shared<TH1D>(*m_prior);
      m_unfolded->SetName("m_unfolded");
    }
    m_unfolded->Reset();

    for (int i = -1; i < m_unfolded->GetNbinsX() + 1; i++) {
      double num = 0, den = 0, err = 0;

      for (int j = -1; j < m_bayes_matrix->GetNbinsY() + 1; j++) {
        num += m_bayes_matrix->GetBinContent(m_bayes_matrix->GetBin(i + 1, j + 1)) * m_measured->GetBinContent(j + 1);
        den += m_resolution_matrix->GetBinContent(m_resolution_matrix->GetBin(i + 1, j + 1));
        err += pow(m_bayes_matrix->GetBinContent(m_bayes_matrix->GetBin(i + 1, j + 1)) * m_measured->GetBinError(j + 1),
                   2);
      }

      if (den) {
        m_unfolded->SetBinContent(i + 1, num / den);
        m_unfolded->SetBinError(i + 1, sqrt(err) / den);
      }
    }

    { // plotting this iteration status
      m_measured->SetLineColor(kRed);
      m_measured->Draw("E");
      m_unfolded->Draw("E same");
      canv->Print("debug.pdf");
    }

    if (iter > 0) {
      double conv_check = m_histo_comparator(m_unfolded.get(), old_unfolded);
      fmt::print(" -- Iteration {}, convergence value at {}\n", iter, conv_check);
      if (conv_check < threshold) {
        fmt::print("    Successfully converged\n");
        break;
      }
    }

  } // end iteration

  canv->Print("debug.pdf]");
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
    m_bayes_matrix = std::make_unique<TH2D>(*m_resolution_matrix);
    m_bayes_matrix->SetName("m_bayes_matrix");
  }
  m_bayes_matrix->Reset();

  // precompute denominators in bayes theorem
  std::vector<double> prob_norms(m_bayes_matrix->GetNbinsY() + 2, 0);
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
      // fmt::print("rm({}, {}) = {}\n", ibinx + 1, ibiny + 1, m_resolution_matrix->GetBinContent(globalbin));
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
      fmt::print("{:<10} Resolution matrix is not normalized.\n", "Warning:"); 
      fmt::print("{:<10} Column {} has sum = {}\n", "", i + 1, sum);
      return false;
    }
  }

  return true;
}

void Bayes::InitFlatPrior() {
  m_prior = std::make_unique<TH1D>("m_prior", "", m_resolution_matrix->GetXaxis()->GetXbins()->GetSize() - 1,
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