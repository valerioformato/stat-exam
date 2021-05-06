// our headers
#include "svd_solver.h"

// ROOT headers
#include <TCanvas.h>
#include <TDecompSVD.h>

// ext headers
#include <fmt/format.h>
#include <fmt/ranges.h>

namespace Solvers {
void SVD::Unfold(unsigned int cutoff_index) {
  RescaleVariables();

  // get references to improve code readability
  TVectorD &measured = *m_measured;
  TMatrixD &resolution_matrix = *m_resolution_matrix;

  TMatrixD Cinv = C->Invert();
  TMatrixD ACinv = resolution_matrix * Cinv;

  TDecompSVD decomp{ACinv};
  bool status = decomp.Decompose();
  if (!status) {
    fmt::print("{:<10} Failed SVD decomposition\n", "Error:");
  }

  auto vT = decomp.GetV();
  vT.Transpose(decomp.GetV());
  auto uT = decomp.GetU();
  uT.Transpose(decomp.GetU());

  auto CT = *C;
  CT.Transpose(*C);

  auto s = decomp.GetSig();
  auto d = uT * measured;

  PlotSingVal(s, d);

  if (cutoff_index > s.GetNrows()) {
    throw std::runtime_error(fmt::format("The cutoff index is out of range (last index: {})", s.GetNrows() - 1));
  }

  // here we do the regularization part
  TVectorD z_tau = d;
  TVectorD z_tau_err = d;
  double tau = s(cutoff_index) * s(cutoff_index);
  for (unsigned int i = 0; i < z_tau.GetNrows(); i++) {
    z_tau(i) *= s(i) / (s(i) * s(i) + tau);
    z_tau_err(i) = pow(s(i) / (s(i) * s(i) + tau), 2);
  }
  TMatrixD z_tau_corr{z_tau.GetNrows(), z_tau.GetNrows()};
  TMatrixDDiag{z_tau_corr} = z_tau_err;

  auto w = Cinv * (decomp.GetV() * z_tau);
  TMatrixD w_err_cov = Cinv * (decomp.GetV() * (z_tau_corr * (vT * CT)));

  RescaleOutput(w, w_err_cov);
}

void SVD::PlotSingVal(const TVectorD &s, const TVectorD &d) {
  auto testS = std::make_unique<TH1D>("testS", ";k;s", s.GetNrows(), 0, s.GetNrows());
  auto testD = std::make_unique<TH1D>("testD", ";k;d", d.GetNrows(), 0, d.GetNrows());

  for (unsigned int i = 0; i < s.GetNrows(); i++) {
    testS->SetBinContent(i + 1, s(i));
  }
  for (unsigned int i = 0; i < d.GetNrows(); i++) {
    testD->SetBinContent(i + 1, std::fabs(d(i)));
  }

  auto c = std::make_unique<TCanvas>("c", "", 0, 0, 1200, 800);
  c->Divide(2, 1);
  c->cd(1);
  testD->Draw();
  gPad->SetLogy();
  c->cd(2);
  testS->Draw();
  gPad->SetLogy();
  c->Print("SVD_DS_plots.pdf");
}

void SVD::RescaleVariables() {
  // get references to improve code readability
  TVectorD &measured = *m_measured;
  TMatrixD &resolution_matrix = *m_resolution_matrix;

  m_xini = std::make_unique<TVectorD>(m_resolution_matrix->GetNcols());
  TVectorD &xini = *m_xini;
  for (unsigned int j = 0; j < m_xini->GetNrows(); j++) {
    xini(j) = 0;
    for (unsigned int i = 0; i < m_resolution_matrix->GetNcols(); i++) {
      xini(j) += resolution_matrix(i, j);
    }
  }

  m_measured_err = std::make_unique<TVectorD>(m_measured->GetNrows());
  TVectorD &measured_err = *m_measured_err;
  for (unsigned int i = 0; i < m_resolution_matrix->GetNcols(); i++) {
    if (measured(i) == 0)
      continue;

    // assume error on measured events is poissonian
    double err = std::sqrt(measured(i));
    measured_err(i) = err;

    measured(i) /= err;

    for (unsigned int j = 0; j < m_resolution_matrix->GetNrows(); j++) {
      resolution_matrix(i, j) /= err;
    }
  }
}

void SVD::RescaleOutput(const TVectorD &w, TMatrixD &w_err_cov) {
  TVectorD &xini = *m_xini;
  TVectorD &measured_err = *m_measured_err;

  TMatrixDDiag w_diag{w_err_cov};

  for (unsigned int i = 0; i < w.GetNrows(); i++) {
    m_unfolded->SetBinContent(i, w(i) * xini(i));
    m_unfolded->SetBinError(i, std::sqrt(w_diag(i)) * xini(i));
  }

  auto c = std::make_unique<TCanvas>("c", "", 0, 0, 1200, 800);
  m_unfolded->Draw();
  c->Print("debug.pdf");
}

} // namespace Solvers