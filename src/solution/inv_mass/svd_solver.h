// our headers
#include "utils.h"

// ROOT headers
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>

// c++ headers
#include <functional>
#include <memory>

namespace Solvers {
class SVD {
public:
  // in this case we want the resolution matrix to be NON NORMALIZED
  // i.e.: each column contains all the generated events (including under/overflow)
  SVD(std::shared_ptr<TH1D> measured, std::shared_ptr<TH2D> resolution_matrix)
      : m_measured{Utils::HistoToVector(measured.get())}, m_resolution_matrix{
                                                              Utils::HistoToMatrix(resolution_matrix.get())} {
    m_unfolded = std::make_shared<TH1D>("m_unfolded_svd", "", resolution_matrix->GetNbinsX(),
                                        resolution_matrix->GetXaxis()->GetXbins()->GetArray());

    C = std::make_unique<TMatrixD>(m_measured->GetNrows(), m_measured->GetNrows());
    for (unsigned int i = 0; i < m_measured->GetNrows(); i++) {
      // small regularization to allow inversion. See https://arxiv.org/pdf/hep-ph/9509307.pdf
      double xi = 1e-4;

      if (i == 0) {
        (*C)(i, i) = -1 + 1e-4;
        (*C)(i + 1, i) = 1;
      } else if (i == m_measured->GetNrows() - 1) {
        (*C)(i, i) = -1 + 1e-4;
        (*C)(i - 1, i) = 1;
      } else {
        (*C)(i, i) = -2 + xi;
        (*C)(i - 1, i) = 1;
        (*C)(i + 1, i) = 1;
      }
    }
  };

  void Unfold(unsigned int cutoff_index);

  std::shared_ptr<TH1D> GetUnfolded() { return m_unfolded; }

private:
  std::unique_ptr<TMatrixD> C = nullptr;

  std::unique_ptr<TMatrixD> m_resolution_matrix = nullptr;
  std::unique_ptr<TVectorD> m_measured = nullptr;
  std::unique_ptr<TVectorD> m_measured_err = nullptr;
  std::unique_ptr<TVectorD> m_xini = nullptr;

  std::shared_ptr<TH1D> m_unfolded = nullptr;

  void BuildSVDDecomposition();
  void RescaleVariables();
  void RescaleOutput(const TVectorD &w, TMatrixD &w_err_cov);

  void PlotSingVal(const TVectorD &s, const TVectorD &d);
};
} // namespace Solvers