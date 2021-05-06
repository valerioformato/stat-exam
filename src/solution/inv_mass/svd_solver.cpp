// our headers
#include "svd_solver.h"

// ext headers
#include <fmt/format.h>
#include <fmt/ranges.h>

namespace Solvers {
void SVD::Unfold() { RescaleVariables(); }

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

  for (unsigned int i = 0; i < m_resolution_matrix->GetNcols(); i++) {
    if (measured(i) == 0)
      continue;

    // assume error on measured events is poissonian
    double err = std::sqrt(measured(i));

    measured(i) /= err;

    for (unsigned int j = 0; j < m_resolution_matrix->GetNrows(); j++) {
      resolution_matrix(i, j) /= err;
    }
  }
}

} // namespace Solvers