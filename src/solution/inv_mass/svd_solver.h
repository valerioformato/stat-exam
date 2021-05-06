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
                                                              Utils::HistoToMatrix(resolution_matrix.get())} {};

  void Unfold();

  std::shared_ptr<TH1D> GetUnfolded() { return m_unfolded; }

private:
  std::unique_ptr<TMatrixD> m_resolution_matrix = nullptr;
  std::unique_ptr<TVectorD> m_measured = nullptr;
  std::unique_ptr<TVectorD> m_xini = nullptr;

  std::shared_ptr<TH1D> m_unfolded = nullptr;

  void BuildSVDDecomposition();
  void RescaleVariables();
};
} // namespace Solvers