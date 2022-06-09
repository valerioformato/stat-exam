// our headers
#include "utils.h"

// ROOT headers
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>

// c++ headers
#include <functional>
#include <memory>

namespace Solvers {
class ForwardFolding {
public:
  ForwardFolding(std::shared_ptr<TH1D> measured, std::shared_ptr<TH2D> resolution_matrix)
      : m_measured{Utils::HistoToVector(measured.get())}, m_resolution_matrix{
                                                              Utils::HistoToMatrix(resolution_matrix.get())} {
    m_unfolded = std::make_shared<TH1D>("m_unfolded_forward", "", resolution_matrix->GetNbinsX(),
                                        resolution_matrix->GetXaxis()->GetXbins()->GetArray());

    m_range.first = m_unfolded->GetBinLowEdge(1);
    m_range.second = m_unfolded->GetBinLowEdge(m_unfolded->GetNbinsX() + 1);
  }

  void Unfold();

  std::shared_ptr<TH1D> GetUnfolded() { return m_unfolded; }

  void SetTF1Ansatz(const std::string &formula, const std::vector<float> &initial_params);

  bool IsResMatrixNormalized();

  std::vector<double> GetResult() {return m_result; };

private:
  std::unique_ptr<TMatrixD> m_resolution_matrix = nullptr;
  std::unique_ptr<TVectorD> m_measured = nullptr;

  std::shared_ptr<TH1D> m_unfolded = nullptr;

  std::pair<float, float> m_range;
  std::unique_ptr<TF1> m_model = nullptr;

  std::vector<double> m_result;

  void NormalizeResMatrix();

};
} // namespace Solvers