// ROOT headers
#include <TH1D.h>
#include <TH2D.h>

// c++ headers
#include <functional>
#include <memory>

namespace Solvers {
class Bayes {
public:
  Bayes(std::shared_ptr<TH1D> measured, std::shared_ptr<TH2D> resolution_matrix)
      : m_measured{measured}, m_resolution_matrix{resolution_matrix} {};

  bool IsResMatrixNormalized();
  void Unfold();

  void SetCustomSmoother(std::function<void(TH1D *)> smoother) { m_smoother = smoother; }
  void SetPrior(TH1D *prior) { m_prior = std::make_shared<TH1D>(*prior); }
  void InitFlatPrior();

  std::shared_ptr<TH1D> GetUnfolded() { return m_unfolded; }

private:
  std::shared_ptr<TH2D> m_resolution_matrix = nullptr;
  std::shared_ptr<TH1D> m_measured = nullptr;
  std::shared_ptr<TH1D> m_unfolded = nullptr;
  std::shared_ptr<TH1D> m_prior = nullptr;

  std::shared_ptr<TH2D> m_bayes_matrix = nullptr;

  std::function<void(TH1D *)> m_smoother;

  bool CheckBinning();

  void BuildBayesMatrix();
  void NormalizeResMatrix();
  void SmoothPrior();
};
} // namespace Solvers