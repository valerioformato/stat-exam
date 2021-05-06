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
      : m_measured{std::make_unique<TH1D>(*measured)}, m_resolution_matrix{
                                                           std::make_unique<TH2D>(*resolution_matrix)} {};

  bool IsResMatrixNormalized();
  void Unfold(unsigned int max_iterations = 5, double threshold = 1);

  void SetCustomSmoother(std::function<void(TH1D *)> smoother) { m_smoother = smoother; }
  void SetCustomConvergenceCheck(std::function<double(TH1D *, TH1D *)> comparator) { m_histo_comparator = comparator; }
  void SetPrior(TH1D *prior) { m_prior = std::make_unique<TH1D>(*prior); }
  void InitFlatPrior();

  std::shared_ptr<TH1D> GetUnfolded() { return m_unfolded; }

private:
  std::unique_ptr<TH2D> m_resolution_matrix = nullptr;
  std::unique_ptr<TH1D> m_measured = nullptr;
  std::shared_ptr<TH1D> m_unfolded = nullptr;
  std::unique_ptr<TH1D> m_prior = nullptr;

  std::unique_ptr<TH2D> m_bayes_matrix = nullptr;

  std::function<void(TH1D *)> m_smoother;
  std::function<double(TH1D *, TH1D *)> m_histo_comparator;

  bool CheckBinning();

  void BuildBayesMatrix();
  void NormalizeResMatrix();
  void SmoothPrior();
};
} // namespace Solvers