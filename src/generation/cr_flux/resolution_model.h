#include <TH1D.h>
#include <TF1.h>
#include <TList.h>

namespace Utils {
class ResolutionModel {
public:
  ResolutionModel(TList *parList);
  ~ResolutionModel() = default;

  double operator()(double trueRigidity);

private:
  void LoadParameters(double trueRigidity);

  double model(double x, std::array<double, 7> p);
  std::unique_ptr<TF1> m_fun;

  TH1D m_mean;
  TH1D m_sigma1;
  TH1D m_sigma2;
  TH1D m_frac1;
  TH1D m_frac2;
  TH1D m_alpha;
  TH1D m_beta;
};
} // namespace Utils