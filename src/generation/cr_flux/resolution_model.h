#include <TF1.h>
#include <TH1D.h>
#include <TList.h>

namespace Utils {
class ResolutionModel {
public:
  ResolutionModel(TList *parList);
  ~ResolutionModel() = default;

  double operator()(double trueRigidity);

private:
  void LoadParameters(double trueRigidity);
  double GetRandom();

  std::array<double, 7> m_pars;
  std::unique_ptr<TF1> m_expoFun;

  TH1D m_mean;
  TH1D m_sigma1;
  TH1D m_sigma2;
  TH1D m_frac1;
  TH1D m_frac2;
  TH1D m_alpha;
  TH1D m_beta;
};
} // namespace Utils