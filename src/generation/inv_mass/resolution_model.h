#include <TF1.h>
#include <TH1D.h>
#include <TList.h>

namespace Utils {
class ResolutionModel {
public:
  double operator()(double trueMomentum);

private:
  std::array<double, 3> m_pars{0.5, 0.004753, -0.00001133}; //  interpolated from https://arxiv.org/pdf/1412.6352.pdf
};
} // namespace Utils