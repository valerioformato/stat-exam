// our headers
#include "spectrum.h"

namespace Utils {
constexpr float amu = 0.93146;

constexpr double GetMass(double A, double Z) {
  double mass = A * amu;
  if (Z == 1) {
    if (A == 1)
      mass = 0.9382723;
    if (A == 2)
      mass = 1.875613;
  } else if (Z == 2) {
    if (A == 3)
      mass = 2.808;
    if (A == 4)
      mass = 3.7273799;
  } else if (Z == 5) {
    if (A == 10)
      mass = 9.3245;
    if (A == 11)
      mass = 10.2526;
  } else if (Z == 6) {
    if (A == 12)
      mass = 11.175;
  }

  return mass;
}

constexpr double E2R(float EA, double A, double Z) {

  /*------------------------------------------------*\
  |     Conversion from KinEnergy/N to Rigidity      |
  \*------------------------------------------------*/

  double R =
      ((double)1 / Z) * sqrt(pow(A * EA, 2) + 2 * GetMass(A, Z) * A * EA);

  return R;
}

constexpr double R2E(float R, double A, double Z) {

  /*------------------------------------------------*\
  |     Conversion from Rigidity to KinEnergy/N      |
  \*------------------------------------------------*/

  double mass = GetMass(A, Z);
  double EA = ((double)1 / A) * (sqrt(pow(Z * R, 2) + pow(mass, 2)) - mass);

  return EA;
}

constexpr double R2B(float R, double A, double Z) {

  /*------------------------------------------------*\
  |         Conversion from Rigidity to beta         |
  \*------------------------------------------------*/

  double B = pow(1 + pow(GetMass(A, Z) / (Z * R), 2), -0.5);

  return B;
}

constexpr double B2R(float B, double A, double Z) {

  /*------------------------------------------------*\
  |         Conversion from beta to Rigidity         |
  \*------------------------------------------------*/

  double R = GetMass(A, Z) * B / sqrt(1 - B * B) / Z;

  return R;
}

constexpr double B2E(float B, double A, double Z) {

  /*------------------------------------------------*\
  |          Conversion from beta to Energy          |
  \*------------------------------------------------*/

  return (double)GetMass(A, Z) / A * (1. / sqrt(1 - B * B) - 1);
}

constexpr double E2B(float EA, double A, double Z) {

  /*------------------------------------------------*\
  |          Conversion from Energy to beta          |
  \*------------------------------------------------*/

  double bsq = 1 - 1. / pow(A / GetMass(A, Z) * EA + 1, 2);
  return sqrt(bsq);
}

std::unique_ptr<TF1> GetSpectrum(unsigned int Z, unsigned int A) {
  auto modelFun = [A, Z](double *_x, double *_p) -> double {
    const double R = _x[0];

    const double Mass = GetMass(A, Z);
    const double x = R2E(R, A, Z);
    const double enm = x * Mass + Mass;
    const double en = x * Mass + Mass + Z * _p[7];
    // const double xb = R2E(_p[5], A, Z);

    // by-eye parametrization of primary fluxes
    // exp = p0 - gamma*log(R) - p1/(p2 R)^0.6
    double exponent = 0.0 - _p[1] * TMath::Log(R) - _p[2] / pow(R * _p[3], 0.6);

    // add spectral index break as for AMS papers
    exponent += _p[4] * TMath::Log(1 + TMath::Power(R / _p[5], _p[6] / _p[4]));

    // add normalization and get flux
    double flux = _p[0] * TMath::Exp(exponent);

    // E -> R jacobian
    flux /= pow((float)A / Z, 2) * ((x + 0.938) / R);

    // apply simple force-field model
    flux *= (enm * enm - Mass * Mass) / (en * en - Mass * Mass);

    return flux < 0 ? 0 : flux;
  };

  auto fun = std::make_unique<TF1>("spectrum", modelFun, 1, 6000, 8);
  fun->SetParameter(0, 1);     // norm
  fun->SetParameter(1, 2.7);   // gamma
  fun->SetParameter(2, 1.5);   // p1
  fun->SetParameter(3, 1);     // p2
  fun->SetParameter(4, 0.027); // s
  fun->SetParameter(5, 230);   // R break
  fun->SetParameter(6, 0.1);   // delta gamma
  fun->SetParameter(7, 0.4);   // phi
  fun->SetNpx(100000);
  return std::move(fun);
}

} // namespace Utils