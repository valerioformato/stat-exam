// our headers
#include "forward_solver.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

#include <fmt/format.h>

namespace Solvers {
void ForwardFolding::Unfold() {
  NormalizeResMatrix();

  const unsigned int total_events = m_measured->Sum();

  // prepare likelihood FCN
  auto logLikelihood = [this, total_events](double const *p) {
    double result{0};

    m_model->SetParameters(p);
    double normalization = m_model->Integral(m_range.first, m_range.second);

    int start_bin = m_unfolded->FindBin(0.5);
    int end_bin = m_unfolded->FindBin(9.5);

    TVectorD real_counts(m_measured->GetNoElements());
    for (int ibin = start_bin; ibin < end_bin + 1; ++ibin) {
      const double m_start = m_unfolded->GetBinLowEdge(ibin + 1);
      const double m_end = m_unfolded->GetBinLowEdge(ibin + 2);

      real_counts[ibin] = std::max(0.0, normalization * m_model->Integral(m_start, m_end)) * total_events;
    }

    TVectorD expected = (*m_resolution_matrix) * real_counts;

    for (int ibin = start_bin; ibin < end_bin + 1; ++ibin) {
      const double meas_counts = (*m_measured)[ibin];

      // fixed normalization
      double expected_counts = expected[ibin];

      // free normalization
      // const double expected_counts = std::max(0.0, m_model->Integral(m_start, m_end));

      double partial_ll = meas_counts * std::log(expected_counts) - expected_counts;
      if (meas_counts > 10) {
        partial_ll -= meas_counts * std::log(meas_counts) - meas_counts + 0.5 * std::log(2 * TMath::Pi() * meas_counts);
      } else {
        partial_ll -= std::log(TMath::Factorial(meas_counts));
      }

      // fmt::print("{}: {} vs {}  ({})\n", ibin, expected_counts, meas_counts, partial_ll);

      result += partial_ll;
    }

    return -1.0 * result;
  };

  // build minimizer
  ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimizer->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
  minimizer->SetMaxIterations(10000);        // for GSL
  minimizer->SetTolerance(1e-14);
  minimizer->SetPrintLevel(1);

  ROOT::Math::Functor fun(logLikelihood, m_model->GetNpar());
  minimizer->SetFunction(fun);

  constexpr double initial_step = 0.01;

  // setup minimizer variables
  for (int ipar = 0; ipar < m_model->GetNpar(); ++ipar) {
    minimizer->SetLowerLimitedVariable(ipar, m_model->GetParName(ipar), m_model->GetParameter(ipar), initial_step, 0.0);
  }

  minimizer->Minimize();

  const double *min_val = minimizer->X();
  for (int ipar = 0; ipar < m_model->GetNpar(); ++ipar) {
    fmt::print("{} = {}\n", m_model->GetParName(ipar), min_val[ipar]);
    m_result.push_back(min_val[ipar]);
  }
}

void ForwardFolding::SetTF1Ansatz(const std::string &formula, const std::vector<float> &initial_params) {
  m_model = std::make_unique<TF1>("forward_model", formula.c_str(), m_range.first, m_range.second);
  for (unsigned int ipar = 0; ipar < initial_params.size(); ++ipar) {
    m_model->SetParameter(ipar, initial_params[ipar]);
  }

  return;

  TF1 temp{"temp", formula.c_str(), m_range.first, m_range.second};
  auto n_normpar = temp.GetNpar();

  auto norm_formula = fmt::format("[{}] * ({})", n_normpar, formula);
  m_model = std::make_unique<TF1>("forward_model", norm_formula.c_str(), m_range.first, m_range.second);
  for (unsigned int ipar = 0; ipar < initial_params.size(); ++ipar) {
    m_model->SetParameter(ipar, initial_params[ipar]);
  }

  // temp value
  m_model->SetParameter(n_normpar, 1);

  auto normalization = m_model->Integral(m_range.first, m_range.second);
  if (std::fabs(normalization - 1) > 1e-6) {
    m_model->SetParameter(n_normpar, 1.0f / normalization);
  }

  fmt::print("{}, norm = {}\n", norm_formula, normalization);
}

bool ForwardFolding::IsResMatrixNormalized() {
  // check if any of the columns has sum greater than 1

  double sum = 0;
  for (int i = 0; i < m_resolution_matrix->GetNcols(); i++) {
    sum = 0;
    for (int j = 0; j < m_resolution_matrix->GetNrows(); j++)
      sum += (*m_resolution_matrix)(i, j);
    if (sum > 1.0000001) {
      fmt::print("{:<10} Resolution matrix is not normalized.\n", "Warning:");
      fmt::print("{:<10} Column {} has sum = {}\n", "", i + 1, sum);
      return false;
    }
  }

  return true;
}

void ForwardFolding::NormalizeResMatrix() {
  if (IsResMatrixNormalized())
    return;

  fmt::print("Normalizing resolution matrix...");

  for (int i = 0; i < m_resolution_matrix->GetNcols(); i++) {
    double norm = 0;

    for (int j = 0; j < m_resolution_matrix->GetNrows(); j++) {
      norm += (*m_resolution_matrix)(i, j);
    }

    if (!norm)
      continue;

    for (int j = 0; j < m_resolution_matrix->GetNrows(); j++) {
      (*m_resolution_matrix)(i, j) = (*m_resolution_matrix)(i, j) / norm;
    }
  }

  fmt::print(" done\n");
}
} // namespace Solvers