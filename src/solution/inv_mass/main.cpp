// our headers
#include "bayes_solver.h"
#include "localization.h"
#include "utils.h"

// ROOT headers
#include <TFile.h>
#include <TTree.h>

// external headers
#include <fmt/format.h>

// c++ headers
#include <memory>

int main(int argc, char const *argv[]) {

  Utils::InitConstants();

  auto inFile_data = std::make_unique<TFile>(
      fmt::format("{}/datafiles/solution/inv_mass/im_data.root", localization::PROJECT_SOURCE_DIR).c_str(), "open");
  TTree *tree_data = inFile_data->Get<TTree>("Data");

  auto inFile_mc = std::make_unique<TFile>(
      fmt::format("{}/datafiles/solution/inv_mass/im_mc.root", localization::PROJECT_SOURCE_DIR).c_str(), "open");
  TTree *tree_mc = inFile_mc->Get<TTree>("Data");

  auto data = std::make_shared<TH1D>("data", "m (GeV);Counts", Utils::binning.size() - 1, Utils::binning.data());

  auto smearing_matrix =
      std::make_shared<TH2D>("smearing_matrix", "m_{gen} (GeV); m_{meas} (GeV)", Utils::binning.size() - 1,
                             Utils::binning.data(), Utils::binning.size() - 1, Utils::binning.data());

  // Bayes
  Solvers::Bayes bayes_solver{data, smearing_matrix};
  bayes_solver.InitFlatPrior();
  bayes_solver.Unfold();

  return 0;
}
