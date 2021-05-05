// our headers
#include "bayes_solver.h"
#include "localization.h"
#include "utils.h"

// ROOT headers
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

// external headers
#include <fmt/format.h>
#include <fmt/ranges.h>

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

  // 4-vectors are stored with the convention
  // p = (px, py, pz, E)
  float muon1[4], muon2[4];
  tree_data->SetBranchAddress("muon1", &muon1);
  tree_data->SetBranchAddress("muon2", &muon2);
  for (unsigned long long iEv = 0; iEv < tree_data->GetEntries(); iEv++) {
    tree_data->GetEntry(iEv);
    TLorentzVector p1{muon1};
    TLorentzVector p2{muon2};

    double mass = (p1 + p2).Mag();
    data->Fill(mass);
  }

  float gen_muon1[4], gen_muon2[4];
  tree_mc->SetBranchAddress("meas_muon1", &muon1);
  tree_mc->SetBranchAddress("meas_muon2", &muon2);
  tree_mc->SetBranchAddress("gen_muon1", &gen_muon1);
  tree_mc->SetBranchAddress("gen_muon2", &gen_muon2);
  for (unsigned long long iEv = 0; iEv < tree_mc->GetEntries(); iEv++) {
    tree_mc->GetEntry(iEv);
    TLorentzVector p1{muon1}, p2{muon2}, gp1{gen_muon1}, gp2{gen_muon2};

    double mass = (p1 + p2).Mag(), gen_mass = (gp1 + gp2).Mag();
    smearing_matrix->Fill(gen_mass, mass);
  }

  // Bayes
  Solvers::Bayes bayes_solver{data, smearing_matrix};
  bayes_solver.InitFlatPrior();
  bayes_solver.Unfold();

  auto result = bayes_solver.GetUnfolded();

  auto outFile = std::make_unique<TFile>("test.root", "recreate");
  outFile->cd();
  data->Write("hMeasured");
  result->Write("hUnfolded_bayes");
  outFile->Close();

  return 0;
}
