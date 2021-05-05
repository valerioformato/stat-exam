// our headers
#include "localization.h"
#include "resolution_model.h"
#include "spectrum.h"

// ROOT headers
#include <TFile.h>
#include <TH2D.h>
#include <TList.h>
#include <TRandom.h>
#include <TTree.h>

// external headers
#include <fmt/format.h>

// c++ headers
#include <iostream>
#include <memory>

int main(int argc, char const *argv[]) {
  auto outFile_data = std::make_unique<TFile>("cr_data.root", "recreate");
  TTree *data_tree = new TTree("Data", "CR Data Tree");

  auto outFile_mc = std::make_unique<TFile>("cr_mc.root", "recreate");
  TTree *mc_tree = new TTree("Data", "CR MC Tree");

  float meas_rigidity = 0, gen_rigidity = 0;
  data_tree->Branch("rigidity", &meas_rigidity);
  mc_tree->Branch("meas_rigidity", &meas_rigidity);
  mc_tree->Branch("gen_rigidity", &gen_rigidity);

  auto spectrum = Utils::GetSpectrum(1, 1);

  auto livetimeFile =
      std::make_unique<TFile>((localization::PROJECT_SOURCE_DIR + "/datafiles/cr_flux/livetime.root").c_str(), "open");
  TH1D *livetime = static_cast<TH2D *>(livetimeFile->Get<TList>("ListLiveTime_rbinv5")->At(1))->ProjectionY();
  livetime->Scale(1.0 / livetime->GetBinContent(livetime->GetNbinsX()));

  auto matFile = std::make_unique<TFile>(
      (localization::PROJECT_SOURCE_DIR + "/datafiles/cr_flux/res_matrix.root").c_str(), "open");
  Utils::ResolutionModel model{matFile->Get<TList>("parList_FS")};

  unsigned long long nEvents_data = 1000000, nEvents_mc = 5000000;

  for (unsigned long long iEv = 0; iEv < nEvents_data; iEv++) {
    // generate rigidity according to some spectrum
    gen_rigidity = spectrum->GetRandom();

    if (gRandom->Uniform() > livetime->Interpolate(gen_rigidity)) {
      iEv--;
      continue;
    }

    meas_rigidity = model(gen_rigidity);

    // fmt::print("Gen = {}, Meas = {}\n", gen_rigidity, meas_rigidity);

    data_tree->Fill();
  }

  for (unsigned long long iEv = 0; iEv < nEvents_mc; iEv++) {
    // generate rigidity according to some spectrum
    gen_rigidity = spectrum->GetRandom();
    meas_rigidity = model(gen_rigidity);

    // fmt::print("Gen = {}, Meas = {}\n", gen_rigidity, meas_rigidity);

    mc_tree->Fill();
  }

  outFile_data->cd();
  livetime->Write("LiveTime");
  data_tree->Write();

  outFile_mc->cd();
  mc_tree->Write();

  return 0;
}
