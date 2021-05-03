// our headers
#include "resolution_model.h"
#include "spectrum.h"

// ROOT headers
#include <TFile.h>
#include <TTree.h>

// c++ headers
#include <memory>

int main(int argc, char const *argv[]) {
  auto outFile = std::make_unique<TFile>("cr_sample.root", "recreate");

  TTree *sample_tree = new TTree("Data", "CR Sample Tree");

  float rigidity = 0;
  sample_tree->Branch("rigidity", &rigidity);

  auto spectrum = Utils::GetSpectrum(1, 1);

  auto matFile = std::make_unique<TFile>("res_matrix.root", "open");
  Utils::ResolutionModel model{matFile->Get<TList>("parList_FS")};

  unsigned long long nEvents = 10000000;
  for (unsigned long long iEv = 0; iEv < nEvents; iEv++) {
    // generate rigidity according to some spectrum
    float trueRigidity = spectrum->GetRandom();

    sample_tree->Fill();
  }

  outFile->cd();
  spectrum.get()->Write();
  sample_tree->Write();
  return 0;
}
