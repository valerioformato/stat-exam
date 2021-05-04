// our headers
#include "resolution_model.h"

// ROOT headers
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TH2D.h>
#include <TList.h>
#include <TRandom.h>
#include <TTree.h>

// external headers
#include <fmt/format.h>

// c++ headers
#include <iostream>
#include <memory>
#include <random>

double AngRes(double E) {
  // in mrad
  double t1 = 18.25 / sqrt(E);
  double t2 = 1.1;

  // in rad
  return 1e-3 * sqrt(t1 * t1 + t2 * t2);
}

int main(int argc, char const *argv[]) {
  TH1D *hMass = new TH1D("hMass", "", 200, 0, 10);
  TH2D *hMass_mc = new TH2D("hMass_mc", "", 200, 0, 10, 200, 0, 10);

  auto outFile_data = std::make_unique<TFile>("im_data.root", "recreate");
  TTree *data_tree = new TTree("Data", "CR Data Tree");

  auto outFile_mc = std::make_unique<TFile>("im_mc.root", "recreate");
  TTree *mc_tree = new TTree("Data", "CR MC Tree");

  float meas_muon1[4], meas_muon2[4];
  float gen_muon1[4], gen_muon2[4];

  data_tree->Branch("muon1", meas_muon1, "muon1[4]/F");
  data_tree->Branch("muon2", meas_muon2, "muon2[4]/F");

  mc_tree->Branch("meas_muon1", meas_muon1, "meas_muon1[4]/F");
  mc_tree->Branch("meas_muon2", meas_muon2, "meas_muon2[4]/F");
  mc_tree->Branch("gen_muon1", gen_muon1, "gen_muon1[4]/F");
  mc_tree->Branch("gen_muon2", gen_muon2, "gen_muon2[4]/F");

  float BMass = 5.6;
  float sqrts = 1300.;
  float pMass = 0.938;
  float muMass = 0.1;
  float phRes = 0.02;
  float jRes = 0.3;
  int nJetMean = 2;

  Utils::ResolutionModel model;

  unsigned long long nEvents_data = 10000, nEvents_mc = 100000;
  float bg_frac_data = 0.7;
  // float bg_frac_data = 0.0;

  int nBkg = gRandom->Poisson(nEvents_data * bg_frac_data);
  int nSig = gRandom->Poisson(nEvents_data * (1 - bg_frac_data));

  fmt::print("{} background events\n", nBkg);
  fmt::print("{} signal events\n", nSig);

  TLorentzVector p1, p2;
  p1.SetXYZM(0, 0, sqrts / 2, pMass);
  p2.SetXYZM(0, 0, -sqrts / 2, pMass);
  TLorentzVector W = p1 + p2;

  auto evGen = std::make_unique<TGenPhaseSpace>();

  struct Event {
    std::vector<TLorentzVector> jets;
    std::vector<TLorentzVector> muons;
  };

  std::vector<Event> events;

  // background events
  for (unsigned long long iEv = 0; iEv < nBkg; iEv++) {
    Event event;

    int nJets = gRandom->Poisson(nJetMean);
    if (nJets == 0) {
      iEv--;
      continue;
    }
    event.jets.resize(nJets);

    const int nP = 2 + nJets;
    double masses[nP];
    masses[0] = masses[1] = muMass;
    for (int ij = 0; ij < nJets; ij++)
      masses[2 + ij] = 3.;

    evGen->SetDecay(W, nP, masses);
    evGen->Generate();

    event.muons.push_back(*(evGen->GetDecay(0)));
    event.muons.push_back(*(evGen->GetDecay(1)));

    if ((event.muons[0] + event.muons[1]).Mag() > 10) {
      iEv--;
      continue;
    }

    for (int ij = 0; ij < nJets; ij++)
      event.jets.push_back(*(evGen->GetDecay(2 + ij)));

    for (auto &muon : event.muons) {
      float meas_momentum = model(muon.Vect().Mag());
      // muon.SetTheta(muon.Theta() * gRandom->Gaus(1, AngRes(muon.E())));
      muon.SetE(sqrt(muMass * muMass + meas_momentum * meas_momentum));
    }
    for (auto &jet : event.jets) {
      jet.SetE(jet.E() * gRandom->Gaus(1, jRes));
    }

    events.push_back(event);

    hMass->Fill((event.muons[0] + event.muons[1]).Mag());
  }

  // signal events
  for (unsigned long long iEv = 0; iEv < nSig; iEv++) {
    Event event;

    int nJets = gRandom->Poisson(nJetMean);
    if (nJets == 0) {
      iEv--;
      continue;
    }
    event.jets.resize(nJets);

    const int nP = 1 + nJets;
    double masses[nP];
    masses[0] = BMass;
    for (int ij = 0; ij < nJets; ij++)
      masses[1 + ij] = 3.;

    evGen->SetDecay(W, nP, masses);
    evGen->Generate();

    for (int ij = 0; ij < nJets; ij++)
      event.jets.push_back(*(evGen->GetDecay(1 + ij)));

    double dummy_mass[2] = {muMass, muMass};
    evGen->SetDecay(*(evGen->GetDecay(0)), 2, dummy_mass);
    evGen->Generate();

    event.muons.push_back(*(evGen->GetDecay(0)));
    event.muons.push_back(*(evGen->GetDecay(1)));

    if ((event.muons[0] + event.muons[1]).Mag() > 10) {
      iEv--;
      continue;
    }

    for (auto &muon : event.muons) {
      float meas_momentum = model(muon.Vect().Mag());

      // fmt::print("{} -> {}\n", muon.Vect().Mag(), meas_momentum);

      // muon.SetTheta(muon.Theta() * gRandom->Gaus(1, AngRes(muon.E())));
      muon.SetE(sqrt(muMass * muMass + meas_momentum * meas_momentum));
    }
    for (auto &jet : event.jets) {
      jet.SetE(jet.E() * gRandom->Gaus(1, jRes));
    }

    events.push_back(event);

    hMass->Fill((event.muons[0] + event.muons[1]).Mag());
  }

  std::shuffle(begin(events), end(events), std::mt19937{std::random_device{}()});

  for (const auto &event : events) {
    event.muons.at(0).GetXYZT(meas_muon1);
    event.muons.at(1).GetXYZT(meas_muon2);
    data_tree->Fill();
  }

  // do MC
  events.clear();
  std::vector<Event> gen_events;

  for (unsigned long long iEv = 0; iEv < nEvents_mc; iEv++) {
    Event event;

    int nJets = gRandom->Poisson(nJetMean);
    if (nJets == 0) {
      iEv--;
      continue;
    }
    event.jets.resize(nJets);

    const int nP = 2 + nJets;
    double masses[nP];
    masses[0] = masses[1] = muMass;
    for (int ij = 0; ij < nJets; ij++)
      masses[2 + ij] = 3.;

    evGen->SetDecay(W, nP, masses);
    evGen->Generate();

    event.muons.push_back(*(evGen->GetDecay(0)));
    event.muons.push_back(*(evGen->GetDecay(1)));

    double gen_mass = (event.muons[0] + event.muons[1]).Mag();

    if (gen_mass > 10) {
      iEv--;
      continue;
    }

    for (int ij = 0; ij < nJets; ij++)
      event.jets.push_back(*(evGen->GetDecay(2 + ij)));

    gen_events.push_back(event);

    for (auto &muon : event.muons) {
      float meas_momentum = model(muon.Vect().Mag());
      // muon.SetTheta(muon.Theta() * gRandom->Gaus(1, AngRes(muon.E())));
      muon.SetE(sqrt(muMass * muMass + meas_momentum * meas_momentum));
    }
    for (auto &jet : event.jets) {
      jet.SetE(jet.E() * gRandom->Gaus(1, jRes));
    }

    double meas_mass = (event.muons[0] + event.muons[1]).Mag();

    events.push_back(event);

    hMass_mc->Fill(gen_mass, meas_mass);
  }

  for (unsigned int iEv = 0; iEv < events.size(); iEv++) {
    events[iEv].muons.at(0).GetXYZT(meas_muon1);
    events[iEv].muons.at(1).GetXYZT(meas_muon2);
    gen_events[iEv].muons.at(0).GetXYZT(gen_muon1);
    gen_events[iEv].muons.at(1).GetXYZT(gen_muon2);
    mc_tree->Fill();
  }

  outFile_data->cd();
  data_tree->Write();
  // hMass->Write();

  outFile_mc->cd();
  mc_tree->Write();
  // hMass_mc->Write();

  return 0;
}
