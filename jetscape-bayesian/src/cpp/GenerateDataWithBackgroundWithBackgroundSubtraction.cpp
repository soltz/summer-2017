#include "Background.h"
#include "Event.h"
#include "Utilities.h"
#include "HIJetFinder.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include "HepMC/GenEvent.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

// g++ GenerateDataWithBackgroundWithBackgroundSubtraction.cpp HIJetFinder.cc Background.cpp Event.cpp Utilities.cpp -lPythia8 -std=c++11 -lfastjet -lHepMC -o bgsub.out

class EventData {

private:

  double quench;
  int mult;
  double xjg;

public:

  EventData(double a, int b, double c) {
    quench = a;
    mult = b;
    xjg = c;
  }

  std::string to_string() {
    return std::to_string(quench) + " " + std::to_string(mult) + " " + std::to_string(xjg) + "\n";
  }

};

int main(void) {
  // Set options
  std::vector<std::string> options;
  options.push_back("PromptPhoton:all = off");
  options.push_back("HardQCD:all = off");
  options.push_back("PromptPhoton:qg2qgamma = on");
  options.push_back("PromptPhoton:qqbar2ggamma = on");
  options.push_back("TimeShower:QEDshowerByGamma = off");
  options.push_back("PhaseSpace:pTHatMin = 15.0");
  // Init event generator
  EventGenerator c(200.0, options);
  BackgroundGenerator a("Au", "Au", 200.0);
  // Loop over number of events
  save_pdg_data();
  std::ofstream f("bg_data3.txt");

  // Create HepMC obj for data
  HepMC::GenEvent* jet_evt = new HepMC::GenEvent();
  HepMC::GenEvent* bg_evt = new HepMC::GenEvent();

  // HIJet init
  AtlasHI::HIJetFinder HI;
  AtlasHI::UE_t u;
  u.setEtaBinning(50, -2.5, 2.5);
  u.setPhiBinning(64,-1 * M_PI, M_PI);                                                                                                                                                                                                                                                                    
  u.setFlowEtaRange(2, 2.5);
  HI.configureUE(u);

  //radiius for jet clustering, radius for “seeds”, radius around seeds to
  // exclude from background, minimum value of max tower to be a seed,
  // minimum value of max/mean tower to be a seed, minimum jet pT to be a
  // seed after first round of subtraction
  HI.setParameters(0.5, 0.2, 0.6, 1, 3, 7);

  for (double quench = 1.0; quench >= 0.13; quench -= 0.05) {
    std::cout << "quench = " << quench << "\n";
    std::vector<EventData> d;
    d.reserve(100000);
    for (int i = 0; i < 100000; i++) {
      if (i % 10000 == 0) {
        std::cout << "i = " << i << "\n";
      }
      c.set_quench(quench);
      
      // Get data
      c.generate(jet_evt);
      a.generate(bg_evt);
      // Prep event for past jet
      std::vector<fastjet::PseudoJet> p = merge_and_prep(jet_evt, bg_evt);
      int mult = p.size();
      // Get max pt photon
      // double pure_pt = tag_max_pt_photon(p);
      fastjet::PseudoJet photon = get_max_pt_photon(p);
      double pure_pt = photon.pt();
      while ((pure_pt < 20.0) || (pure_pt > 25.0)) {
        jet_evt->clear();
        c.generate(jet_evt);
        p = merge_and_prep(jet_evt, bg_evt);
        mult = p.size();
        photon = get_max_pt_photon(p);
        pure_pt = photon.pt();
      }
      std::vector<fastjet::PseudoJet> reco_jets;
      HI.processEvent(p, reco_jets);
      std::vector<fastjet::PseudoJet> sortedJets = fastjet::sorted_by_pt(reco_jets);
      int photon_location = -1;
      double delta_r = 400.0;
      for (int j = 0; j < sortedJets.size(); j++) {
        double temp = photon.delta_R(sortedJets[j]);
        if (temp < delta_r) {
          delta_r = temp;
          photon_location = j;
        }
      }
      // std::cout << delta_r << " " << photon_location << std::endl;

      // Sort jets by PT
      // std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(5.0));

      // for (int i = 0; i < sortedJets[0].constituents().size(); i++) {
      //   std::cout << sortedJets[0].constituents()[i].user_index() << "\n"; 
      // }

      // Locate highest PT non-photon jet
      // Search for photon in max-Pt jet
      // bool found = contains_tagged_photon(clust_seq, inclusive_jets[0]);
      bool found = true;
      if (photon_location != 0) {
        found = false;
      }
      double quench_pt = 0.0;
      if (found) {
        quench_pt = sortedJets[1].pt();
      } else {
        quench_pt = sortedJets[0].pt();
      }
      double xjg;
      if (sortedJets.size() < 2) {
        xjg = 0.0;
      } else {
        xjg = quench_pt / pure_pt;
      }
      d.push_back(EventData(quench, mult, xjg));
      jet_evt->clear();
      bg_evt->clear();
    }
    for (int i = 0; i < d.size(); i++) {
      f << d[i].to_string();
    }
  }
  f.close();
}
