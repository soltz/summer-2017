#include "Event.h"
#include "Utilities.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include "HepMC/GenEvent.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

// g++ GenerateDataOriginal.cpp Background.cpp Event.cpp Utilities.cpp -lPythia8 -std=c++11 -lfastjet -lHepMC -o orig.out

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
  options.push_back("PhaseSpace:pTHatMin = 20.0");
  options.push_back("PhaseSpace:pTHatMax = 25.0");
  // Init event generator
  EventGenerator c(200.0, options);
  // Loop over number of events
  std::ofstream f("orig_data.txt");

  // Create HepMC obj for data
  HepMC::GenEvent* jet_evt = new HepMC::GenEvent();
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
      // Prep event for past jet
      std::vector<fastjet::PseudoJet> p = prep_single_event(jet_evt);
      int mult = p.size();
      // Get max pt photon
      double pure_pt = tag_max_pt_photon(p);
      // while ((pure_pt < 20.0) || (pure_pt > 25.0)) {
      //   jet_evt->clear();
      //   c.generate(jet_evt);
      //   p = prep_single_event(jet_evt);
      //   mult = p.size();
      //   pure_pt = tag_max_pt_photon(p);
      // }
      // Set up fastjet
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.5);
      // Run fastjet
      fastjet::ClusterSequence clust_seq(p, jet_def);
      // Sort jets by PT
      std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(5.0));

      // Locate highest PT non-photon jet
      // Search for photon in max-Pt jet
      bool found = contains_tagged_photon(clust_seq, inclusive_jets[0]);
      double quench_pt = 0.0;
      if (found) {
        quench_pt = inclusive_jets[1].pt();
      } else {
        quench_pt = inclusive_jets[0].pt();
      }
      double xjg;
      if (inclusive_jets.size() < 2) {
        xjg = 0.0;
      } else {
        xjg = quench_pt / pure_pt;
      }
      d.push_back(EventData(quench, mult, xjg));
      jet_evt->clear();
    }
    for (int i = 0; i < d.size(); i++) {
      f << d[i].to_string();
    }
  }
  f.close();
}
