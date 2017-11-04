#include "Background.h"
#include "Utilities.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include "HepMC/GenEvent.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

int main(void) {
  // Init background generator
  BackgroundGenerator a("Au", "Au", 200.0);

  save_pdg_data();
  std::ofstream f("bg_data.txt");

  // Create HepMC obj for background only
  HepMC::GenEvent* bg_evt = new HepMC::GenEvent();
  
  // Loop over number of events
  for (int i = 0; i < 1000; i++) {
    if (i % 100 == 0) {
      std::cout << "i = " << i << "\n";
    }
      
    // Get data
    // Note, to save full HEPMC event, call without arguments
    //    a.generate();
    a.generate(bg_evt);

    // Prep event for past jet
    std::vector<fastjet::PseudoJet> p = prep_single_event(bg_evt);

    // Loop over particles in event
    std::cout << p.size() << std::endl;
    for (int j = 0; j < p.size(); j++) {
      // double theta = 2 * atan(exp(-1 * p[j].pseudorapidity()));
      f << p[j].phi_std() << " " << p[j].rap() << " " << p[j].pt() << std::endl;
    }
    bg_evt->clear();
  }
  f.close();
}
