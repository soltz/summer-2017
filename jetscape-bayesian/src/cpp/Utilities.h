#include <vector>
#include "HepMC/GenEvent.h"
#include <fastjet/PseudoJet.hh>

// Methods to convert HepMC event to a PseudoJet vector for fastjet
std::vector<fastjet::PseudoJet> prep_single_event(HepMC::GenEvent* evt);
std::vector<fastjet::PseudoJet> merge_and_prep(HepMC::GenEvent* evt1,
                                               HepMC::GenEvent* evt2);

// Methods to tag and locate unquenched hard photon
double tag_max_pt_photon(std::vector<fastjet::PseudoJet>& p);
