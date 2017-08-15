#include "Utilities.h"
#include <vector>
#include "HepMC/GenEvent.h"
#include "HepMC/SimpleVector.h"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

std::vector<fastjet::PseudoJet> prep_single_event(HepMC::GenEvent* evt) {

  std::vector<fastjet::PseudoJet> particles;

  // Iterate over all particles in evt
  for(HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
      p != evt->particles_end(); ++p) {

    // Push final state particles into vector
    if ((*p)->is_undecayed()) {
      HepMC::FourVector fv = (*p)->momentum();
      fastjet::PseudoJet q = fastjet::PseudoJet(fv.px(), fv.py(),
                                                fv.pz(), fv.e());
      // Set user_index to particle_id for later use
      q.set_user_index((*p)->pdg_id());
      particles.push_back(fastjet::PseudoJet(q));
    }
  }

  return particles;
}

std::vector<fastjet::PseudoJet> merge_and_prep(HepMC::GenEvent* evt1,
                                               HepMC::GenEvent* evt2) {

  std::vector<fastjet::PseudoJet> particles;

  // Iterate over all particles in first event
  for(HepMC::GenEvent::particle_const_iterator p = evt1->particles_begin();
      p != evt1->particles_end(); ++p) {

    // Push final state particles into vector
    if ((*p)->is_undecayed()) {
      HepMC::FourVector fv = (*p)->momentum();
      fastjet::PseudoJet q = fastjet::PseudoJet(fv.px(), fv.py(),
                                                fv.pz(), fv.e());
      // Set user_index to particle_id for later use
      q.set_user_index((*p)->pdg_id());
      particles.push_back(fastjet::PseudoJet(q));
    }
  }

  // Iterate over all particles in second event
  for(HepMC::GenEvent::particle_const_iterator p = evt2->particles_begin();
      p != evt2->particles_end(); ++p) {

    // Push final state particles into vector
    if ((*p)->is_undecayed()) {
      HepMC::FourVector fv = (*p)->momentum();
      fastjet::PseudoJet q = fastjet::PseudoJet(fv.px(), fv.py(),
                                                fv.pz(), fv.e());
      // Set user_index to particle_id for later use
      q.set_user_index((*p)->pdg_id());
      particles.push_back(fastjet::PseudoJet(q));
    }
  }

  return particles;
}

double tag_max_pt_photon(std::vector<fastjet::PseudoJet>& p) {
  int index = -1;
  double pt = 0;

  // Iterate over PseudoJets in vector
  for (int i = 0; i < p.size(); i++) {
    // If harder photon located, tag it
    if ((p[i].user_index() == 22) && (p[i].pt() > pt)) {
      // Reset previous max photon user index
      if (index != -1) {
        p[index].set_user_index(0);
      }
      // Flag current max photon user index
      p[i].set_user_index(-1);
      index = i;
      pt = p[i].pt();
    } else {
      // If not photon, reset its user index
      p[i].set_user_index(0);
    }
  }

  // Return hard photon pT or 0 if none found
  if (index != -1) {
    return pt;
  } else {
    return 0.0;
  }
}

bool contains_tagged_photon(const fastjet::ClusterSequence& clust,
                            const fastjet::PseudoJet& jet) {
  std::vector<fastjet::PseudoJet> jet_const = clust.constituents(jet);
  bool found = false;
  for (int j = 0; j < jet_const.size(); j++) {
    if (jet_const[j].user_index() == -1) {
      found = true;
    }
  }
  return found;
}

int locate_tagged_photon(const fastjet::ClusterSequence& clust,
                         const std::vector<fastjet::PseudoJet>& jets) {
  int index = -1;
  for (int i = 0; i < jets.size(); i++) {
    if (contains_tagged_photon(clust, jets[0])) {
      index = i;
      break;
    }
  }
  return index;
}
