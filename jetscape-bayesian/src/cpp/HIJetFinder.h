#ifndef ATLASHI_HIJETFINDER_H__
#define ATLASHI_HIJETFINDER_H__


#include "UnderlyingEvent.h"
#include <fastjet/PseudoJet.hh>
#include <vector>
#include <iostream>

namespace AtlasHI
{
  typedef UnderlyingEvent UE_t;

  class HIJetFinder
  {
  public:

    HIJetFinder();
    void processEvent(const std::vector<fastjet::PseudoJet>& fjInputs, std::vector<fastjet::PseudoJet>& outputs);

    void applySubtractionToJet(fastjet::PseudoJet& j, const UE_t& u);
    void buildSeedsPtThreshold(const std::vector<fastjet::PseudoJet>& jets, std::vector<fastjet::PseudoJet>& seeds);
    void buildSeedsMaxOverMean(const std::vector<fastjet::PseudoJet>& jets, std::vector<fastjet::PseudoJet>& seeds);
    void updateUEFromSeeds( UE_t& u, const std::vector<fastjet::PseudoJet>& seeds, const std::vector<fastjet::PseudoJet>& the_towers);
    float getUEFromParticle(const fastjet::PseudoJet& p);
    void configureUE(const UE_t& u){m_UE=u;};
    void setParameters(float cluster_R, float seed_R, float exclusion_DR, float cmax, float cdiscrim, float iter_pt)
    {
      m_cluster_Rval=cluster_R;
      m_seed_Rval=seed_R;
      m_seed_DR=exclusion_DR;
      m_max_over_mean_cut=cdiscrim;
      m_max_constit_cut=cmax;
      m_iteration_seed_pt_cut=iter_pt;
    }
    void print() const;
    enum TowerScheme{FOUR_VECTOR=0, MASSLESS};
    enum UEScheme{UE_ET=0, UE_E, UE_PT};
    enum TowerDRMode{PHYSICAL=0, FIXED, BINWISE};
  
  private:

    TowerScheme m_tower_scheme;
    UEScheme m_ue_scheme;
    TowerDRMode m_tower_DR_mode;
    float m_pt_min;
    float m_cluster_Rval;
    float m_seed_Rval;
    float m_seed_DR;
    float m_iteration_seed_pt_cut;
    float m_max_over_mean_cut;
    float m_max_constit_cut;
    UE_t m_UE;
    float m_ghost_scale;
    bool m_exclude_empty_towers;
  };

  //keep track of per jet subtracted quantities using
  //fastjet's UserInfoBase
  //requires implementation of helper class
  //not strictly needed for the algorithm, but can give useful additional info
  class HIJetUserInfo : public fastjet::PseudoJet::UserInfoBase
  {
  public:
    HIJetUserInfo(){};
    ~HIJetUserInfo(){};
    fastjet::PseudoJet unsubtracted;
  };
}
#endif
