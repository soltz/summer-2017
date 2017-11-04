#include "HIJetFinder.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SharedPtr.hh>
#include <iomanip>

namespace AtlasHI
{
  HIJetFinder::HIJetFinder() : m_tower_scheme(MASSLESS),
			       m_ue_scheme(UE_ET),
			       m_tower_DR_mode(PHYSICAL),
			       m_pt_min(5.),
			       m_cluster_Rval(0.4),
			       m_seed_Rval(0.2),
			       m_seed_DR(0.6),
			       m_iteration_seed_pt_cut(20.),
			       m_max_over_mean_cut(3.),
			       m_max_constit_cut(3.),
			       m_UE(),
			       m_ghost_scale(1e-12),
			       m_exclude_empty_towers(false)
  {

  }

  void HIJetFinder::processEvent(const std::vector<fastjet::PseudoJet>& fjInputs, std::vector<fastjet::PseudoJet>& outputs)
  {
    //set up the UE binning
    UE_t UE_unsub(m_UE);

    //fill the towers and UE from the input particles
    unsigned int num_ep_bins=m_UE.getNumEtaPhiBins();
    std::vector<fastjet::PseudoJet> the_towers(num_ep_bins,fastjet::PseudoJet(m_ghost_scale,m_ghost_scale,m_ghost_scale,4*m_ghost_scale));
    std::vector<unsigned int> tower_mult(num_ep_bins,0);
    std::vector<float> sum_w(num_ep_bins,0);
    std::vector<float> sum_eta(num_ep_bins,0);
    std::vector<float> sum_phi(num_ep_bins,0);
    if(m_tower_scheme==MASSLESS)
    {
      sum_eta.assign(num_ep_bins,0);
      sum_phi.assign(num_ep_bins,0);
    }
    for(auto p : fjInputs)
    {
      float eta=p.eta();
      float phi=p.phi();
      unsigned int epb=UE_unsub.findBinEtaPhi(eta,phi);
      //      if(epb < 0 || epb >= num_ep_bins) continue;
      if(epb >= num_ep_bins) continue;
      the_towers[epb]+=p;
      tower_mult[epb]++;
      if(m_tower_scheme==MASSLESS)
      {
	float weight=getUEFromParticle(p);
	sum_w.at(epb)+=weight;
	sum_eta.at(epb)+=eta*weight;
	sum_phi.at(epb)+=phi*weight;
      }

    }
    //degeneracy in DR can cause problems for fastjet
    //FP rounding issues may make clustering sequence irreproducible or may hurt performance
    //empty towers must be kept for indexing, but they have no effect on results
    //fix used here is arbitrary but ensures unique DR for zero energy towers and is reproducible
    float eta_shift_factor=1e-3*UE_unsub.getBinSizeEta()/static_cast<float>(UE_unsub.getNumEtaPhiBins());

    std::vector<fastjet::PseudoJet> the_filled_towers;
    the_filled_towers.reserve(the_towers.size());
    for(unsigned int itow=0; itow<the_towers.size(); itow++)
    {
      fastjet::PseudoJet& tower=the_towers[itow];
      if(tower_mult[itow]==0)
      {
	float eta=UE_unsub.getBinCenterEta(itow / UE_unsub.numPhiBins());
	float phi=UE_unsub.getBinCenterPhi(itow % UE_unsub.numPhiBins());
	//avoid degeneracy by shifting eta by a small (<1e-3 x bin width), unique amount
	eta+=eta_shift_factor*itow;
	float ET=m_ghost_scale/std::cosh(eta);
	tower.reset_momentum(ET*std::cos(phi),ET*std::sin(phi),ET*std::sinh(eta),m_ghost_scale);

      }
      else if(m_tower_scheme==MASSLESS && sum_w[itow] > 0.)
      {
	float energy=tower.e();
	float weight=sum_w[itow];
	float eta=sum_eta[itow]/weight;
	float phi=sum_phi[itow]/weight;
	float ET=energy/std::cosh(eta);
	tower.reset_momentum(ET*std::cos(phi),ET*std::sin(phi),ET*std::sinh(eta),energy);
      }
      if(tower_mult[itow]!=0 && m_exclude_empty_towers) the_filled_towers.push_back(tower);
      UE_unsub.update(tower.eta(),tower.phi(),getUEFromParticle(tower));
    }

    //build the seed container
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, m_seed_Rval,fastjet::E_scheme,fastjet::Best);
    std::vector<fastjet::PseudoJet>& pj_vec=(m_exclude_empty_towers) ? the_filled_towers : the_towers;
    fastjet::ClusterSequence clustSeq(pj_vec, jetDef);
    auto incjets=clustSeq.inclusive_jets(m_pt_min);
    std::vector<fastjet::PseudoJet> initial_jets=fastjet::sorted_by_pt(incjets);
  
    //first iteration
    //using max/mean
    UE_t UE_sub1(UE_unsub);

    std::vector<fastjet::PseudoJet> seeds1;
    buildSeedsMaxOverMean(initial_jets,seeds1);
    updateUEFromSeeds(UE_sub1,seeds1,the_towers);
    //subtract
    std::vector<fastjet::PseudoJet> sub1_jets(initial_jets);
    for(auto itr=sub1_jets.begin(); itr!=sub1_jets.end(); itr++) applySubtractionToJet(*itr,UE_sub1);
 
    //second iteration
    std::vector<fastjet::PseudoJet> seeds2;
    buildSeedsPtThreshold(sub1_jets,seeds2);
    UE_t UE_sub2(UE_unsub);
    updateUEFromSeeds(UE_sub2,seeds2,the_towers);

    fastjet::JetDefinition jetDef_out(fastjet::antikt_algorithm, m_cluster_Rval,fastjet::E_scheme,fastjet::Best);
    fastjet::ClusterSequence clustSeq_out(pj_vec, jetDef_out);

    if(m_seed_Rval==m_cluster_Rval) outputs.assign(initial_jets.begin(),initial_jets.end());
    else outputs=fastjet::sorted_by_pt(clustSeq_out.inclusive_jets(m_pt_min));
    for(auto itr=outputs.begin(); itr!=outputs.end(); itr++) applySubtractionToJet(*itr,UE_sub2);
  }

  void HIJetFinder::applySubtractionToJet(fastjet::PseudoJet& j, const UE_t& u)
  {

    HIJetUserInfo* info_p=new HIJetUserInfo();
    j.set_user_info(info_p);
    fastjet::PseudoJet& unsubtracted=info_p->unsubtracted;
    j.reset_momentum(0,0,0,0);

    for(auto p : j.constituents())
    {
      float ue_contrib=u.getUEDensity(p.eta(),p.phi());
      float particle_contrib=getUEFromParticle(p);
      float sf=1.0-ue_contrib/particle_contrib;
      fastjet::PseudoJet subtracted_constit(p);
      //special fix for neg E?
      unsubtracted+=p;
      subtracted_constit*=sf;
      j+=subtracted_constit;
    }
  }
  void HIJetFinder::buildSeedsPtThreshold(const std::vector<fastjet::PseudoJet>& jets, std::vector<fastjet::PseudoJet>& seeds)
  {
    seeds.reserve(jets.size());
    for(auto j : jets)
    {
      if(j.pt() > m_iteration_seed_pt_cut) seeds.push_back(j);
    }
  }
  void HIJetFinder::buildSeedsMaxOverMean(const std::vector<fastjet::PseudoJet>& jets, std::vector<fastjet::PseudoJet>& seeds)
  {
    seeds.reserve(jets.size());
    for(auto j : jets)
    {
      float max_constit_weight=0;
      float mean_constit_weight=0;
      auto constit=j.constituents();
      for(auto p : constit)
      {
	float ue_contrib=getUEFromParticle(p);
	mean_constit_weight+=ue_contrib;
	if(ue_contrib > max_constit_weight) max_constit_weight=ue_contrib;
      }
      if(constit.size()==0) continue;
      mean_constit_weight=mean_constit_weight/static_cast<float>(constit.size());
      if(max_constit_weight<m_max_constit_cut) continue;
      if(max_constit_weight/mean_constit_weight < m_max_over_mean_cut) continue;

      //if it survives, its a seed
      seeds.push_back(j);
    }
  }
  void HIJetFinder::updateUEFromSeeds(UE_t& u, const std::vector<fastjet::PseudoJet>& seeds, const std::vector<fastjet::PseudoJet>& the_towers)
  {

    std::set<unsigned int> excluded_towers;
    int Nr=m_seed_DR/m_UE.getBinSizeEta();
    float jsf=(m_tower_DR_mode==PHYSICAL)?  m_UE.getBinSizePhi()/m_UE.getBinSizeEta() : 1.0;
    for(auto s : seeds)
    {
      unsigned int i_center=m_UE.findBinEta(s.eta());
      unsigned int j_center=m_UE.findBinPhi(s.phi_std());
      for(int j=-Nr; j<=Nr; j++)
      {
    	int j_prime=j_center+j;
    	//handle wrap around
    	if(j_prime < 0) j_prime=m_UE.numPhiBins()+j_prime;
    	else if(j_prime >= static_cast<int>(m_UE.numPhiBins())) j_prime=j_prime-m_UE.numPhiBins();
    	int i_max=(m_tower_DR_mode==FIXED) ? Nr-std::abs(j) : std::sqrt(Nr*Nr-j*j*jsf);

    	for(int i=-i_max; i<=i_max;i++)
    	{
    	  int i_prime=i_center+i;
    	  if(i_prime < 0) continue;
    	  if(i_prime >= static_cast<int>(m_UE.numEtaBins())) continue;

    	  unsigned int index=i_prime*m_UE.numPhiBins()+j_prime;

    	  if(excluded_towers.emplace(index).second) //returns true and inserts if index was not already in set
    	    u.update(the_towers[index].eta(),the_towers[index].phi(),-1.0*getUEFromParticle(the_towers[index]));

    	}
      }
    }

    // for(auto s : seeds)
    // {
    //   for(unsigned int index=0; index<the_towers.size();index++)
    //   {
    // 	const fastjet::PseudoJet& tower=the_towers[index];
    // 	float dphi=tower.delta_phi_to(s);
    // 	float deta=tower.eta()-s.eta();
    // 	if(deta*deta+dphi*dphi < m_seed_DR*m_seed_DR)
    // 	{
    // 	  std::cout << std::setw(12) << s.pt()
    // 		    << std::setw(12) << s.eta()
    // 		    << std::setw(12) << s.phi_std() << "|"
    // 		    << std::setw(12) << tower.pt()
    // 		    << std::setw(12) << tower.eta()
    // 		    << std::setw(12) << tower.phi_std() << "|"
    // 		    << std::setw(8)  << index
    // 		    << std::setw(12)  << std::sqrt(deta*deta+dphi*dphi)
    // 		    << std::endl;
    // 	  if(excluded_towers.emplace(index).second) //returns true and inserts if index was not already in set
    // 	    u.update(the_towers[index].eta(),the_towers[index].phi(),-1.0*getUEFromParticle(the_towers[index]));
    // 	}

    //   }
    // }
    //remodulate

    std::vector<float> mod_factors(u.numEtaBins(),0);
    std::vector<float> mod_counts(u.numEtaBins(),0);

    for(auto index : excluded_towers)
    {
      unsigned int phi_bin=index % u.numPhiBins();
      unsigned int eta_bin=index / u.numPhiBins();
      mod_counts[eta_bin]++;
      mod_factors[eta_bin]+=(u.getModulation(u.getBinCenterPhi(phi_bin))-1.);
    }
    double nphibins=u.numPhiBins();
    for(unsigned int eb=0; eb<u.numEtaBins(); eb++)
    {
      double neff=nphibins-mod_counts[eb];
      double cf=neff/(neff-mod_factors[eb]);
      //check on value of cf;
      if(cf < 0.) cf =1;
      if(cf > 2.) cf=2;
      u.remodulateSlice(eb,cf);
    }
  }

  float HIJetFinder::getUEFromParticle(const fastjet::PseudoJet& p)
  {
    if(m_ue_scheme==UE_ET) return p.Et();
    else if(m_ue_scheme==UE_E) return p.E();
    return p.pt();
  }

  void HIJetFinder::print() const
  {
    std::cout << "Configuring HIJetFinder with settings "  << std::endl;
    std::cout << "eta: " << m_UE.numEtaBins()
	      << " bins on range "
	      << m_UE.etaMin() << " - " << m_UE.etaMax() 
	      << "; "
	      << "phi: " << m_UE.numPhiBins()
	      << " bins on range "
	      << m_UE.phiMin() << " - " << m_UE.phiMax() 
	      << std::endl;
    if(m_UE.getHarmonics().size()==0) std::cout << "No harmonic modulation" << std::endl;
    else
    {
      std::cout << "Using harmonics:";
      for(auto h : m_UE.getHarmonics()) std::cout << " " << h;
      std::cout << std::endl;
    }
  }
}
