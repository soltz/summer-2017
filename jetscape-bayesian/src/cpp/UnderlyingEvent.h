#ifndef ATLASHI_UNDERLYINGEVENT_H__
#define ATLASHI_UNDERLYINGEVENT_H__

#include <vector>
#include <cmath>
#include <set>
#include <iostream>
#include <iomanip>

namespace AtlasHI
{
  class UnderlyingEvent
  {

  public:
    inline unsigned int numEtaBins() const {return m_num_eta_bins;}
    inline float etaMin() const {return m_eta_min;}
    inline float etaMax() const {return m_eta_max;}
    inline unsigned int numPhiBins() const {return m_num_phi_bins;}
    inline float phiMin() const {return m_phi_min;}
    inline float phiMax() const {return m_phi_max;}

    inline float getBinSizeEta() const {return (etaMax()-etaMin())/((float)numEtaBins());}
    inline float getBinSizePhi() const {return (phiMax()-phiMin())/((float)numPhiBins());}

    inline float getBinArea() const {return getBinSizeEta()*getBinSizePhi();}
    inline float getTotalArea() const {return (etaMax()-etaMin())*(phiMax()-phiMin());}
    inline  unsigned int getNumEtaPhiBins() const {return numEtaBins()*numPhiBins();}

    inline float getBinLowEdgeEta(unsigned int eb) const {return etaMin()+eb*getBinSizeEta();}
    inline float getBinUpEdgeEta(unsigned int eb) const {return etaMin()+(eb+1)*getBinSizeEta();}
    inline float getBinCenterEta(unsigned int eb) const {return etaMin()+(eb+0.5)*getBinSizeEta();}

    inline float getBinLowEdgePhi(unsigned int pb) const {return phiMin()+pb*getBinSizePhi();}
    inline float getBinUpEdgePhi(unsigned int pb) const {return phiMin()+(pb+1)*getBinSizePhi();}
    inline float getBinCenterPhi(unsigned int pb) const {return phiMin()+(pb+0.5)*getBinSizePhi();}

    inline unsigned int findBinEta(float eta) const {return std::floor((eta-etaMin())/getBinSizeEta());}
    inline unsigned int findBinPhi(float phi) const 
    {
      int pb=std::floor((phi-phiMin()) /getBinSizePhi());
      return pb % numPhiBins();
    }
    inline unsigned int findBinEtaPhi(float eta, float phi) const { return numPhiBins()*findBinEta(eta)+findBinPhi(phi);}

    UnderlyingEvent() : m_num_eta_bins(100),
			m_num_phi_bins(64),
			m_eta_min(-5.),
			m_eta_max(5.),
			m_phi_min(M_PI),
			m_phi_max(-M_PI),
			m_flow_eta_min(3.2),
			m_flow_eta_max(5.),
			m_flow_abs_eta(true)
    {}

    // UnderlyingEvent(const UnderlyingEvent& u)
    // {
    //   this->m_UE.assign(u.m_UE.begin(),u.m_UE.end());
    //   this->m_ntowers.assign(u.m_ntowers.begin(),u.m_ntowers.end());

    //   this->setEtaBinning(u.m_num_eta_bins,u.m_eta_min,u.m_eta_max);
    //   this->setPhiBinning(u.m_num_phi_bins,u.m_phi_min,u.m_phi_max);
    //   this->setFlowEtaRange(u.m_flow_eta_min,u.m_flow_eta_max);
    //   this->setHarmonics(u.m_harmonics);
    //   m_flow_abs_eta=u.m_flow_abs_eta;
    // }

  private:
    std::vector<float> m_UE;
    std::vector<int> m_ntowers;
    std::vector<float> m_Qx;
    std::vector<float> m_Qy;
    float m_Q_weight;
    std::set<unsigned int> m_harmonics;
    unsigned int m_num_eta_bins;
    unsigned int m_num_phi_bins;
    float m_eta_min;
    float m_eta_max;
    float m_phi_min;
    float m_phi_max;

    float m_flow_eta_min;
    float m_flow_eta_max;
    bool m_flow_abs_eta;

  public:
    const std::set<unsigned int>& getHarmonics() const {return m_harmonics;};
    void setEtaBinning(unsigned int neb, float e1, float e2)
    {
      m_num_eta_bins=neb;
      m_eta_min=e1;
      m_eta_max=e2;
      m_UE.assign(neb,0.);
      m_ntowers.assign(neb,0.);
    }

    void setPhiBinning(unsigned int neb, float e1, float e2)
    {
      m_num_phi_bins=neb;
      m_phi_min=e1;
      m_phi_max=e2;
    }

    void setHarmonics(const std::set<unsigned int>& h)
    {
      m_Qx.assign(h.size(),0);
      m_Qy.assign(h.size(),0);
      m_harmonics=h;
    }

    void setFlowEtaRange(float eta1, float eta2)
    {
      m_flow_eta_min=eta1;
      m_flow_eta_max=eta2;
    }

    float getModulation(float phi) const
    {
      if(m_Q_weight > 0.) return 1.;
      float mod=1.;
      for(auto h : m_harmonics)
      {
	unsigned int ih=h-1;
	// cos(n(phi-Psi_n))= cos (nphi) cos(nPsi_n) + sin (nphi) sin(nPsi_n)
	float fh=static_cast<float>(h);
	mod+=2.0*(std::cos(fh*phi)*m_Qx[ih]+std::sin(fh*phi)*m_Qy[ih])/m_Q_weight;
      }
      return mod;
    }

    float getUEDensity(float eta) const
    {
      unsigned int eb=findBinEta(eta);
      if(m_ntowers[eb] < 1) return 0;
      return m_UE[eb]/static_cast<float>(m_ntowers[eb]);
    }
  
    float getUEDensity(float eta, float phi) const
    {
      return getUEDensity(eta)*getModulation(phi);
    }

    void updateSlice(unsigned int eb, float weight)
    {
      m_UE[eb]+=weight;
      int dn=(weight > 0. ) ? 1 : -1;
      m_ntowers[eb]+=dn;
    }
    void update(float eta, float phi, float weight)
    {
      unsigned int eb=findBinEta(eta);
      updateSlice(eb,weight);
      float eta_v=(m_flow_abs_eta) ? std::abs(eta) : eta;
      if(eta_v >= m_flow_eta_min && eta_v <= m_flow_eta_max)
      {
	m_Q_weight+=weight;
	for(unsigned int ih=0; ih<m_Qx.size(); ih++)
	{
	  float fh=static_cast<float>(ih+1);
	  m_Qx[ih]+=weight*std::cos(fh*phi);
	  m_Qy[ih]+=weight*std::sin(fh*phi);
	}
      }
    }
    void remodulateSlice(unsigned int eb, float sf)
    {
      m_UE[eb]*=sf;
    }

    void print()
    {
      for(unsigned int i=0; i< m_UE.size(); i++)
      {
	std::cout << std::setw(10) << i
		  << std::setw(10) << getBinLowEdgeEta(i)
		  << std::setw(10) << getBinUpEdgeEta(i)
		  << std::setw(15) << m_UE[i]
		  << std::setw(10) << m_ntowers[i]
		  << std::endl;

      }
    }
  };
}
#endif
