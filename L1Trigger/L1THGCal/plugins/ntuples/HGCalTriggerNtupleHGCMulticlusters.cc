#include <math.h>
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1THGCal/interface/ClusterShapes.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerNtupleBase.h"



class HGCalTriggerNtupleHGCMulticlusters : public HGCalTriggerNtupleBase
{

  public:
    HGCalTriggerNtupleHGCMulticlusters(const edm::ParameterSet& conf);
    ~HGCalTriggerNtupleHGCMulticlusters(){};
    virtual void initialize(TTree&, const edm::ParameterSet&, edm::ConsumesCollector&&) override final;
    virtual void fill(const edm::Event& e, const edm::EventSetup& es) override final;

  private:
    virtual void clear() override final;


    edm::EDGetToken multiclusters_token_;

    int cl3d_n_ ;
    std::vector<float> cl3d_pt_;
    std::vector<float> cl3d_energy_;
    std::vector<float> cl3d_eta_;
    std::vector<float> cl3d_phi_;
    std::vector< std::vector<float> > cl3d_cl_mipPt_;
    std::vector< std::vector<int> > cl3d_cl_layer_;
    std::vector<float> cl3d_x_;
    std::vector<float> cl3d_y_;
    std::vector<float> cl3d_z_;
    std::vector<float> cl3d_xN_;
    std::vector<float> cl3d_yN_;
    std::vector<int> cl3d_nclu_;

    std::vector<int> cl3d_C2dLayer1stMaxE_;
    std::vector<int> cl3d_C2dLayer2ndMaxE_;
    std::vector<int> cl3d_C2dLayer3rdMaxE_;

    std::vector<float> cl3d_C2dEdensity_L05_;
    std::vector<float> cl3d_C2dEdensity_L10_;
    std::vector<float> cl3d_C2dEdensity_L15_;
    std::vector<float> cl3d_C2dEdensity_L30_;    
    std::vector<float> cl3d_C2dEdensity_L35_;

    std::vector<float> cl3d_C2dMipTdensity_L05_;
    std::vector<float> cl3d_C2dMipTdensity_L10_;
    std::vector<float> cl3d_C2dMipTdensity_L15_;
    std::vector<float> cl3d_C2dMipTdensity_L30_;    
    std::vector<float> cl3d_C2dMipTdensity_L35_;

    std::vector< std::vector<float> > cl3d_C2dEdensities_;
    std::vector< std::vector<float> > cl3d_C2dMipTdensities_;

    std::vector<int> cl3d_firstC2dLayer_;
    std::vector<int> cl3d_lastC2dLayer_;
    std::vector<int> cl3ds_NxEndcap_;
    std::vector< std::vector<float> > cl2d_See_;
    std::vector<float> cl3d_halfDepth_;
    std::vector<float> cl3d_pTxHalfDepth_;
};

DEFINE_EDM_PLUGIN(HGCalTriggerNtupleFactory,
    HGCalTriggerNtupleHGCMulticlusters,
    "HGCalTriggerNtupleHGCMulticlusters" );


HGCalTriggerNtupleHGCMulticlusters::
HGCalTriggerNtupleHGCMulticlusters(const edm::ParameterSet& conf):HGCalTriggerNtupleBase(conf)
{
}

void
HGCalTriggerNtupleHGCMulticlusters::
initialize(TTree& tree, const edm::ParameterSet& conf, edm::ConsumesCollector&& collector)
{
  multiclusters_token_ = collector.consumes<l1t::HGCalMulticlusterBxCollection>(conf.getParameter<edm::InputTag>("Multiclusters"));

  tree.Branch("cl3d_n", &cl3d_n_, "cl3d_n/I");
  tree.Branch("cl3d_pt", &cl3d_pt_);
  tree.Branch("cl3d_energy", &cl3d_energy_);
  tree.Branch("cl3d_eta", &cl3d_eta_);
  tree.Branch("cl3d_phi", &cl3d_phi_);
  tree.Branch("cl3d_cl_mipPt", &cl3d_cl_mipPt_);
  tree.Branch("cl3d_cl_layer", &cl3d_cl_layer_);
  tree.Branch("cl3d_x", &cl3d_x_);
  tree.Branch("cl3d_y", &cl3d_y_);
  tree.Branch("cl3d_z", &cl3d_z_);
  tree.Branch("cl3d_xN", &cl3d_xN_);
  tree.Branch("cl3d_yN", &cl3d_yN_);
  tree.Branch("cl3d_nclu", &cl3d_nclu_);
  tree.Branch("cl3d_firstC2dLayer", &cl3d_firstC2dLayer_);
  tree.Branch("cl3d_lastC2dLayer", &cl3d_lastC2dLayer_);
  tree.Branch("cl3ds_NxEndcap", &cl3ds_NxEndcap_ );
  tree.Branch("cl2d_See", &cl2d_See_);
  tree.Branch("cl3d_halfDepth", &cl3d_halfDepth_);
  tree.Branch("cl3d_pTxHalfDepth", &cl3d_pTxHalfDepth_);
  tree.Branch("cl3d_C2dLayer1stMaxE", &cl3d_C2dLayer1stMaxE_);
  tree.Branch("cl3d_C2dLayer2ndMaxE", &cl3d_C2dLayer2ndMaxE_);
  tree.Branch("cl3d_C2dEdensities", &cl3d_C2dEdensities_);
  tree.Branch("cl3d_C2dEdensity_L05", &cl3d_C2dEdensity_L05_);
  tree.Branch("cl3d_C2dEdensity_L10", &cl3d_C2dEdensity_L10_);
  tree.Branch("cl3d_C2dEdensity_L15", &cl3d_C2dEdensity_L15_);
  tree.Branch("cl3d_C2dEdensity_L30", &cl3d_C2dEdensity_L30_);
  tree.Branch("cl3d_C2dEdensity_L35", &cl3d_C2dEdensity_L35_);
  tree.Branch("cl3d_C2dMipTdensities", &cl3d_C2dMipTdensities_);  
  tree.Branch("cl3d_C2dMipTdensity_L05", &cl3d_C2dMipTdensity_L05_);
  tree.Branch("cl3d_C2dMipTdensity_L10", &cl3d_C2dMipTdensity_L10_);
  tree.Branch("cl3d_C2dMipTdensity_L15", &cl3d_C2dMipTdensity_L15_);
  tree.Branch("cl3d_C2dMipTdensity_L30", &cl3d_C2dMipTdensity_L30_);
  tree.Branch("cl3d_C2dMipTdensity_L35", &cl3d_C2dMipTdensity_L35_);

}

void
HGCalTriggerNtupleHGCMulticlusters::
fill(const edm::Event& e, const edm::EventSetup& es)
{

  // retrieve clusters 3D
  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters_h;
  e.getByToken(multiclusters_token_, multiclusters_h);
  const l1t::HGCalMulticlusterBxCollection& multiclusters = *multiclusters_h;

  // retrieve geometry
  edm::ESHandle<HGCalTriggerGeometryBase> geometry;
  es.get<IdealGeometryRecord>().get(geometry);

  clear();
  for(auto cl3d_itr=multiclusters.begin(0); cl3d_itr!=multiclusters.end(0); cl3d_itr++)
  {
    cl3d_n_++;
    // physical values 
    cl3d_pt_.emplace_back(cl3d_itr->pt());
    cl3d_energy_.emplace_back(cl3d_itr->energy());
    cl3d_eta_.emplace_back(cl3d_itr->eta());
    cl3d_phi_.emplace_back(cl3d_itr->phi());
    cl3d_x_.emplace_back(cl3d_itr->centre().x());
    cl3d_y_.emplace_back(cl3d_itr->centre().y());
    cl3d_z_.emplace_back(cl3d_itr->centre().z());
    cl3d_xN_.emplace_back(cl3d_itr->centreProj().x());
    cl3d_yN_.emplace_back(cl3d_itr->centreProj().y());

    /* compute the distance between the beginning and the barycentre of a C3d */
    edm::PtrVector<l1t::HGCalCluster>::const_iterator first_clu = cl3d_itr->constituents_begin();
    float halfDepth = sqrt( ( (*first_clu)->centre().x() - cl3d_itr->centre().x() ) * ( (*first_clu)->centre().x() - cl3d_itr->centre().x() )
                            + ( (*first_clu)->centre().y() - cl3d_itr->centre().y() ) * ( (*first_clu)->centre().y() - cl3d_itr->centre().y() )
                            + ( (*first_clu)->centre().z() - cl3d_itr->centre().z() ) * ( (*first_clu)->centre().z() - cl3d_itr->centre().z() ) );
    float pTxHalfDepth = cl3d_itr->pt() *  halfDepth;

    cl3d_halfDepth_.emplace_back(halfDepth);
    cl3d_pTxHalfDepth_.emplace_back(pTxHalfDepth);

    /* loop over the 2D-clusters */    
    const edm::PtrVector<l1t::HGCalCluster> pertinentClu = cl3d_itr->constituents();
    cl3d_nclu_.emplace_back(cl3d_itr->constituents().size());

    std::vector<int> layerC2dSet;
    std::vector<float> mipPtC2dSet;
    std::vector<float> thiscl2d_See_;
    std::vector<float > EdensitySet;
    std::vector<float > MipTdensitySet;

    float largest = -1.;
    float second = -1.;

    for( edm::PtrVector<l1t::HGCalCluster>::const_iterator it_clu=pertinentClu.begin(); it_clu<pertinentClu.end(); it_clu++){
        
       /* loop over the trigger cells pertinent to a cluster */    
        const edm::PtrVector<l1t::HGCalTriggerCell> triggerCells = (*it_clu)->constituents();
        unsigned int Ntc = triggerCells.size();

        l1t::ClusterShapes shape;
        for(unsigned int itc=0; itc<Ntc;itc++){
            l1t::HGCalTriggerCell thistc = *triggerCells[itc];
            shape.Add(thistc.energy(),thistc.eta(),thistc.phi(),thistc.position().z());
        }
        
        if( (*it_clu)->energy() > largest) {
            second = largest;
            largest = (*it_clu)->energy();
        } else if ( (*it_clu)->energy() > second) {
            second = (*it_clu)->energy();
        }

        int layerN = -1; 
        if( (*it_clu)->subdetId()==3 ){
            layerN = (*it_clu)->layer();
        }
        else if((*it_clu)->subdetId()==4){
            layerN = (*it_clu)->layer()+28;
        }                    
        
        layerC2dSet.push_back(layerN);
        mipPtC2dSet.push_back( (*it_clu)->mipPt() );
        thiscl2d_See_.emplace_back(shape.SigmaEtaEta());

        EdensitySet.emplace_back( (*it_clu)->energy()/Ntc );
        MipTdensitySet.emplace_back( (*it_clu)->mipPt()/Ntc );

        if(layerN==5)        { 
            cl3d_C2dEdensity_L05_.emplace_back( (*it_clu)->energy()/Ntc ); 
            cl3d_C2dMipTdensity_L05_.emplace_back( (*it_clu)->mipPt()/Ntc ); 
        }
        if(layerN==10){ 
            cl3d_C2dEdensity_L10_.emplace_back( (*it_clu)->energy()/Ntc );
            cl3d_C2dMipTdensity_L10_.emplace_back( (*it_clu)->mipPt()/Ntc ); 
        }
        if(layerN==15){ 
            cl3d_C2dEdensity_L15_.emplace_back( (*it_clu)->energy()/Ntc ); 
            cl3d_C2dMipTdensity_L15_.emplace_back( (*it_clu)->mipPt()/Ntc ); 
        }
        if(layerN==30){ 
            cl3d_C2dEdensity_L30_.emplace_back( (*it_clu)->energy()/Ntc ); 
            cl3d_C2dMipTdensity_L30_.emplace_back( (*it_clu)->mipPt()/Ntc ); 
        }
        if(layerN==35){ 
            cl3d_C2dEdensity_L35_.emplace_back( (*it_clu)->energy()/Ntc ); 
            cl3d_C2dMipTdensity_L35_.emplace_back( (*it_clu)->mipPt()/Ntc ); 
        }

    }

    cl3d_C2dLayer1stMaxE_.emplace_back( largest );
    cl3d_C2dLayer2ndMaxE_.emplace_back( second );

    cl3d_firstC2dLayer_.emplace_back( layerC2dSet.front() );
    cl3d_lastC2dLayer_.emplace_back( layerC2dSet.back() );
    cl3d_C2dEdensities_.emplace_back( EdensitySet );
    cl3d_C2dMipTdensities_.emplace_back( MipTdensitySet );

    cl3d_cl_mipPt_.emplace_back(mipPtC2dSet);
    cl3d_cl_layer_.emplace_back(layerC2dSet);
    cl2d_See_.emplace_back(thiscl2d_See_);


    thiscl2d_See_.clear();

  }

  int zsides[2] = {-1, 1};
  for(int i_z(0); i_z<2; ++i_z){
      int cl3dNxEndcap = 0;
      for(auto cl3d_itr=multiclusters.begin(0); cl3d_itr!=multiclusters.end(0); cl3d_itr++){             
          if( cl3d_itr->zside() == zsides[i_z] ){              
              cl3dNxEndcap++;
          }              
      }
      /* in the event there are clNxlayer */
      
      cl3ds_NxEndcap_.emplace_back(cl3dNxEndcap);
  }
}


void
HGCalTriggerNtupleHGCMulticlusters::
clear()
{
  cl3d_n_ = 0;
  cl3d_pt_.clear();
  cl3d_energy_.clear();
  cl3d_eta_.clear();
  cl3d_phi_.clear();
  cl3d_cl_mipPt_.clear();
  cl3d_cl_layer_.clear();
  cl3d_x_.clear();
  cl3d_y_.clear();
  cl3d_z_.clear();
  cl3d_xN_.clear();
  cl3d_yN_.clear();
  cl3d_nclu_.clear();
  cl3d_firstC2dLayer_.clear();
  cl3d_lastC2dLayer_.clear();
  cl3ds_NxEndcap_.clear();
  cl2d_See_.clear();
  cl3d_halfDepth_.clear();
  cl3d_pTxHalfDepth_.clear();
  cl3d_C2dLayer1stMaxE_.clear();
  cl3d_C2dLayer2ndMaxE_.clear();
  cl3d_C2dEdensities_.clear();
  cl3d_C2dEdensity_L05_.clear();
  cl3d_C2dEdensity_L10_.clear();
  cl3d_C2dEdensity_L15_.clear();
  cl3d_C2dEdensity_L30_.clear();
  cl3d_C2dEdensity_L35_.clear();
  cl3d_C2dMipTdensities_.clear();
  cl3d_C2dMipTdensity_L05_.clear();
  cl3d_C2dMipTdensity_L10_.clear();
  cl3d_C2dMipTdensity_L15_.clear();
  cl3d_C2dMipTdensity_L30_.clear();
  cl3d_C2dMipTdensity_L35_.clear();


}




