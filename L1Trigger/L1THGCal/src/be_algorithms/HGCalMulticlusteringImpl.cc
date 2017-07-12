#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalMulticlusteringImpl.h"
#include "DataFormats/Math/interface/deltaR.h"


HGCalMulticlusteringImpl::HGCalMulticlusteringImpl( const edm::ParameterSet& conf ) :
    dr_(conf.getParameter<double>("dR_multicluster")),
    ptC3dThreshold_(conf.getParameter<double>("minPt_multicluster")),
    calibSF_(conf.getParameter<double>("calibSF_multicluster")),
    dEdX_(conf.getParameter< std::vector<double> >("calibCoeffMtx"))
{    
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster dR for Near Neighbour search: " << dr_;  
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster minimum transverse-momentum: " << ptC3dThreshold_;
    edm::LogInfo("HGCalMulticlusterParameters") << "Multicluster global calibration factor: " << calibSF_;
    for(unsigned i(0); i<dEdX_.size(); i++){
        std::cout << "i " << i << "   a_i = " << dEdX_.at(i) << std::endl; 
    }

}


bool HGCalMulticlusteringImpl::isPertinent( const l1t::HGCalCluster & clu, 
                                            const l1t::HGCalMulticluster & mclu, 
                                            double dR ) const
{
    HGCalDetId cluDetId( clu.detId() );
    HGCalDetId firstClusterDetId( mclu.detId() );
    
    if( cluDetId.zside() != firstClusterDetId.zside() ){
        return false;
    }
    if( ( mclu.centreProj() - clu.centreProj() ).mag() < dR ){
        return true;
    }
    return false;

}


void HGCalMulticlusteringImpl::clusterize( const edm::PtrVector<l1t::HGCalCluster> & clustersPtrs, 
                                           l1t::HGCalMulticlusterBxCollection & multiclusters)
{
           
    std::vector<l1t::HGCalMulticluster> multiclustersTmp;

    int iclu = 0;
    for(edm::PtrVector<l1t::HGCalCluster>::const_iterator clu = clustersPtrs.begin(); clu != clustersPtrs.end(); ++clu, ++iclu){
        
        int imclu=0;
        vector<int> tcPertinentMulticlusters;
        for( const auto& mclu : multiclustersTmp ){
            if( this->isPertinent(**clu, mclu, dr_) ){
                tcPertinentMulticlusters.push_back(imclu);
            }
            ++imclu;
        }
        if( tcPertinentMulticlusters.size() == 0 ){
            multiclustersTmp.emplace_back( *clu );
        }
        else{
            unsigned minDist = 1;
            unsigned targetMulticlu = 0; 
            for( int imclu : tcPertinentMulticlusters ){
                double d = ( multiclustersTmp.at(imclu).centreProj() - (*clu)->centreProj() ).mag() ;
                if( d < minDist ){
                    minDist = d;
                    targetMulticlu = imclu;
                }
            } 

            multiclustersTmp.at( targetMulticlu ).addConstituent( *clu );
            
        }        
    }

    /* making the collection of multiclusters */
    for( unsigned i(0); i<multiclustersTmp.size(); ++i ){

        double calibPt=0.;
        bool MatrixCalib = true;
        if(MatrixCalib){
            const edm::PtrVector<l1t::HGCalCluster> pertinentClu = multiclustersTmp.at(i).constituents();
            for( edm::PtrVector<l1t::HGCalCluster>::const_iterator it_clu=pertinentClu.begin(); it_clu<pertinentClu.end(); it_clu++){
                int layerN = -1; 
                if( (*it_clu)->subdetId()==3 ){
                    layerN = (*it_clu)->layer();
                }
                else if((*it_clu)->subdetId()==4){
                    layerN = (*it_clu)->layer()+28;
                }
//            std::cout << "layer " << layerN << "   a_i = " << dEdX_.at(layerN) << " e-c2d = " << (*it_clu)->pt() << std::endl; 
                calibPt += dEdX_.at(layerN) * (*it_clu)->mipPt();
            }
            
//        std::cout << "SF = "<<  multiclustersTmp.at(i).pt() * calibSF_ << "  matrix-inversion = "<< calibPt << std::endl;
     
        }else if(!MatrixCalib){        
            calibPt = multiclustersTmp.at(i).pt() * calibSF_; 
        }
        math::PtEtaPhiMLorentzVector calibP4( calibPt, 
                                              multiclustersTmp.at(i).eta(), 
                                              multiclustersTmp.at(i).phi(), 
                                              0. );
        
        // overwriting the 4p with the calibrated 4p     
        multiclustersTmp.at(i).setP4( calibP4 );
        if( multiclustersTmp.at(i).pt() > ptC3dThreshold_ ){
            multiclusters.push_back( 0, multiclustersTmp.at(i));  
        }
    }
    
}
