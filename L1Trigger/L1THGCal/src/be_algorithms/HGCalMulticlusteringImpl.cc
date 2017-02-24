#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalMulticlusteringImpl.h"



HGCalMulticlusteringImpl::HGCalMulticlusteringImpl(const edm::ParameterSet& beCodecConfig){
    dR_forC3d_ = beCodecConfig.getParameter<double>("dR_searchNeighbour");
}

void HGCalMulticlusteringImpl::clusterise( std::unique_ptr<l1t::HGCalClusterBxCollection>& clusters_, 
                                           std::unique_ptr<l1t::HGCalMulticlusterBxCollection>& multiclusters_
    ){
   
    std::cout << "------ Max-Distance in the normalized plane to search for next-layer C2d to merge into the C3d: " << dR_forC3d_ << std::endl;  

    if(clusters_->size()>0){
        std::vector<size_t> isMerged;

        size_t seedx=0;
 
        for(l1t::HGCalClusterBxCollection::const_iterator cl = clusters_->begin(); cl != clusters_->end(); ++cl, ++seedx){
//            edm::PtrVector<l1t::HGCalCluster> ClusterCollection;
        
            l1t::HGCalMulticluster multicluster( reco::LeafCandidate::LorentzVector(), 0, 0, 0);//, ClusterCollection);
            double_t tmpEta = 0.;
            double_t tmpPhi = 0.;           
            double_t C3d_pt  = 0.;
            double_t C3d_eta = 0.;
            double_t C3d_phi = 0.;
            uint32_t C3d_hwPtEm = 0;
            uint32_t C3d_hwPtHad = 0;
            uint32_t totLayer = 0;

            bool skip=false;
        
            //std::cout << "In the Cl collection, seed the C3d with this : " << seedx << " - "<< cl->p4().Pt() << " eta: " <<  cl->p4().Eta() << " --> layer " << cl->layer() << "  skip before 2nd loop "<< skip << std::endl;                
        
            size_t idx=0;
            for(l1t::HGCalClusterBxCollection::const_iterator cl_aux = clusters_->begin(); cl_aux != clusters_->end(); ++cl_aux, ++idx){
                //  std::cout << "     loop over Cl again and search for match:" << "   idx: " << idx << "  eta: " << cl_aux->p4().Eta() << std::endl;
                //std::cout << "   before isMerged loop: " << skip<< std::endl;
                for(size_t i(0); i<isMerged.size(); i++){
                    //std::cout <<  isMerged.at(i) << ", ";
                    if(idx==isMerged.at(i)){
                        skip=true;
                        continue;
                    }
                }
                //std::cout << "\n";
                double dR =  deltaR( cl->p4().Eta(), cl_aux->p4().Eta(), cl->p4().Phi(), cl_aux->p4().Phi() ); 
                std::cout << "looping on the cl directly from the collection : "<< cl->p4().pt() << "  --> layer " << cl->layer() << " dR: " << dR << "  SKIP var = " << skip << std::endl;
            
                if(skip){
                    skip=false;
                    //std::cout << "     the cl considered has been already merged!!";
                    continue;
                }
                if( dR < dR_forC3d_*1 ){
                    //    std::cout << "     The idx "<< idx << " Cl has been matched and kept for 3D to the " << seedx 
                    //          << " - "<< cl_aux->p4().Pt() << " eta: " <<  cl_aux->p4().Eta() 
                    //          << " --> layer " << cl_aux->layer() << std::endl;             
                    isMerged.push_back(idx);
                    tmpEta+=cl_aux->p4().Eta() * cl_aux->p4().Pt();
                    tmpPhi+=cl_aux->p4().Phi() * cl_aux->p4().Pt();
                    C3d_pt+=cl_aux->p4().Pt();
                    C3d_hwPtEm+=cl_aux->hwPtEm();
                    C3d_hwPtHad+=cl_aux->hwPtHad();
                    totLayer++;
                    //ClusterCollection.push_back(&(*cl));
                }
            }
        
            std::cout <<"STO PER ENTRARE NEL MULTICLUSTERING" << std::endl;
            if( totLayer > 2){
                // edm::PtrVector<l1t::HGCalCluster> ClusterCollection; //push_back()
                multicluster.setNtotLayer(totLayer);
                multicluster.setHwPtEm(C3d_hwPtEm);
                multicluster.setHwPtHad(C3d_hwPtHad);
                C3d_eta=tmpEta/C3d_pt;
                C3d_phi=tmpPhi/C3d_pt;                
                math::PtEtaPhiMLorentzVector calib3dP4(C3d_pt, C3d_eta, C3d_phi, 0 );
                multicluster.setP4(calib3dP4);                    
                std::cout << "  A MULTICLUSTER has been built with pt, eta, phi = " << C3d_pt << ", " << C3d_eta << ", "<< C3d_phi <<  std::endl;
                multiclusters_->push_back(0,multicluster);
            }                    
        }
    }
} 
