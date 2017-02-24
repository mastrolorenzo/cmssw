
#ifndef __L1Trigger_L1THGCal_HGCalClusteringImpl_h__
#define __L1Trigger_L1THGCal_HGCalClusteringImpl_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalTriggerCellCalibration.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"


class HGCalClusteringImpl{

public:
  
    HGCalClusteringImpl( const edm::ParameterSet &conf);    
    void clusterise( std::unique_ptr<l1t::HGCalTriggerCellBxCollection>& trgcells_, 
                     std::unique_ptr<l1t::HGCalClusterBxCollection> & clusters_);

private:
    
    double seedThr_;
    double tcThr_;

};

#endif
