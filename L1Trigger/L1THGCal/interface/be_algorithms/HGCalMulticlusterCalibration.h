#ifndef __L1Trigger_L1THGCal_HGCalMulticlusterCalibration_h__
#define __L1Trigger_L1THGCal_HGCalMulticlusterCalibration_h__

#include <stdint.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"

class HGCalMulticlusterCalibration{

public:
  
    HGCalMulticlusterCalibration(const edm::ParameterSet &conf);    
    void calibrateC3d(l1t::HGCalMulticluster&);
    void matrixInversionCalibC3d(l1t::HGCalMulticluster&);

private:
    
    double globalCalib_;

};

#endif
