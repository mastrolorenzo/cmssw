#include "L1Trigger/L1THGCal/interface/be_algorithms/HGCalMulticlusterCalibration.h"

//class constructor
HGCalMulticlusterCalibration::HGCalMulticlusterCalibration(const edm::ParameterSet& beCodecConfig){    
    globalCalib_ = beCodecConfig.getParameter<double>("calibSF_multicluster");
    edm::LogInfo("HGCalMulticlusterCalibration") << "Multicluster global calibration factor: " << globalCalib_;  
}

void HGCalMulticlusterCalibration::calibrateC3d(l1t::HGCalMulticluster& multicluster)
{

    /* calibrate the momentum of the multicluster */
    math::PtEtaPhiMLorentzVector calibP4(  multicluster.pt() * globalCalib_, 
                                           multicluster.eta(), 
                                           multicluster.phi(), 
                                           0. );
    /* overwriting the 4p with the calibrated 4p */     
    multicluster.setP4( calibP4 );

} 

void HGCalMulticlusterCalibration::matrixInversionCalibC3d(l1t::HGCalMulticluster& multicluster)
{
    /* Apply the calibration coefficient estimated with the Matrix Inversion Method */
} 


