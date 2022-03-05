// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief apply cuts to get photons from table aod::StoredV0Datas produced by o2-analysis-lf-lambdakzerobuilder and produce plots.
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/StrangenessTables.h"


#include <TH1.h>
#include <TH1F.h>

#include <cmath>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct gammaConversionsMcTruthOnly {
  OutputObj<TH1F> hEventCounter{TH1F("hEventCounter", "hEventCounter", 2, 0.f, 2.f)};
  OutputObj<TH1F> hGammaProdInAcc{TH1F("hGammaProdInAcc", "hGammaProdInAcc", 800, 0.f, 25.f)};
  OutputObj<TH1F> hGammaMoreThanTwoDaughtersPt{TH1F("hGammaMoreThanTwoDaughtersPt", "hGammaMoreThanTwoDaughtersPt", 800, 0.f, 25.f)};
  OutputObj<TH1F> hGammaConvertedPt{TH1F("hGammaConvertedPt", "hGammaConvertedPt", 800, 0.f, 25.f)};
  OutputObj<TH1F> hGammaConvertedR{TH1F("hGammaConvertedR", "hGammaConvertedR", 1600, 0.f, 500.f)};
  OutputObj<TH2F> hGammaConvertedRPt{TH2F("hGammaConvertedRPt", "hGammaConvertedRPt", 400, 0.f, 250.f, 400, 0.f, 25.f)};
  
  // loop over MC truth McCollisions
  void process(aod::McCollision                             const& mcCollision, 
               soa::SmallGroups<soa::Join<aod::McCollisionLabels, 
                                          aod::Collisions>> const& collisions,
               aod::McParticles                             const& theMcParticles)
  {    
    hEventCounter->Fill(0.5);
    for (auto& collision : collisions) {
      hEventCounter->Fill(1.5);
      for (auto& lMcParticle : theMcParticles) {  
        
        if ((lMcParticle.pdgCode() == 22) && 
            (std::abs(lMcParticle.eta()) < 0.8)){ // SFS todo: track v0 eta??
          
          hGammaProdInAcc->Fill(lMcParticle.pt());
          
          size_t lBothElectrons = 0;
          if (lMcParticle.has_daughters()) {
            for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
              lBothElectrons += std::abs(lDaughter.pdgCode()) == 11;
            }
            if (lBothElectrons == 2){
              for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
                float lConversionRadius = std::sqrt(std::pow(lDaughter.vx(),2) + std::pow(lDaughter.vy(),2));
                hGammaConvertedR->Fill(lConversionRadius);
                hGammaConvertedRPt->Fill(lConversionRadius, lMcParticle.pt());
                
                if (lConversionRadius > 5. && lConversionRadius < 180.){
                  hGammaConvertedPt->Fill(lMcParticle.pt());  
                }
                break;
              }
            }
            if (lBothElectrons > 2){
              hGammaMoreThanTwoDaughtersPt->Fill(lMcParticle.pt());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<gammaConversionsMcTruthOnly>(cfgc)};
}
