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

// SFS todo: update description
/// \brief extract relevant mc truth information that allows to compute efficiency, purity and more quantities for the photon conversion analysis.
/// dependencies: none
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "gammaTables.h"

#include "TVector3.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using tracksAndMcLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;

struct gammaConversionsMcTruthOnlyConsumer {

  Configurable<float> fEtaMax{"fEtaMax", 0.8, "aMaximum photon eta"};

  HistogramRegistry registry{
    "registry",
    {
      {"hEtaDiff", "hEtaDiff", {HistType::kTH1F, {{400, -2.f, 2.f}}}},
      
      {"hNDaughters", "hNDaughters", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
      {"hPdgCodeDaughters", "hPdgCodeDaughters", {HistType::kTH1F, {{2000, -1000.f, 1000.f}}}},
      {"hNElectrons", "hNElectrons", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
      
      {"hGammaProdInEtaAccP", "hGammaProdInEtaAccP", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      
      {"hGammaConvertedR", "hGammaConvertedR", {HistType::kTH1F, {{1600, 0.f, 500.f}}}},
      {"hGammaConvertedRselP", "hGammaConvertedRselP", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      
      {"hGammaConvertedEtaP", "hGammaConvertedEtaP", {HistType::kTH2F, {{400, -2.f, 2.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedEtaR", "hGammaConvertedEtaR", {HistType::kTH2F, {{400, -2.f, 2.f}, {400, 0.f, 250.f}}}},
      {"hGammaConvertedEtaZ", "hGammaConvertedEtaZ", {HistType::kTH2F, {{400, -2.f, 2.f}, {400, -250.f, 250.f}}}},
      
      {"hGammaConvertedRP", "hGammaConvertedRP", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedRZ", "hGammaConvertedRZ", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, -250.f, 250.f}}}},

      {"hGammaConvertedZP", "hGammaConvertedZP", {HistType::kTH2F, {{400, -250.f, 250.f}, {400, 0.f, 25.f}}}},
      
      {"hGammaProdInEtaAccPt", "hGammaProdInEtaAccPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaMoreThanTwoDaughtersPt", "hGammaMoreThanTwoDaughtersPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      
      {"hGammaConvertedRPt", "hGammaConvertedRPt", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedRselPt", "hGammaConvertedRselPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hPeculiarOccurences", "hPeculiarOccurences", {HistType::kTH1F, {{50, -25.f, 25.f}}}},
    },
  };

  // loop over MC truth McCollisions
  void process(aod::McCollision const &theMcCollision,
               aod::MCGammas const  &theMcGammas,
               aod::MCGammaDaughters const &theMcGammaDaughters)
  {
    for (auto &lMcGamma : theMcGammas) {

      if (std::abs(lMcGamma.eta()) < fEtaMax) { // SFS todo: track v0 eta??

        registry.fill(HIST("hGammaProdInEtaAccP"), lMcGamma.p());
        registry.fill(HIST("hGammaProdInEtaAccPt"), lMcGamma.pt());
        
        // todo: look at theMcGammaDaughters, slice, check number of daughters, that they are electrons,
        // etc 
        int const lNDaughters = lMcGamma.nDaughters();
        if (lNDaughters){
          auto lDaughters = theMcGammaDaughters.sliceBy(aod::truthOnly2::motherId, lMcGamma.gamma());
          
          if (lDaughters.size() != lNDaughters){
            LOGF(warning, "SFS differing number of daughters. This should never happen. %d vs %d", lDaughters.size(), lNDaughters);
            registry.fill(HIST("hPeculiarOccurences"), -0.5);
          }
          
          size_t lNElectrons = 0;
          for (auto &lDaughter : lDaughters) {
            registry.fill(HIST("hPdgCodeDaughters"), 0.5 + lDaughter.pdgCode());
            lNElectrons += std::abs(lDaughter.pdgCode()) == 11;
          }
          registry.fill(HIST("hNElectrons"), 0.5 + lNElectrons);
          
          // "regular" conversion
          //~ if (lNElectrons == 2) {
          if (lNElectrons >= 2) {
            
            if (lNDaughters != 2) {
              registry.fill(HIST("hPeculiarOccurences"), 0.5);
            }
            
            // access first daughter to get conversion point
            auto const &lDaughter0 = lDaughters.begin();
            float lConversionRadius = std::sqrt(std::pow(lDaughter0.vx(), 2) + std::pow(lDaughter0.vy(), 2));
            registry.fill(HIST("hGammaConvertedEtaP"), lMcGamma.eta(), lMcGamma.p());
            registry.fill(HIST("hGammaConvertedEtaR"), lMcGamma.eta(), lConversionRadius);
            registry.fill(HIST("hGammaConvertedEtaZ"), lMcGamma.eta(), lDaughter0.vz());
            registry.fill(HIST("hGammaConvertedR"), lConversionRadius);
            registry.fill(HIST("hGammaConvertedRP"), lConversionRadius, lMcGamma.p());
            registry.fill(HIST("hGammaConvertedRZ"), lConversionRadius, lDaughter0.vz());
            registry.fill(HIST("hGammaConvertedZP"), lDaughter0.vz(), lMcGamma.p());

            registry.fill(HIST("hGammaConvertedRPt"), lConversionRadius, lMcGamma.pt());
            
            TVector3 lDaughter0Vtx(lDaughter0.vx(),lDaughter0.vy(), lDaughter0.vz());
            
            float_t lEtaDiff = lDaughter0Vtx.Eta() - lMcGamma.eta();
            registry.fill(HIST("hEtaDiff"), lEtaDiff);
            
            

            if (lConversionRadius > 5. && lConversionRadius < 180.) {
              registry.fill(HIST("hGammaConvertedRselP"), lMcGamma.p());
              registry.fill(HIST("hGammaConvertedRselPt"), lMcGamma.pt());
            }
          }
          if (lNElectrons > 2) {
            registry.fill(HIST("hGammaMoreThanTwoDaughtersPt"), lMcGamma.pt());
          }
        }
        registry.fill(HIST("hNDaughters"), 0.5 + lNDaughters);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<gammaConversionsMcTruthOnlyConsumer>(cfgc)};
}
