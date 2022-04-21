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

/// \brief extract relevant mc truth information that allows to compute efficiency, purity and more quantities for the photon conversion analysis.
/// dependencies: none
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using tracksAndMcLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;

struct gammaConversionsMcTruthOnlyConsumer {

  HistogramRegistry registry{
    "registry",
    {
      {"hNDaughters", "hNDaughters", {HistType::kTH1F, {{50, 0.f, 50.f}}}},
      {"hCollisionZ", "hCollisionZ", {HistType::kTH1F, {{800, -10.f, 10.f}}}},
      {"hGammaProdInEtaAccP", "hGammaProdInEtaAccP", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaProdInEtaAccPt", "hGammaProdInEtaAccPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaMoreThanTwoDaughtersPt", "hGammaMoreThanTwoDaughtersPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      {"hGammaConvertedRP", "hGammaConvertedRP", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedR", "hGammaConvertedR", {HistType::kTH1F, {{1600, 0.f, 500.f}}}},
      {"hGammaConvertedRselP", "hGammaConvertedRselP", {HistType::kTH1F, {{800, 0.f, 25.f}}}},
      
      {"hGammaConvertedRPt", "hGammaConvertedRPt", {HistType::kTH2F, {{400, 0.f, 250.f}, {400, 0.f, 25.f}}}},
      {"hGammaConvertedRselPt", "hGammaConvertedRselPt", {HistType::kTH1F, {{800, 0.f, 25.f}}}}
      },
  };

  // loop over MC truth McCollisions
  void process(aod::McCollision const &theMcCollision,
               soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                          aod::Collisions>> const &theCollisions,
               aod::MCGammas const  &theMcGammas,
               aod::MCGammaTracks const &theMcGammaTracks)
  {
    for (auto &lCollision : theCollisions) {

      registry.fill(HIST("hCollisionZ"), lCollision.posZ());
      for (auto &lMcParticle : theMcParticles) {

        if ((lMcParticle.pdgCode() == 22) &&
            (std::abs(lMcParticle.eta()) < 0.8)) { // SFS todo: track v0 eta??

          registry.fill(HIST("hGammaProdInEtaAccP"), lMcParticle.p());
          registry.fill(HIST("hGammaProdInEtaAccPt"), lMcParticle.pt());
          
          size_t lNDaughters = 0;
          size_t lBothElectrons = 0;
          if (lMcParticle.has_daughters()) {
            for (auto &lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
              ++lNDaughters;
              lBothElectrons += std::abs(lDaughter.pdgCode()) == 11;
            }
            registry.fill(HIST("hNDaughters"), 0.5 + lNDaughters);
            if (lBothElectrons == 2) {
              if (lNDaughters != 2){
                LOGF(info, "SFS ALARM");
              }
              for (auto &lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
                float lConversionRadius = std::sqrt(std::pow(lDaughter.vx(), 2) + std::pow(lDaughter.vy(), 2));
                registry.fill(HIST("hGammaConvertedR"), lConversionRadius);
                registry.fill(HIST("hGammaConvertedRP"), lConversionRadius, lMcParticle.p());
                registry.fill(HIST("hGammaConvertedRPt"), lConversionRadius, lMcParticle.pt());

                if (lConversionRadius > 5. && lConversionRadius < 180.) {
                  registry.fill(HIST("hGammaConvertedRselP"), lMcParticle.p());
                  registry.fill(HIST("hGammaConvertedRselPt"), lMcParticle.pt());
                }
                break;
              }
            }
            if (lBothElectrons > 2) {
              registry.fill(HIST("hGammaMoreThanTwoDaughtersPt"), lMcParticle.pt());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<gammaConversionsMcTruthOnlyConsumer>(cfgc)};
}
