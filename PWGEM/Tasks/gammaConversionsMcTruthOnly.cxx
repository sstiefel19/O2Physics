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

#include <TH1.h>
#include <TH1F.h>

#include <cmath>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"


// Loop over MCColisions and get corresponding collisions (there can be more than one)
// For each of them get the corresponding tracks
// Note the use of "SmallGroups" template, that allows to handle both Run 2, where
// we have exactly 1-to-1 correspondence between collisions and mc collisions, and
// Run 3, where we can have 0, 1, or more collisions for a given mc collision
struct LoopOverMcMatched {
  OutputObj<TH1F> etaDiff{TH1F("etaDiff", ";eta_{MC} - eta_{Rec}", 100, -2, 2)};
  void process(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions,
               soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles const& mcParticles)
  {
    // access MC truth information with mcCollision() and mcParticle() methods
    LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
    for (auto& collision : collisions) {
      LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());

      // NOTE this will be replaced by a improved grouping in the future
      auto groupedTracks = tracks.sliceBy(aod::track::collisionId, collision.globalIndex());
      LOGF(info, "  which has %d tracks", groupedTracks.size());
      for (auto& track : groupedTracks) {
        if (!track.has_mcParticle()) {
          LOGF(warning, "No MC particle for track, skip...");
          continue;
        }
        etaDiff->Fill(track.mcParticle().eta() - track.eta());
      }
    }
  }
};

// Simple access to collision
struct gammaConversionsMcTruthOnly {
  OutputObj<TH1F> hGammaPt{TH1F("hGammaPt", "hGammaPt", 800, 0.f, 25.f)};
  OutputObj<TH1F> hGammaConvertedPt{TH1F("hGammaConvertedPt", "hGammaConvertedPt", 800, 0.f, 25.f)};
  OutputObj<TH1F> hGammaConvertedR{TH1F("hGammaConvertedR", "hGammaConvertedR", 800, 0.f, 250.f)};


  // loop over MC truth McCollisions
  void process(aod::McCollision const& theMcCollision,
               aod::McParticles const& theMcParticles)
  {    
    for (auto& lMcParticle : theMcParticles) {
      
      if (lMcParticle.pdgCode() == 22){
        hGammaPt->Fill(lMcParticle.pt());
        
        size_t lBothElectrons = 0;
        if (lMcParticle.has_daughters()) {
          for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
            lBothElectrons += std::abs(lDaughter.pdgCode()) == 11;
          }
          if (lBothElectrons == 2){
            hGammaConvertedPt->Fill(lMcParticle.pt());
            
            // SFS todo: this counts every conversion twice. correct
            for (auto& lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
              hGammaConvertedR->Fill(std::sqrt(std::pow(lDaughter.vx(),2) + std::pow(lDaughter.vy(),2)));
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
