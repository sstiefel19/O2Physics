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

/// \brief write relevant information for photon conversion analysis to an AO2D.root file. This file is then the only necessary input to perform
/// pcm analysis.
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

// runme like: o2-analysis-trackselection -b --aod-file ${sourceFile} --aod-writer-json ${writerFile} | o2-analysis-timestamp -b | o2-analysis-trackextension -b | o2-analysis-lf-lambdakzerobuilder -b | o2-analysis-pid-tpc -b | o2-analysis-em-skimmermc -b

// todo: remove reduantant information in GammaConversionsInfoTrue
#include "gammaTables.h"

#include "TVector3.h"


#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using tracksAndMcLabels = soa::Join<aod::Tracks, aod::McTrackLabels>;


// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
//~ using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

struct skimmerGammaConversionsTruthOnlyMc {

  Produces<aod::StoredMcGammasInfoTrue> fFuncTableMcGammas;
  Produces<aod::MCGammaDaughters> fFuncTableMcGammaDaughters;
  
  HistogramRegistry registry{
    "registry",
    {
      {"hMcCollisionZ", "hMcCollisionZ", {HistType::kTH1F, {{800, -10.f, 10.f}}}},
      {"hEtaDiff", "hEtaDiff", {HistType::kTH1F, {{400, -2.f, 2.f}}}},

    },
  };
  
  // loop over MC truth McCollisions
  void process(aod::McCollision const& theMcCollision,
               aod::McParticles const& theMcParticles)
  {
    registry.fill(HIST("hMcCollisionZ"), theMcCollision.posZ());
    for (auto &lMcParticle : theMcParticles) {
      if (lMcParticle.pdgCode() == 22) {
        size_t lNDaughters = 0;
        if (lMcParticle.has_daughters()) {
          for (auto &lDaughter : lMcParticle.daughters_as<aod::McParticles>()) {
            ++lNDaughters;
            
            TVector3 lDaughterVtx(lDaughter.vx(),lDaughter.vy(), lDaughter.vz());
            if (lMcParticle.isPhysicalPrimary()){
              float_t lEtaDiff = lDaughterVtx.Eta() - lMcParticle.eta();
              registry.fill(HIST("hEtaDiff"), lEtaDiff);
            }
            
            fFuncTableMcGammaDaughters(lMcParticle.mcCollisionId(),
                                       lMcParticle.globalIndex(),
                                       //~ lDaughter.mothersIds().size() ? lDaughter.mothersIds()[0] : -1000, // SFS this is potentially unsafe, what if there are more mothers?, or could this still point to 0 even though it is a daughter? todo: make this cleaner
                                       lDaughter.mothersIds().size(),
                                       lDaughter.pdgCode(),
                                       lDaughter.vx(), lDaughter.vy(), lDaughter.vz(),
                                       lDaughter.eta(),
                                       lDaughter.p(),
                                       lDaughter.phi(),
                                       lDaughter.pt());
          }
        }
        fFuncTableMcGammas(
            lMcParticle.mcCollisionId(),
            lMcParticle.globalIndex(),
            -1,
            lMcParticle.statusCode(),
            lMcParticle.flags(),
            lMcParticle.px(), lMcParticle.py(), lMcParticle.pz(), 
            lMcParticle.vx(), lMcParticle.vy(), lMcParticle.vz(), lMcParticle.vt(), 
            lMcParticle.template daughters_as<aod::McParticles>().size());
            
        //~ fFuncTableMcGammas(lMcParticle.mcCollisionId(),
                           //~ lMcParticle.globalIndex(),
                           //~ lMcParticle.isPhysicalPrimary(),
                           //~ lNDaughters,
                           //~ lMcParticle.eta(),
                           //~ lMcParticle.p(),
                           //~ lMcParticle.phi(),
                           //~ lMcParticle.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversionsTruthOnlyMc>(cfgc)};
}
