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

/// \brief write relevant information for photon conversion analysis to AO2D.root file. This file is then to the used as the only input to perform
/// the actual pcm analysis.
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/StrangenessTables.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h" // for BigTracks

#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTPC.h"

//~ #include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

namespace o2::aod
{

namespace gammatrackreco
{
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, tpcFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(IsFromConversionPhoton, isFromConversionPhoton, bool);
}

DECLARE_SOA_TABLE(GammaConversionTracks, "AOD", "V0TRACKS",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  gammatrackreco::IsFromConversionPhoton,
                  gammatrackreco::TpcFoundOverFindableCls,
                  gammatrackreco::TpcCrossedRowsOverFindableCls,
                  track::P,
                  track::TPCSignal,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaPi);

namespace gammamctrue
{
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
} // gammamctrue

// SFS todo: need to add some sort of indexing here to get from recos to true table
DECLARE_SOA_TABLE(GammaConversionsInfoTrue, "AOD", "V0INFOTRUE",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  v0data::X, 
                  v0data::Y, 
                  v0data::Z,
                  gammamctrue::Eta,
                  gammamctrue::Phi,
                  gammamctrue::Pt,
                  v0data::V0Radius<v0data::X, v0data::Y>);
} // namespace o2::aod

struct SkimmerMc {
  // ============================ DEFINITION OF CUT VARIABLES =============================================
  
  // ============================ DEFINITION OF HISTOGRAMS ================================================
    
  // ============================ TABLES TO WRITTEN TO AO2D.root ==========================================
  Produces<aod::GammaConversionsInfoTrue> fFuncTableV0InfoTrue;
  Produces<aod::GammaConversionTracks> fFuncTableGammaTracks;
  
  // ============================ FUNCTION DEFINITIONS ====================================================
  void process(aod::Collisions::iterator  const &theCollision,
               aod::V0s                  const& V0s,
               aod::V0Datas               const &theV0s,
               tracksAndTPCInfoMC         const &theTracks,
               aod::McParticles           const &theMcParticles)
  {
    auto fillTrackTable = [&](auto &theV0, auto &theTrack, bool theIsFromConversionPhoton){
      fFuncTableGammaTracks(
        theV0.v0(),
        theIsFromConversionPhoton,
        theTrack.tpcFoundOverFindableCls(),
        theTrack.tpcCrossedRowsOverFindableCls(),
        theTrack.p(),
        theTrack.tpcSignal(),
        theTrack.tpcNSigmaEl(),
        theTrack.tpcNSigmaPi());
    };
    
    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter    
      
      bool lIsConversionPhoton = isConversionPhoton(lV0,
                                                    lTrackPos,
                                                    lTrackNeg,
                                                    theMcParticles);
    
      fillTrackTable(lV0, lTrackPos, lIsConversionPhoton);
      fillTrackTable(lV0, lTrackNeg, lIsConversionPhoton);
    }
  }
  
  template <typename TV0, typename TTRACK, typename TMC>
  bool isConversionPhoton(const TV0& theV0, const TTRACK& theTrackPos, const TTRACK& theTrackNeg, const TMC& theMcParticles)
  {
    bool result = false;
    // todo: verify it is enough to check only mother0 being equal

    /* example from https://aliceo2group.github.io/analysis-framework/docs/tutorials/indexTables.html?highlight=_as:
     * track0.collision_as<myCol>().mult() : access multiplicity of collission associated with track0
     */

    //~ // use & here?
    auto lMcPos = theTrackPos.template mcParticle_as<aod::McParticles_001>();
    auto lMcNeg = theTrackNeg.template mcParticle_as<aod::McParticles_001>();
    
    //! Mother tracks (possible empty) array. Iterate over mcParticle.mothers_as<aod::McParticles>())
    std::vector<int> lMothers;
    // SFS todo: remove all those mothers_as<aod::McParticles_001>, are the loops even necesarry?
    for (auto& mP : lMcPos.template mothers_as<aod::McParticles_001>()) {
      LOGF(info, "   mother index mP: %d", mP.globalIndex());
      lMothers.push_back(mP.globalIndex());
    }

    if (lMothers.size() > 0) {
      for (auto& mN : lMcNeg.template mothers_as<aod::McParticles_001>()) {
        LOGF(info, "   mother index mN: %d", mN.globalIndex());
        lMothers.push_back(mN.globalIndex());
      }
    }
    
    // SFS verify theyre all the same category and remove
    int lSame = 0;
    {
      if (lMothers.size() == 2) {
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size2: 01");
          lSame = 1;
        }
      }

      if (lMothers.size() == 3) {
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size3: 01");
          lSame = 2;
        }
        if (lMothers[0] == lMothers[2]) {
          LOGF(info, "size2: 02");
          lSame = 3;
        }
        if (lMothers[1] == lMothers[2]) {
          LOGF(info, "size2: 12");
          lSame = 4;
        }
      }

      if (lMothers.size() == 4) {
        if (lMothers[0] == lMothers[2]) {
          LOGF(info, "size4 02");
          lSame = 4;
        }
        if (lMothers[1] == lMothers[3]) {
          LOGF(info, "size4 13");
          lSame = 5;
        }
        if (lMothers[0] == lMothers[3]) {
          LOGF(info, "size4 03");
          lSame = 6;
        }
        if (lMothers[1] == lMothers[2]) {
          LOGF(info, "size4 12");
          lSame = 7;
        }
        if (lMothers[0] == lMothers[1]) {
          LOGF(info, "size4 01");
          lSame = 8;
        }
        if (lMothers[2] == lMothers[3]) {
          LOGF(info, "size4 23");
          lSame = 9;
        }
      }
    }

    if (lSame) {
      // SFS todo: actually no loop required here, for this
      for (auto& lMother : lMcNeg.template mothers_as<aod::McParticles_001>()) {

        if ((result = lMother.pdgCode()==22)) {
          fFuncTableV0InfoTrue(
            theV0.v0(),
            lMcPos.vx(),
            lMcPos.vy(), 
            lMcPos.vz(),
            lMother.eta(),
            lMother.phi(),
            lMother.pt());
        }
      }
    }
    return result;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SkimmerMc>(cfgc)};
}
