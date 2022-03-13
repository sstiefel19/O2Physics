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
/*namespace gammatrackspos
{
DECLARE_SOA_COLUMN(Pt, pospt, float);
DECLARE_SOA_COLUMN(Eta, posmaxpt, float);
DECLARE_SOA_COLUMN(P, posp, float);
DECLARE_SOA_COLUMN(TpcNSigmaEl, postpcnsigmael, float);
DECLARE_SOA_COLUMN(TpcNSigmaPi, postpcnsigmapi, float);
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, postpcfoundoverfindablecls, float);
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, postpccrossedrowsoverfindablecls, float);
} // namespace gammatrackspos
namespace gammatracksneg
{
DECLARE_SOA_COLUMN(Pt, negpt, float);
DECLARE_SOA_COLUMN(Eta, negmaxpt, float);
DECLARE_SOA_COLUMN(P, negp, float);
DECLARE_SOA_COLUMN(TpcNSigmaEl, negtpcnsigmael, float);
DECLARE_SOA_COLUMN(TpcNSigmaPi, negtpcnsigmapi, float);
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, negtpcfoundoverfindablecls, float);
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, negtpccrossedrowsoverfindablecls, float);
} // namespace gammatracksneg 

DECLARE_SOA_TABLE(GammaConversionsTrackDataPos, "AOD", "GACOTRPOS",
                  gammatrackspos::Pt,
                  gammatrackspos::Eta,
                  gammatrackspos::P,
                  gammatrackspos::TpcNSigmaEl,
                  gammatrackspos::TpcNSigmaPi,
                  gammatrackspos::TpcFoundOverFindableCls,
                  gammatrackspos::TpcCrossedRowsOverFindableCls);

DECLARE_SOA_TABLE(GammaConversionsTrackDataNeg, "AOD", "GACOTRNEG",                  
                  gammatracksneg::Pt,
                  gammatracksneg::Eta,
                  gammatracksneg::P,
                  gammatracksneg::TpcNSigmaEl,
                  gammatracksneg::TpcNSigmaPi,
                  gammatracksneg::TpcFoundOverFindableCls,
                  gammatracksneg::TpcCrossedRowsOverFindableCls);
*/

namespace gammatracksreco
{
DECLARE_SOA_COLUMN(PosPt, pospt, float);
DECLARE_SOA_COLUMN(PosEta, posmaxpt, float);
DECLARE_SOA_COLUMN(PosP, posp, float);
DECLARE_SOA_COLUMN(PosTpcNSigmaEl, postpcnsigmael, float);
DECLARE_SOA_COLUMN(PosTpcNSigmaPi, postpcnsigmapi, float);
DECLARE_SOA_COLUMN(PosTpcFoundOverFindableCls, postpcfoundoverfindablecls, float);
DECLARE_SOA_COLUMN(PosTpcCrossedRowsOverFindableCls, postpccrossedrowsoverfindablecls, float);

DECLARE_SOA_COLUMN(NegPt, negpt, float);
DECLARE_SOA_COLUMN(NegEta, negmaxpt, float);
DECLARE_SOA_COLUMN(NegP, negp, float);
DECLARE_SOA_COLUMN(NegTpcNSigmaEl, negtpcnsigmael, float);
DECLARE_SOA_COLUMN(NegTpcNSigmaPi, negtpcnsigmapi, float);
DECLARE_SOA_COLUMN(NegTpcFoundOverFindableCls, negtpcfoundoverfindablecls, float);
DECLARE_SOA_COLUMN(NegTpcCrossedRowsOverFindableCls, negtpccrossedrowsoverfindablecls, float);

DECLARE_SOA_COLUMN(IsConversionPhoton, isconversionphoton, bool);
} // namespace gammatracksreco

namespace gammamctrue
{
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
} // gammamctrue

DECLARE_SOA_TABLE(GammaConversionsTracksReco, "AOD", "V0TRACKSRECO",
                  gammatracksreco::PosPt,
                  gammatracksreco::PosEta,
                  gammatracksreco::PosP,
                  gammatracksreco::PosTpcNSigmaEl,
                  gammatracksreco::PosTpcNSigmaPi,
                  gammatracksreco::PosTpcFoundOverFindableCls,
                  gammatracksreco::PosTpcCrossedRowsOverFindableCls, 
                  gammatracksreco::NegPt,
                  gammatracksreco::NegEta,
                  gammatracksreco::NegP,
                  gammatracksreco::NegTpcNSigmaEl,
                  gammatracksreco::NegTpcNSigmaPi,
                  gammatracksreco::NegTpcFoundOverFindableCls,
                  gammatracksreco::NegTpcCrossedRowsOverFindableCls,
                  gammatracksreco::IsConversionPhoton);

DECLARE_SOA_TABLE(GammaConversionsTrue, "AOD", "V0STRUE",
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
  Produces<aod::GammaConversionsTracksReco> fFuncTableTracksRecoData;
  Produces<aod::GammaConversionsTrue> fFuncTableGammasTrueData;

  // ============================ FUNCTION DEFINITIONS ====================================================
  void process(aod::Collisions::iterator  const &theCollision,
               aod::V0Datas               const &theV0s,
               tracksAndTPCInfoMC         const &theTracks,
               aod::McParticles           const &theMcParticles)
  {
    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter
      
      fFuncTableTracksRecoData(lTrackPos.pt(), 
                               lTrackPos.eta(),
                               lTrackPos.p(),
                               lTrackPos.tpcNSigmaEl(),
                               lTrackPos.tpcNSigmaPi(),
                               lTrackPos.tpcFoundOverFindableCls(),
                               lTrackPos.tpcCrossedRowsOverFindableCls(),
                               lTrackNeg.pt(), 
                               lTrackNeg.eta(),
                               lTrackNeg.p(),
                               lTrackNeg.tpcNSigmaEl(),
                               lTrackNeg.tpcNSigmaPi(),
                               lTrackNeg.tpcFoundOverFindableCls(),
                               lTrackNeg.tpcCrossedRowsOverFindableCls(),
                               isPhoton(lV0,
                                        lTrackPos,
                                        lTrackNeg,
                                        theMcParticles));
    }
  }
  
  template <typename TV0, typename TTRACK, typename TMC>
  bool isPhoton(const TV0& theV0, const TTRACK& theTrackPos, const TTRACK& theTrackNeg, const TMC& theMcParticles)
  {
    bool lIsPhoton = false;
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

        if ((lIsPhoton = lMother.pdgCode()==22)) {
          fFuncTableGammasTrueData(lMcPos.vx(),
                                   lMcPos.vy(), 
                                   lMcPos.vz(),
                                   lMother.eta(),
                                   lMother.phi(),
                                   lMother.pt());
        }
      }
    }
    return lIsPhoton;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SkimmerMc>(cfgc)};
}
