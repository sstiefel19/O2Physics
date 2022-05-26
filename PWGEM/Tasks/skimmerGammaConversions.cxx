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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
using tracksAndTPCInfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi>;
using tracksAndTPCInfoMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended, aod::pidTPCEl, aod::pidTPCPi, aod::McTrackLabels>;

struct skimmerGammaConversions {

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hCollisionZ_MCRec", "hCollisionZ_MCRec", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_all_MCTrue", "hCollisionZ_all_MCTrue", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hCollisionZ_MCTrue", "hCollisionZ_MCTrue", {HistType::kTH1F, {{800, -50.f, 50.f}}}},
      {"hMcParticlesSize", "hMcParticlesSize", {HistType::kTH1F, {{100, 0.f, 1000000.f}}}},

      //~ {"hMotherSizes", "hMotherSizes", {HistType::kTH1F, {{14, 0.f, 14.f}}}},
      {"hPeculiarOccurences", "hPeculiarOccurences", {HistType::kTH1F, {{50, -25.f, 25.f}}}},
    },
  };

  std::shared_ptr<TH1> fMotherSizesHisto{};
      //~ std::shared_ptr<TH1> lHisto = getTH<TH1>(theMap, theName);


  Produces<aod::V0DaughterTracks> fFuncTableV0DaughterTracks;
  Produces<aod::McGammasTrue> fFuncTableMcGammasFromConfirmedV0s;

  void init(InitContext const&)
  {
    //~ HistPtr lHistPtr = fRegistry.add("hMotherSizes", "hMotherSizes", {HistType::kTH1F, {{14, 0.f, 14.f}}});
    fMotherSizesHisto = std::get<std::shared_ptr<TH1>>(fRegistry.add("hMotherSizes", "hMotherSizes", {HistType::kTH1F, {{28, -14.f, 14.f}}}));
    //~ fMotherSizesHisto = std::get<std::shared_ptr<TH1F*>>(fRegistry.add("hMotherSizes", "hMotherSizes", {HistType::kTH1F, {{14, 0.f, 14.f}}}));
  }

  // ============================ FUNCTION DEFINITIONS ====================================================
  void processRec(aod::Collisions::iterator const& theCollision,
                  aod::V0s,
                  aod::V0Datas const& theV0s,
                  tracksAndTPCInfo const& theTracks)
  {
    auto fillTrackTable = [&](auto& theV0, auto& theTrack, bool theIsPositive, bool theIsFromConversionPhoton) {
      fFuncTableV0DaughterTracks(
        theV0.v0Id(),
        theIsFromConversionPhoton,
        theTrack.dcaXY(),
        theTrack.eta(),
        theTrack.p(),
        theTrack.phi(),
        theTrack.pt(),
        theIsPositive,
        theTrack.tpcCrossedRowsOverFindableCls(),
        theTrack.tpcFoundOverFindableCls(),
        theTrack.tpcNClsCrossedRows(),
        theTrack.tpcNSigmaEl(),
        theTrack.tpcNSigmaPi(),
        theTrack.tpcSignal());
    };

    fRegistry.fill(HIST("hCollisionZ"), theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfo>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfo>(); // negative daughter

      fillTrackTable(lV0, lTrackPos, true, false);
      fillTrackTable(lV0, lTrackNeg, false, false);
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processRec, "process reconstructed info only", true);

  void processMc(aod::McCollision const& theMcCollision,
                 soa::SmallGroups<soa::Join<aod::McCollisionLabels,
                                            aod::Collisions>> const& theCollisions,
                 aod::V0s,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const& theTracks,
                 aod::McParticles const& theMcParticles)
  {
    auto fillTrackTable = [&](auto& theV0, auto& theTrack, bool theIsPositive, bool theIsFromConversionPhoton) {
      fFuncTableV0DaughterTracks(
        theV0.v0Id(),
        theIsFromConversionPhoton,
        theTrack.dcaXY(),
        theTrack.eta(),
        theTrack.p(),
        theTrack.phi(),
        theTrack.pt(),
        theIsPositive,
        theTrack.tpcCrossedRowsOverFindableCls(),
        theTrack.tpcFoundOverFindableCls(),
        theTrack.tpcNClsCrossedRows(),
        theTrack.tpcNSigmaEl(),
        theTrack.tpcNSigmaPi(),
        theTrack.tpcSignal());
    };

    fRegistry.fill(HIST("hCollisionZ_all_MCTrue"), theMcCollision.posZ());
    if (theCollisions.size() == 0) {
      return;
    }

    fRegistry.fill(HIST("hCollisionZ_MCTrue"), theMcCollision.posZ());
    fRegistry.fill(HIST("hMcParticlesSize"), theMcParticles.size());

    for (auto& lCollision : theCollisions) {
      fRegistry.fill(HIST("hCollisionZ_MCRec"), lCollision.posZ());
LOGF(info, "collision");
      // todo: replace by sliceByCached
      auto lGroupedV0s = theV0s.sliceBy(aod::v0data::collisionId, lCollision.globalIndex());
      for (auto& lV0 : lGroupedV0s) {

        auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
        auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter

        bool lIsConversionPhoton = isConversionPhoton(lV0,
                                                      lTrackPos,
                                                      lTrackNeg);

        fillTrackTable(lV0, lTrackPos, true, lIsConversionPhoton);
        fillTrackTable(lV0, lTrackNeg, false, lIsConversionPhoton);
      }
    }
  }
  PROCESS_SWITCH(skimmerGammaConversions, processMc, "process reconstructed and mc info ", false);

  // SFS todo: make pretty and short
  template <typename TV0, typename TTRACK>
  bool isConversionPhoton(TV0 const& theV0,
                          TTRACK const& theTrackPos,
                          TTRACK const& theTrackNeg)
  {
    auto trackHasMcParticle = [&](auto const& theRecTrack){
      if (!theRecTrack.has_mcParticle()){
        LOGF(info,
             "SFS V0 daughter track %d is a track without mc particle. Can't be a confirmable v0.",
             theRecTrack.globalIndex());
        fMotherSizesHisto->Fill(-0.5);
        return false;
      }
      return true;
    };

    auto getMotherIndeces = [&](auto const& theMcParticle){
      std::vector<int> lMothersIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(info, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      fMotherSizesHisto->Fill(0.5 + (float)lMothersIndeces.size());
      return std::move(lMothersIndeces);
    };
    // ================== end lambda definitions =========================

    // one of the tracks doesnt have a mcParticle, has been accounted for in fMotherSizesHisto
    if (!(trackHasMcParticle(theTrackPos) && trackHasMcParticle(theTrackNeg))){
      return false;
    }

    // get mcParticles
    auto lMcPos = theTrackPos.mcParticle();
    auto lMcNeg = theTrackNeg.mcParticle();

    // get indeces of mcMother of tracks
    std::vector<int> lMothersIndecesPos = getMotherIndeces(lMcPos);
    std::vector<int> lMothersIndecesNeg = getMotherIndeces(lMcNeg);

    // none of tracks has a mother, has been accounted for in fMotherSizesHisto
    if (!(lMothersIndecesPos.size() || lMothersIndecesNeg.size())){
      return false;
    }

    // exactly one has a mother
    if ((lMothersIndecesPos.size() + lMothersIndecesNeg.size()) == 1){
      fMotherSizesHisto->Fill(-1.5);
      return false;
    }

    // we know now both tracks have at least one mother
    // check if it is the same
    if (lMothersIndecesPos[0] != lMothersIndecesNeg[0]){
      fMotherSizesHisto->Fill(-2.5);
      return false;
    }

    // both tracks have the first mother
    // SFS todo: actually no loop required here, for this
    bool lResult = false;
    for (auto& lMcMother : lMcNeg.template mothers_as<aod::McParticles>()) {

      // case we have a confirmed photon conversion and have the mother mc particle
      // todo: store if it was not a photon but another mother particle
      if ((lResult = lMcMother.pdgCode() == 22)) {

        // get mc daughter in order to compute true conversion point
        auto lDaughters = lMcMother.template daughters_as<aod::McParticles>();
        if (lDaughters.begin() == lDaughters.end()) {
          // mc converted mother has no mc daughters, should never happen
          fMotherSizesHisto->Fill(-3.5);
          return false;
        }
        auto lDaughter0 = lDaughters.begin();
        float lDaughter0Vx = lDaughter0.vx();
        float lDaughter0Vy = lDaughter0.vy();
        float lDaughter0Vz = lDaughter0.vz();
        float lV0Radius = sqrt(pow(lDaughter0Vx, 2) + pow(lDaughter0Vy, 2));

        fFuncTableMcGammasFromConfirmedV0s(
          lMcMother.mcCollisionId(),
          lMcMother.globalIndex(),
          theV0.v0Id(),
          lMcMother.statusCode(),
          lMcMother.flags(),
          lMcMother.px(), lMcMother.py(), lMcMother.pz(),
          lMcMother.vx(), lMcMother.vy(), lMcMother.vz(), lMcMother.vt(),
          lDaughters.size(),
          lMcMother.eta(), lMcMother.phi(), lMcMother.p(), lMcMother.pt(), lMcMother.y(),
          lDaughter0Vx, lDaughter0Vy, lDaughter0Vz,
          lV0Radius);
      }
      break; // because we only want to look at the first mother. If there are more it will show up in fMotherSizesHisto
    }

    return lResult;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerGammaConversions>(cfgc)};
}
