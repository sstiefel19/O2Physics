// SFS todo: use existing columns instead creating new ones!

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

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/StrangenessTables.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h" // for BigTracks

#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTPC.h"

#include <TH1.h>
#include <TH1F.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>

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

namespace gammatracks
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
} // namespace gammatracks

/* SFS todo: + index auf gammasinfotrue table for true indices*/

DECLARE_SOA_TABLE(GammaConversionsTrackData, "AOD", "GACONVTR",
                  gammatracks::PosPt,
                  gammatracks::PosEta,
                  gammatracks::PosP,
                  gammatracks::PosTpcNSigmaEl,
                  gammatracks::PosTpcNSigmaPi,
                  gammatracks::PosTpcFoundOverFindableCls,
                  gammatracks::PosTpcCrossedRowsOverFindableCls, 
                  gammatracks::NegPt,
                  gammatracks::NegEta,
                  gammatracks::NegP,
                  gammatracks::NegTpcNSigmaEl,
                  gammatracks::NegTpcNSigmaPi,
                  gammatracks::NegTpcFoundOverFindableCls,
                  gammatracks::NegTpcCrossedRowsOverFindableCls);
} // namespace o2::aod

#include "Framework/runDataProcessing.h"

struct SkimmerMc {
  // ============================ DEFINITION OF CUT VARIABLES ==========================================================
  Configurable<bool> fDoEventSel{"fDoEventSel", 0, "demand sel7 for events"};
  Configurable<float> fCentMin{"fCentMin", 0.0, "lower bound of centrality selection"};
  Configurable<float> fCentMax{"fCentMax", 100.0, "upper bound of centrality selection"};
  Configurable<float> fSinglePtCut{"fSinglePtCut", 0.04, "minimum daughter track pt"};
  Configurable<float> fEtaCut{"fEtaCut", 0.8, "accepted eta range"};
  Configurable<float> fMinR{"fMinR", 5., "minimum conversion radius of the V0s"};
  Configurable<float> fMaxR{"fMaxR", 180., "maximum conversion radius of the V0s"};
  Configurable<float> fPIDnSigmaBelowElectronLine{"fPIDnSigmaBelowElectronLine", -3., "minimum sigma electron PID for V0 daughter tracks"};
  Configurable<float> fPIDnSigmaAboveElectronLine{"fPIDnSigmaAboveElectronLine", 3., "maximum sigma electron PID for V0 daughter tracks"};
  Configurable<float> fPIDnSigmaAbovePionLine{"fPIDnSigmaAbovePionLine", 3., "minimum sigma to be over the pion line for low momentum tracks"}; //case 4: 3.0sigma, 1.0 sigma at high momentum
  Configurable<float> fPIDnSigmaAbovePionLineHighP{"fPIDnSigmaAbovePionLineHighP", 1., "minimum sigma to be over the pion line for high momentum tracks"};
  Configurable<float> fPIDMinPnSigmaAbovePionLine{"fPIDMinPnSigmaAbovePionLine", 0.4, "minimum track momentum to apply any pion rejection"}; //case 7:  // 0.4 GeV
  Configurable<float> fPIDMaxPnSigmaAbovePionLine{"fPIDMaxPnSigmaAbovePionLine", 8., "border between low and high momentum pion rejection"}; //case 7:  // 8. GeV
  Configurable<float> fMinTPCFoundOverFindableCls{"fMinTPCNClsFoundOverFindable", 0.6, "minimum ratio found tpc clusters over findable"}; //case 9:  // 0.6
  Configurable<float> fMinTPCCrossedRowsOverFindableCls{"fMinTPCCrossedRowsOverFindableCls", 0.0, "minimum ratio TPC crossed rows over findable clusters"};
  Configurable<float> fQtPtMax{"fQtPtMax", 0.11, "up to fQtMax, multiply the pt of the V0s by this value to get the maximum qt "};
  Configurable<float> fQtMax{"fQtMax", 0.040, "maximum qt"};
  Configurable<float> fMaxPhotonAsymmetry{"fMaxPhotonAsymmetry", 0.95, "maximum photon asymetry"};
  Configurable<float> fPsiPairCut{"fPsiPairCut", 0.1, "maximum psi angle of the track pair"};
  Configurable<float> fCosPAngleCut{"fCosPAngleCut", 0.85, "mimimum cosinus of the pointing angle"}; // case 4


  // ============================ DEFINITION OF HISTOGRAMS ==========================================================
  
  // SFS todo: struct MyHistogramSpec : public HistogramSpec {
  struct MyHistogramSpec {
    MyHistogramSpec(char const* const name_, char const* const title_, HistogramConfigSpec config_, bool callSumw2_ = false)
      : name(name_),
        hash(compile_time_hash(name_)),
        title(title_),
        config(config_),
        callSumw2(callSumw2_)
    {
    }

    MyHistogramSpec(char const* const name_, HistogramConfigSpec config_, bool callSumw2_ = false)
      : name(name_),
        hash(compile_time_hash(name_)),
        title(name_),
        config(config_),
        callSumw2(callSumw2_)
    {
    }

    MyHistogramSpec()
      : name(""),
        hash(0),
        config()
    {
    }
    MyHistogramSpec(MyHistogramSpec const& other) = default;
    MyHistogramSpec(MyHistogramSpec&& other) = default;

    std::string name{};
    uint32_t hash{};
    std::string title{};
    HistogramConfigSpec config{};
    bool callSumw2{}; // wether or not hist needs heavy error structure produced by Sumw2()
  };

  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};

  // track histograms
  std::vector<MyHistogramSpec> fTrackHistoDefinitions{
    {"hTPCFoundOverFindableCls", {HistType::kTH1F, {{800, 0.9f, 1.01f}}}},
    {"hTPCCrossedRowsOverFindableCls", {HistType::kTH1F, {{800, 0.8f, 1.5f}}}},

    {"hTPCdEdx", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, 0.f, 200.f}}}},
    {"hTPCdEdxSigEl", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, -10.f, 10.f}}}},
    {"hTPCdEdxSigPi", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, -10.f, 10.f}}}}};

  // v0 histograms
  std::vector<MyHistogramSpec> fV0HistoDefinitions{
    {"hPtRec", {HistType::kTH1F, {{800, 0.0f, 25.0f}}}},
    {"hEtaRec", {HistType::kTH1F, {{800, -2.f, 2.f}}}},
    {"hPhiRec", {HistType::kTH1F, {{800, 0.f, 2.f * M_PI}}}},
    {"hConvPointR", {HistType::kTH1F, {{800, 0.f, 200.f}}}},
    {"hArmenteros", {HistType::kTH2F, {{800, -1.f, 1.f}, {800, 0.f, 0.25f}}}},
    {"hPsiPtRec", {HistType::kTH2F, {{800, -2.f, 2.f}, {800, 0.f, 10.f}}}},
    {"hCosPAngle", {HistType::kTH1F, {{800, 0.99f, 1.005f}}}},

    {"IsPhotonSelected", {HistType::kTH1F, {{13, -0.0f, 12.5f}}}} // only in afterCuts
  };

  // only in mc
  // resolution histos
  std::vector<MyHistogramSpec> fV0ResolutionHistoDefinitions{
    {"hPtRes", "hPtRes_Rec-MC", {HistType::kTH1F, {{800, -0.5f, 0.5f}}}},
    {"hEtaRes", "hEtaRes_Rec-MC", {HistType::kTH1F, {{800, -0.05f, 0.05f}}}},
    {"hPhiRes", "hPhiRes_Rec-MC", {HistType::kTH1F, {{800, -0.05f, 0.05f}}}},
    {"hConvPointRRes", "hConvPointRRes_Rec-MC", {HistType::kTH1F, {{800, -200.f, 200.f}}}},
    {"hConvPointAbsoluteDistanceRes", "hConvPointAbsoluteDistanceRes", {HistType::kTH1F, {{800, -0.0f, 200.f}}}},
  };

  // Reconstructed info of MC validated V0s histos fV0ReconstructedInfoMCValidatedHistos
  std::vector<MyHistogramSpec> fV0ReconstructedInfoMCValidatedHistoDefinitions{
    {"hValidatedPtRec", {HistType::kTH1F, {{800, 0.0f, 25.f}}}},
    {"hValidatedRRec", {HistType::kTH1F, {{800, 0.0f, 250.f}}}},
  };

  // MC info only histos
  std::vector<MyHistogramSpec> fV0MCInfoOnlyHistoDefinitions{
    {"hValidatedPtTrue", {HistType::kTH1F, {{800, -0.0f, 25.f}}}},
  };

  // ============================ CONTAINERS TO HOLD POINTERS TO HISTOGRAMS ==========================================
  std::map<std::string, HistPtr> fTrackHistos;
  std::map<std::string, HistPtr> fV0Histos;
  std::map<std::string, HistPtr> fV0ResolutionHistos;
  std::map<std::string, HistPtr> fV0ReconstructedInfoMCValidatedHistos;
  std::map<std::string, HistPtr> fV0MCInfoOnlyHistos;

  std::string fPrefixReconstructedInfoHistos{"reconstructedInformationOnly/"};
  std::string fPrefixMCInfoNeededHistos{"MCinformationNeeded/"};
  std::string fFullNameIsPhotonSelectedHisto{fPrefixReconstructedInfoHistos + "v0/afterCuts/IsPhotonSelected"};
  
  std::map<std::string, size_t> fPhotonCutIndeces{
    {"kPhotonIn", 0},
    {"kTrackEta", 1},
    {"kTrackPt", 2},
    {"kElectronPID", 3},
    {"kPionRejLowMom", 4},
    {"kPionRejHighMom", 5},
    {"kTPCFoundOverFindableCls", 6},
    {"kTPCCrossedRowsOverFindableCls", 7},
    {"kV0Radius", 8},
    {"kArmenteros", 9},
    {"kPsiPair", 10},
    {"kCosinePA", 11},
    {"kPhotonOut", 12}};

  Produces<aod::GammaConversionsTrackData> fFuncTableTracks;

  //Produces<aod::GammaConversionsData> fFuncConversionData;
  //Produces<aod::GammaConversionsTrueData> fFuncConversionTrueData;
  //Produces<aod::StoredV0Datas> fFuncConversionTrueData;
  /*Produces<aod::GammaConversionsTrackDataPos> fFuncTableTracksPos;
  Produces<aod::GammaConversionsTrackDataNeg> fFuncTableTracksNeg;*/

  void init(InitContext const&)
  {
    for (auto bac : std::vector<std::string>{"beforeCuts/", "afterCuts/"}) {
      
      // reconstructed data histograms
      {
        // track histograms
        {
          std::string lPath(fPrefixReconstructedInfoHistos + "track/" + bac);
          for (auto& tHisto : fTrackHistoDefinitions) {
            std::string lFullName(lPath + tHisto.name);
            LOGF(debug, "adding %s", lFullName);
            fTrackHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }

        // v0 histograms
        {
          std::string lPath(fPrefixReconstructedInfoHistos + "v0/" + bac);
          for (auto tHisto : fV0HistoDefinitions) {
            if (tHisto.name == "IsPhotonSelected" && bac == "beforeCuts/") {
              continue;
            }
            std::string lFullName(lPath + tHisto.name);
                        LOGF(debug, "adding %s", lFullName);

            fV0Histos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }
      }

      // MC information histos
      {
        // v0 Resolution histos
        {
          std::string lPath(fPrefixMCInfoNeededHistos + "v0/resolutions/" + bac);
          for (auto tHisto : fV0ResolutionHistoDefinitions) {
            std::string lFullName(lPath + tHisto.name);
            LOGF(debug, "adding %s", lFullName);
            fV0ResolutionHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }

        // reconstructed info of MC validated V0s histos
        {
          {
            std::string lPath(fPrefixMCInfoNeededHistos + "v0/reconstructedInfoOfMCvalidated/" + bac);
            for (auto tHisto : fV0ReconstructedInfoMCValidatedHistoDefinitions) {
              std::string lFullName(lPath + tHisto.name);
              LOGF(debug, "adding %s", lFullName);
              fV0ReconstructedInfoMCValidatedHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
            }
          }
        }

        // v0 Info only histos
        {
          std::string lPath(fPrefixMCInfoNeededHistos + "v0/MCinformationOnly/" + bac);
          for (auto tHisto : fV0MCInfoOnlyHistoDefinitions) {
            std::string lFullName(lPath + tHisto.name);
            LOGF(debug, "adding %s", lFullName);
            fV0MCInfoOnlyHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }
      }
    }

    {
      // do some labeling
      auto lIsPhotonSelectedHisto = fV0Histos.find(fFullNameIsPhotonSelectedHisto);
      if (lIsPhotonSelectedHisto != fV0Histos.end()) {
        TAxis* lXaxis = std::get<std::shared_ptr<TH1>>(lIsPhotonSelectedHisto->second)->GetXaxis();
        for (auto& lPairIt : fPhotonCutIndeces) {
          lXaxis->SetBinLabel(lPairIt.second + 1, lPairIt.first.data());
        }
      }
    }
  }

  // SFS todo: think about if this is actually too expensive. Going the other way round with the indices as keys wouldnt require lookups at inserting but pbly produce a but of code duplication at the definition of the cut names
  size_t gMax_size = (size_t)-1;
  size_t getPhotonCutIndex(std::string const& theKey)
  {
    auto lPairIt = fPhotonCutIndeces.find(theKey);
    if (lPairIt != fPhotonCutIndeces.end()) {
      return lPairIt->second;
    }
    return gMax_size;
  }

template <typename TCOLL, typename TV0, typename TTRACK>
bool processPhoton(TCOLL const& theCollision, const TV0& theV0, const TTRACK& theTrackPos, const TTRACK& theTrackNeg)
{
  
  // apply track cuts
  if (!(trackPassesCuts(theTrackPos) && trackPassesCuts(theTrackNeg))) {
    return kFALSE;
  }
  
  // apply photon cuts
  float lV0CosinePA = theV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());
  if (!passesPhotonCuts(theV0, lV0CosinePA)) {
    return kFALSE;
  }
  
  fFuncTableTracks(theTrackPos.pt(), 
                   theTrackPos.eta(),
                   theTrackPos.p(),
                   theTrackPos.tpcNSigmaEl(),
                   theTrackPos.tpcNSigmaPi(),
                   theTrackPos.tpcFoundOverFindableCls(),
                   theTrackPos.tpcCrossedRowsOverFindableCls(),
                   theTrackNeg.pt(), 
                   theTrackNeg.eta(),
                   theTrackNeg.p(),
                   theTrackNeg.tpcNSigmaEl(),
                   theTrackNeg.tpcNSigmaPi(),
                   theTrackNeg.tpcFoundOverFindableCls(),
                   theTrackNeg.tpcCrossedRowsOverFindableCls());
  /*
  fFuncConversionData(theV0.pt(),
                      theV0.eta(),
                      theV0.phi(),
                      theV0.v0radius(),
                      lV0CosinePA,
                      theV0.alpha(),
                      theV0.qtarm(),
                      theV0.psipair());
  */
  fillReconstructedInfoHistograms("afterCuts/",
                                  theV0,
                                  theTrackPos,
                                  theTrackNeg,
                                  lV0CosinePA);
  return kTRUE;
}

template <typename TV0, typename TTRACK, typename TMC>
void processTruePhoton(std::string theBAC, const TV0& theV0, const TTRACK& theTrackPos, const TTRACK& theTrackNeg, const TMC& theMcParticles)
{
  // todo: verify it is enough to check only mother0 being equal

  /* example from https://aliceo2group.github.io/analysis-framework/docs/tutorials/indexTables.html?highlight=_as:
   * track0.collision_as<myCol>().mult() : access multiplicity of collission associated with track0
   */

  //~ // use & here?
  auto lMcPos = theTrackPos.template mcParticle_as<aod::McParticles_001>();
  auto lMcNeg = theTrackNeg.template mcParticle_as<aod::McParticles_001>();

}

  void processMC(aod::Collisions::iterator  const& theCollision,
                 aod::V0Datas const& theV0s,
                 tracksAndTPCInfoMC const& theTracks,
                 aod::McParticles const& theMcParticles)
  {
    for (auto& lV0 : theV0s) {

      auto lTrackPos = lV0.template posTrack_as<tracksAndTPCInfoMC>(); // positive daughter
      auto lTrackNeg = lV0.template negTrack_as<tracksAndTPCInfoMC>(); // negative daughter


      if (!processPhoton(theCollision, lV0, lTrackPos, lTrackNeg)) {
        continue;
      }
      //~ processTruePhoton("afterCuts/",
                        //~ lV0,
                        //~ lTrackPos,
                        //~ lTrackNeg,
                        //~ theMcParticles);
    }
  }

  PROCESS_SWITCH(SkimmerMc, processMC, "Process MC", true);

  template <typename T>
  std::shared_ptr<T> getTH(std::map<std::string, HistPtr> const& theMap, std::string const& theName)
  {
    auto lPairIt = theMap.find(theName);
    if (lPairIt == theMap.end()) {
      LOGF(warning, "SFS getTH(): No element with key %s in map.", theName.data());
      return std::shared_ptr<T>(); // SFS todo verify this is correct
    }
    if (!std::holds_alternative<std::shared_ptr<T>>(lPairIt->second)) {
      LOGF(warning, "SFS getTH(): No shared_ptr<T> in map for key %s.", theName.data());
      return std::shared_ptr<T>();
    }
    std::shared_ptr<T> lHisto = std::get<std::shared_ptr<T>>(lPairIt->second);
    if (lHisto == nullptr) {
      LOGF(warning, "SFS getTH(): Found shared_ptr<T> for key %s but it is a nullptr.", theName.data());
      return std::shared_ptr<T>();
    }
    return lHisto;
  }

  void fillTH1(std::map<std::string, HistPtr> const& theMap, std::string const& theName, float theValue)
  {
    std::shared_ptr<TH1> lHisto = getTH<TH1>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValue);
    }
  }

  void fillTH2(std::map<std::string, HistPtr> const& theMap, std::string const& theName, float theValueX, float theValueY)
  {
    std::shared_ptr<TH2> lHisto = getTH<TH2>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValueX, theValueY);
    }
  }

  template <typename TTRACK>
  void fillTrackHistograms(std::string theHistoPath, const TTRACK& theTrackPos, const TTRACK& theTrackNeg)
  {
    auto fillTrackHistogramsI = [&](const TTRACK& theTrack) {
      fillTH1(fTrackHistos, theHistoPath + "hTPCFoundOverFindableCls", theTrack.tpcFoundOverFindableCls());
      fillTH1(fTrackHistos, theHistoPath + "hTPCCrossedRowsOverFindableCls", theTrack.tpcCrossedRowsOverFindableCls());
      fillTH2(fTrackHistos, theHistoPath + "hTPCdEdxSigEl", theTrack.p(), theTrack.tpcNSigmaEl());
      fillTH2(fTrackHistos, theHistoPath + "hTPCdEdxSigPi", theTrack.p(), theTrack.tpcNSigmaPi());
      fillTH2(fTrackHistos, theHistoPath + "hTPCdEdx", theTrack.p(), theTrack.tpcSignal());
    };
    fillTrackHistogramsI(theTrackPos);
    fillTrackHistogramsI(theTrackNeg);
  }

  template <typename TV0>
  void fillV0Histograms(std::string theHistoPath, const TV0& theV0, float theV0CosinePA)
  {
    fillTH1(fV0Histos, theHistoPath + "hPtRec", theV0.pt());
    fillTH1(fV0Histos, theHistoPath + "hEtaRec", theV0.eta());
    fillTH1(fV0Histos, theHistoPath + "hPhiRec", theV0.phi());
    fillTH1(fV0Histos, theHistoPath + "hConvPointR", theV0.v0radius());
    fillTH1(fV0Histos, theHistoPath + "hCosPAngle", theV0CosinePA);
    fillTH2(fV0Histos, theHistoPath + "hArmenteros", theV0.alpha(), theV0.qtarm());
    fillTH2(fV0Histos, theHistoPath + "hPsiPtRec", theV0.psipair(), theV0.pt());
  }

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // single track eta cut
    if (TMath::Abs(theTrack.eta()) > fEtaCut) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTrackEta"));
      return kFALSE;
    }

    // single track pt cut
    if (theTrack.pt() < fSinglePtCut) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTrackPt"));
      return kFALSE;
    }

    if (!(selectionPIDTPC_track(theTrack))) {
      return kFALSE;
    }

    if (theTrack.tpcFoundOverFindableCls() < fMinTPCFoundOverFindableCls) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTPCFoundOverFindableCls"));
      return kFALSE;
    }

    if (theTrack.tpcCrossedRowsOverFindableCls() < fMinTPCCrossedRowsOverFindableCls) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTPCCrossedRowsOverFindableCls"));
      return kFALSE;
    }

    return kTRUE;
  }

  template <typename T>
  bool passesPhotonCuts(const T& theV0, float theV0CosinePA)
  {
    if (theV0.v0radius() < fMinR || theV0.v0radius() > fMaxR) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kV0Radius"));
      return kFALSE;
    }

    if (!ArmenterosQtCut(theV0.alpha(), theV0.qtarm(), theV0.pt())) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kArmenteros"));
      return kFALSE;
    }

    if (TMath::Abs(theV0.psipair()) > fPsiPairCut) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPsiPair"));
      return kFALSE;
    }

    if (theV0CosinePA < fCosPAngleCut) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kCosinePA"));
      return kFALSE;
    }

    return kTRUE;
  }

  template <typename TV0, typename TTRACK>
  void fillReconstructedInfoHistograms(std::string theBAC, const TV0& theV0, const TTRACK& theTrackPos, const TTRACK& theTrackNeg, float theV0CosinePA)
  {
    fillTrackHistograms(fPrefixReconstructedInfoHistos + "track/" + theBAC, theTrackPos, theTrackNeg);
    fillV0Histograms(fPrefixReconstructedInfoHistos + "v0/" + theBAC, theV0, theV0CosinePA);

    if (theBAC == "beforeCuts/") {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPhotonIn"));
    } else if (theBAC == "afterCuts/") {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPhotonOut"));
    }
  }

  Bool_t ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt)
  {
    // in AliPhysics this is the cut for if fDo2DQt && fDoQtGammaSelection == 2
    Float_t lQtMaxPtDep = fQtPtMax * thePt;
    if (lQtMaxPtDep > fQtMax) {
      lQtMaxPtDep = fQtMax;
    }
    if (!(TMath::Power(theAlpha / fMaxPhotonAsymmetry, 2) + TMath::Power(theQt / lQtMaxPtDep, 2) < 1)) {
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename T>
  bool selectionPIDTPC_track(const T& theTrack)
  {
    // TPC Electron Line
    if (theTrack.tpcNSigmaEl() < fPIDnSigmaBelowElectronLine || theTrack.tpcNSigmaEl() > fPIDnSigmaAboveElectronLine) {
      fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kElectronPID"));
      return kFALSE;
    }

    // TPC Pion Line
    if (theTrack.p() > fPIDMinPnSigmaAbovePionLine) {
      // low pt Pion rej
      if (theTrack.p() < fPIDMaxPnSigmaAbovePionLine) {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLine) {
          fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPionRejLowMom"));
          return kFALSE;
        }
      }
      // High Pt Pion rej
      else {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineHighP) {
          fillTH1(fV0Histos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPionRejHighMom"));
          return kFALSE;
        }
      }
    }
    return kTRUE;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SkimmerMc>(cfgc)};
}
