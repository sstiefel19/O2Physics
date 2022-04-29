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

#include "Framework/AnalysisDataModel.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

// todo: declare more columns in this file dynamic or expression, atm I save a lot of redundant information
namespace o2::aod
{
namespace gammatrackreco
{
DECLARE_SOA_COLUMN(IsFromConversionPhoton, isFromConversionPhoton, bool);                //! Whether this track is from a MC confirmed photon 
DECLARE_SOA_COLUMN(Eta, eta, float);                                                     //! Pseudorapidity
DECLARE_SOA_COLUMN(P, p, float);                                                         //! Total momentum in GeV/c
DECLARE_SOA_COLUMN(Phi, phi, float);                                                     //! Azimuthal angle
DECLARE_SOA_COLUMN(Pt, pt, float);                                                       //! Transversal momentum in GeV/c
DECLARE_SOA_COLUMN(PositivelyCharged, positivelyCharged, bool);                                //! True for positively charged track
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float); //! Ratio  crossed rows over findable clusters
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, tpcFoundOverFindableCls, float);             //! Ratio of found over findable clusters
DECLARE_SOA_COLUMN(TpcNClsCrossedRows, tpcNClsCrossedRows, float);                       //! Number of crossed TPC Rows
} // namespace gammatrackreco

DECLARE_SOA_TABLE(GammaConversionTracks, "AOD", "V0TRACKS",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  gammatrackreco::IsFromConversionPhoton,
                  track::DcaXY,
                  gammatrackreco::Eta,
                  gammatrackreco::P,
                  gammatrackreco::Phi,
                  gammatrackreco::Pt,
                  gammatrackreco::PositivelyCharged,
                  gammatrackreco::TpcCrossedRowsOverFindableCls,
                  gammatrackreco::TpcFoundOverFindableCls,
                  gammatrackreco::TpcNClsCrossedRows,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaPi,
                  track::TPCSignal);

namespace gammamctrue
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);    //! Collision to which this McGamma belongs
DECLARE_SOA_COLUMN(Gamma, gamma, int64_t);    //! just a number that can be used to associate daughter particles with a gamma
DECLARE_SOA_COLUMN(NDaughters, nDaughters, int); // SFS use unsigned!
//~ DECLARE_SOA_COLUMN(ConversionRadius, conversionRadius, float);  // !
//~ DECLARE_SOA_COLUMN(Alpha, alpha, float);  // !
//~ DECLARE_SOA_COLUMN(QtArm, qtArm, float);  // !
//~ DECLARE_SOA_COLUMN(PsiPair, psiPair, float);  // !
}

DECLARE_SOA_TABLE(StoredMcGammasInfoTrue, "AOD", "MCGAINFOTRUE",
                  o2::soa::Index<>,
                  mcparticle::McCollisionId,
                  gammamctrue::Gamma,
                  v0data::V0Id, // reference to reconstructed v0
                  mcparticle::StatusCode, 
                  mcparticle::Flags,
                  mcparticle::Px, mcparticle::Py, mcparticle::Pz,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz, mcparticle::Vt,
                  gammamctrue::NDaughters,
                  
                  // Dynamic columns
                  mcparticle::ProducedByGenerator<mcparticle::Flags>,
                  mcparticle::FromBackgroundEvent<mcparticle::Flags>,
                  mcparticle::GetGenStatusCode<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::GetProcess<mcparticle::Flags, mcparticle::StatusCode>,
                  mcparticle::IsPhysicalPrimary<mcparticle::Flags>);
                  
DECLARE_SOA_EXTENDED_TABLE(McGammasInfoTrue, StoredMcGammasInfoTrue, "MCGAINFOTRUE", //! Basic MC particle properties
                           mcparticle::Phi,
                           mcparticle::Eta,
                           mcparticle::Pt,
                           mcparticle::P,
                           mcparticle::Y);

namespace gammamctrue
{
DECLARE_SOA_INDEX_COLUMN_FULL(Mother, mother, int, McGammasInfoTrue, ""); // SFS do I need to worry of overflows?
DECLARE_SOA_COLUMN(NMothers, nMothers, int);
DECLARE_SOA_DYNAMIC_COLUMN(R, r,   //! vertex 2d r
      [](float x, float y) -> float { return RecoDecay::sqrtSumOfSquares(x, y); });
}

// table to hold daughter particles of MC gammas
DECLARE_SOA_TABLE(MCGammaDaughters, "AOD", "MCGADAUGHTERS",
                  o2::soa::Index<>,
                  mcparticle::McCollisionId, // SFS might be superflous since there is already a pointer to MCGammas which point to mccollision
                  gammamctrue::MotherId,
                  gammamctrue::NMothers,
                  mcparticle::PdgCode,
                  mcparticle::Vx, mcparticle::Vy, mcparticle::Vz,
                  mcparticle::Eta,
                  mcparticle::P,
                  mcparticle::Phi,
                  mcparticle::Pt,
                  gammamctrue::R<mcparticle::Vx, mcparticle::Vy>); //                       v0data::V0Radius<v0data::X, v0data::Y>,
} // namespace o2::aod
