#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/StrangenessTables.h"

#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTPC.h"

namespace o2::aod
{
namespace gammatrackreco
{    
DECLARE_SOA_COLUMN(PositiveCharge, positiveCharge, bool);    
DECLARE_SOA_COLUMN(TpcFoundOverFindableCls, tpcFoundOverFindableCls, float);
DECLARE_SOA_COLUMN(TpcCrossedRowsOverFindableCls, tpcCrossedRowsOverFindableCls, float);
DECLARE_SOA_COLUMN(IsFromConversionPhoton, isFromConversionPhoton, bool);
}

DECLARE_SOA_TABLE(GammaConversionTracks, "AOD", "V0TRACKS",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  gammatrackreco::PositiveCharge,
                  gammatrackreco::IsFromConversionPhoton,
                  gammatrackreco::TpcFoundOverFindableCls,
                  gammatrackreco::TpcCrossedRowsOverFindableCls,
                  track::Eta,
                  track::P,
                  track::Pt,
                  track::TPCSignal,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaPi);

namespace gammamctrue
{
// using suffices C here for computed from daughter tracks
//~ DECLARE_SOA_COLUMN(EtaC, etaC, float);
//~ DECLARE_SOA_COLUMN(PhiC, phiC, float);
//~ DECLARE_SOA_COLUMN(PtC, ptC, float);
} // gammamctrue

//todo add true track info table?
// SFS todo: need to add some sort of indexing here to get from recos to true table
DECLARE_SOA_TABLE(GammaConversionsInfoTrue, "AOD", "V0INFOTRUE",
                  o2::soa::Index<>, v0data::V0Id,
                  v0data::X,     v0data::Y,     v0data::Z,
                  v0data::PxPos, v0data::PyPos, v0data::PzPos,
                  v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,

                  // true info from mother
                  track::Eta,
                  track::P,
                  track::Phi,
                  track::Pt,
                  
                  // derived info from true daughter track info
                  // Dynamic columns
                  v0data::V0Radius<v0data::X, v0data::Y>,
                  v0data::DistOverTotMom<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                  v0data::V0CosPA<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                  v0data::DCAV0ToPV<v0data::X, v0data::Y, v0data::Z, v0data::Px, v0data::Py, v0data::Pz>,
                  v0data::Alpha<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::QtArm<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::PsiPair<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  
                  // Longitudinal
                  v0data::YK0Short<v0data::Px, v0data::Py, v0data::Pz>,
                  v0data::YLambda<v0data::Px, v0data::Py, v0data::Pz>,
                  v0data::NegativePt<v0data::PxNeg, v0data::PyNeg>,
                  v0data::PositivePt<v0data::PxPos, v0data::PyPos>,
                  v0data::NegativeEta<v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::NegativePhi<v0data::PxNeg, v0data::PyNeg>,
                  v0data::PositiveEta<v0data::PxPos, v0data::PyPos, v0data::PzPos>,
                  v0data::PositivePhi<v0data::PxPos, v0data::PyPos>);
                  
                  //~ // changed the name here so the names originally defined for StoredV0Datas can be used for the true mother info
                  //~ v0data::EtaC<v0data::Px, v0data::Py, v0data::Pz>,
                  //~ v0data::PhiC<v0data::Px, v0data::Py>,
                  //~ v0data::PtC<v0data::PxPos, v0data::PyPos, v0data::PxNeg, v0data::PyNeg>);
} // o2::aod
