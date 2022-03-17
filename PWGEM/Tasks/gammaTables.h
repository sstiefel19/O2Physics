#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/StrangenessTables.h"

#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTPC.h"

//todo: declare more columns in this file dynamic or expression, atm I save a lot of redundant information
namespace o2::aod
{
namespace gammatrackreco
{    
DECLARE_SOA_COLUMN(PositiveCharge, positiveCharge, bool);    
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
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
                  gammatrackreco::Eta,
                  gammatrackreco::P,
                  gammatrackreco::Phi,
                  gammatrackreco::Pt,
                  track::TPCSignal,
                  pidtpc::TPCNSigmaEl,
                  pidtpc::TPCNSigmaPi);

// need to declare this again(these columns also exist in v0data::) since some of them are dynamical there
namespace gammamctrue
{
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
// todo: declare those expression or dynamic columns from Px,Py, Pz
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
} // gammamctrue

DECLARE_SOA_TABLE(GammaConversionsInfoTrue, "AOD", "V0INFOTRUE",
                  o2::soa::Index<>,
                  v0data::V0Id,
                  // obtained from daughter tracks
                  v0data::X, v0data::Y, v0data::Z, // origin positive daughter
                  v0data::PxPos, v0data::PyPos, v0data::PzPos,
                  v0data::PxNeg, v0data::PyNeg, v0data::PzNeg,
                  // true info from mother
                  gammamctrue::Px, gammamctrue::Py, gammamctrue::Pz,
                  gammamctrue::Eta,
                  gammamctrue::P,
                  gammamctrue::Phi,
                  gammamctrue::Pt,

                  // Dynamic columns
                  v0data::V0Radius<v0data::X, v0data::Y>,
                  v0data::DistOverTotMom<v0data::X, v0data::Y, v0data::Z, gammamctrue::Px, gammamctrue::Py, gammamctrue::Pz>,
                  v0data::V0CosPA<       v0data::X, v0data::Y, v0data::Z, gammamctrue::Px, gammamctrue::Py, gammamctrue::Pz>,
                  v0data::DCAV0ToPV<     v0data::X, v0data::Y, v0data::Z, gammamctrue::Px, gammamctrue::Py, gammamctrue::Pz>,
                  v0data::Alpha<  v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::QtArm<  v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::PsiPair<v0data::PxPos, v0data::PyPos, v0data::PzPos, v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,

                  // Longitudinal
                  v0data::YK0Short<gammamctrue::Px, gammamctrue::Py, gammamctrue::Pz>,
                  v0data::YLambda<gammamctrue::Px, gammamctrue::Py, gammamctrue::Pz>,
                  v0data::NegativePt<v0data::PxNeg, v0data::PyNeg>,
                  v0data::PositivePt<v0data::PxPos, v0data::PyPos>,
                  v0data::NegativeEta<v0data::PxNeg, v0data::PyNeg, v0data::PzNeg>,
                  v0data::NegativePhi<v0data::PxNeg, v0data::PyNeg>,
                  v0data::PositiveEta<v0data::PxPos, v0data::PyPos, v0data::PzPos>,
                  v0data::PositivePhi<v0data::PxPos, v0data::PyPos>);
} // o2::aod
