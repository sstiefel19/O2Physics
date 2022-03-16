#include "Framework/AnalysisDataModel.h"

#include "Common/DataModel/StrangenessTables.h"

#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/PIDTPC.h"

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
} // o2::aod
