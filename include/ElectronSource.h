#pragma once

#include <vector>

#include "Config.h"
#include "Context.h"

namespace urwell {

// Returns the primary electron seeds to avalanche: either TrackHeed
// clusters from an alpha track (kUseTrackHeed == true) or the manually
// specified kManualElectrons list. Requires ctx.sensor to already be set
// up (see SetupDetector).
std::vector<ElectronSeed> GenerateElectronSeeds(SimContext& ctx);

}  // namespace urwell
