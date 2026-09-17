#pragma once

#include <vector>

#include "Config.h"
#include "Context.h"

namespace urwell {

struct AvalancheResults {
  int neTotal = 0;
  int niTotal = 0;
  int nIonsDrifted = 0;
  int nIonsFailed = 0;
};

// Creates ctx.avalanche and ctx.driftIon (and ctx.viewDrift if
// kStoreDriftLines is set) and wires up signal calculation on both.
// Call once, after SetupDetector and before RunAvalancheLoop.
void SetupAvalanche(SimContext& ctx);

// Avalanches every seed electron, drifts the resulting ions, and prints
// a rate-limited progress bar as it goes. Returns the running totals
// and prints the same summary lines the original monolithic loop did.
AvalancheResults RunAvalancheLoop(SimContext& ctx,
                                   const std::vector<ElectronSeed>& seeds);

}  // namespace urwell
