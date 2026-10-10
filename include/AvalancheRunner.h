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
  int nOutsideMap = 0;    // start point not in a valid drift region (hybrid mode)
  int nLostInDrift = 0;   // attached / never reached kZSwitch (hybrid mode)
  std::vector<int> gains;
  std::vector<int> gainsRaw;  // avalanche size, only for electrons that entered the well
  int nTried = 0;             // electrons actually avalanched (not lost in drift / invalid medium)
  int nEntered = 0;           // of those, how many reached below kWellEntryZ
};

// Creates ctx.avalanche and, when kEnableIonSimulation is true, also
// ctx.driftIon (and ctx.viewDrift if kStoreDriftLines is set). The ion
// drift path can be disabled at runtime with --no-ion-simulation to keep
// only the electron avalanche path.
// Call once, after SetupDetector and before RunAvalancheLoop.
void SetupAvalanche(SimContext& ctx);

// Avalanches every seed electron, drifts the resulting ions, and prints
// a rate-limited progress bar as it goes. Returns the running totals
// and prints the same summary lines the original monolithic loop did.
AvalancheResults RunAvalancheLoop(SimContext& ctx,
                                   const std::vector<ElectronSeed>& seeds);

}  // namespace urwell
