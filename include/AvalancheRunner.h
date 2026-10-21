#pragma once

#include <vector>

#include "Config.h"
#include "Context.h"

namespace urwell {

// One entry per avalanched electron: where it was seeded, where the microscopic
// avalanche actually started (after the hybrid drift), its gain, and whether it
// entered the well. Positions in cm.
struct GainRecord {
  double xSeed, ySeed, zSeed;
  double xStart, yStart, zStart;
  int gain;
  int entered;
};

// Where the electrons of one avalanche ended up. Per avalanche (not per electron), so the file stays small:
//   endZ[b]      = number of electron endpoints with z1 in z-bin b of [kEndZMin, kEndZMax] (cm)
//   endStatus[s] = number of electron endpoints with |status| == s (s clamped to kEndStatusBins-1)
inline constexpr int kEndZBins = 30;
inline constexpr double kEndZMin = -0.02;   // cm (bottom of the sensor area)
inline constexpr double kEndZMax = 0.01;    // cm (= kZSwitch)
inline constexpr int kEndStatusBins = 16;

struct EndpointRecord {
  int gain;
  int entered;
  short endZ[kEndZBins];
  short endStatus[kEndStatusBins];
};

struct AvalancheResults {
  std::vector<EndpointRecord> endRecords;  // per-avalanche endpoint summary (written to endpointTree)
  std::vector<GainRecord> records;  // per-electron position + gain (written to gainPosTree)
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
