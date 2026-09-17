#pragma once

#include <vector>

namespace urwell {

// ---------------------------------------------------------------------
// Performance switches. Flip these off for fast/batch runs where you
// only care about totals and the induced signal, not the drift-line
// picture. Storing/plotting every electron and ion trajectory is the
// main memory/time cost of EnableDriftLines()+EnablePlotting(); the
// mesh canvas is only meaningful when drift lines are stored, so it's
// gated by the same flag.
// ---------------------------------------------------------------------
inline constexpr bool kStoreDriftLines = true;  // EnableDriftLines() on avalanche + ion drift
inline constexpr bool kShowMeshPlot    = true;  // draw cMesh (geometry + drift lines)

// Ion drift-line integration step size [cm]. Finer = more accurate but
// slower; only matters when kStoreDriftLines is true, but also affects
// the granularity/cost of the ion transport integration itself.
inline constexpr double kIonStepSize = 2.e-3;

// ---------------------------------------------------------------------
// Weighting field switch.
//   true  -> load the time-dependent (dynamic) weighting potential from
//            COMSOL via SetDynamicWeightingPotential(), and enable the
//            Sensor's delayed-signal machinery on top of it. This is
//            what's needed to see the resistive-layer/DLC charge
//            spreading contribution, but it is the more expensive path
//            (every drift/avalanche step now also touches the delayed
//            signal), especially for large TrackHeed avalanches.
//   false -> load the static weighting field via SetWeightingField()
//            instead -- prompt-only signal, no delayed component, much
//            cheaper. Good for fast gain/geometry debugging runs where
//            the resistive spreading isn't what you're checking.
// All delayed-signal setup, diagnostics, and plots are gated on this
// single flag, so flipping it is the only thing you need to do.
// ---------------------------------------------------------------------
inline constexpr bool kUseDynamicWeightingField = false;

// ---------------------------------------------------------------------
// Number of time points used for the DELAYED signal calculation only.
// This is independent of nSteps (the prompt-signal bin count, held in
// SimContext::nSteps). Every drift/avalanche step's delayed contribution
// gets distributed across this many time points, so it's the dominant
// cost multiplier for large avalanches (e.g. TrackHeed alpha tracks with
// millions of electrons/ions). We've confirmed the delayed weighting
// potential only varies meaningfully on a ~10-50 ns scale and plateaus
// by ~100 ns, so it doesn't need fine (0.1 ns) binning -- a coarser grid
// here costs essentially no physical accuracy while cutting runtime by
// roughly nSteps / kNumDelayedSteps. Raise this back up if you ever need
// finer time resolution specifically on the delayed/resistive component.
// Only used when kUseDynamicWeightingField is true.
// ---------------------------------------------------------------------
inline constexpr int kNumDelayedSteps = 4;  // was implicitly 2000 (== nSteps)

// ---------------------------------------------------------------------
// Electron source switch.
//   true  -> use TrackHeed to generate a realistic alpha-particle track
//            and take its ionization clusters as the starting electrons.
//   false -> skip TrackHeed entirely and avalanche a manually specified
//            list of electron starting points instead (see
//            kManualElectrons below). Useful for single-electron gain
//            studies, scans over starting position, or fast debugging
//            without paying the cost of TrackHeed cluster generation.
// ---------------------------------------------------------------------
inline constexpr bool kUseTrackHeed = false;

// Print at least every kProgressIntervalSeconds of wall-clock time even
// if the electron-count trigger hasn't fired yet -- important because
// individual avalanches near the well can be expensive (Magboltz
// rebuilding its collision-rate table repeatedly), so a count-based
// trigger alone could mean a long silent stretch with no feedback.
inline constexpr double kProgressIntervalSeconds = 30.;

// x, y, z [cm], t [ns] for each manually seeded electron. Only used
// when kUseTrackHeed is false. Add/remove entries as needed, e.g. to
// scan across the gap or repeat the same point N times for statistics.
struct ElectronSeed { double x, y, z, t; };

inline const std::vector<ElectronSeed> kManualElectrons = {
  {0., 0., 0.16, 0.},
  {0., 0., 0.10, 0.},
  {0., 0., 0.03, 0.},
};

}  // namespace urwell
