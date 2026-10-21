#pragma once

#include <string>
#include <vector>

namespace urwell {

// ---------------------------------------------------------------------
// Global run mode switches.
//   kComputePlots = false -> skip all ROOT plotting and all plotting-based
//   diagnostics. This is the coarse master switch.
//   kComputeGainHistogram = true -> draw the separate amplification
//   histogram window.
//   kComputeDetectorPlots = true -> draw the mesh/signal/field plots.
//   kStoreDriftLines is only respected when the detector plot path is enabled.
//   kEnableIonSimulation defaults to true, but can be overridden from the
//   command line with --no-ion-simulation to keep only the electron
//   avalanche path.
// ---------------------------------------------------------------------
inline constexpr bool kComputePlots = true;
inline constexpr bool kComputeGainHistogram = true;
inline constexpr bool kComputeDetectorPlots = true;
inline constexpr bool kComputeAmplificationDistribution = true;
inline constexpr bool kStoreDriftLines = true;  // EnableDriftLines() on avalanche + ion drift
inline constexpr bool kShowMeshPlot    = true;  // draw cMesh (geometry + drift lines)
inline bool kEnableIonSimulation = true;

// ---------------------------------------------------------------------
// Manual electron generation for gain studies.
// If kUseTrackHeed is false, the code generates kManualElectronCount
// identical or evenly spaced seed electrons around a reference point.
// This is convenient for single-point gain distributions, z-scans, or
// repeated avalanche statistics without needing a full TrackHeed track.
// ---------------------------------------------------------------------
inline constexpr int kManualElectronCount = 1;
inline constexpr double kManualElectronX = 0.;
inline constexpr double kManualElectronY = 0.;
inline constexpr double kManualElectronZ0 = 0.01;
inline constexpr double kManualElectronDeltaZ = 0.;
inline constexpr double kManualElectronT = 0.;

// Ion drift-line integration step size [cm]. Finer = more accurate but
// slower; only matters when kStoreDriftLines is true, but also affects
// the granularity/cost of the ion transport integration itself.
inline constexpr double kIonStepSize = 2.e-3;

// z [cm] a primary electron must get below to count as "entered the well".
// Set it just BELOW the top surface of the top electrode (read it off the
// mesh plot). 0.0 is a placeholder.
inline constexpr double kWellEntryZ = 2.e-10;

// ---------------------------------------------------------------------
// Hybrid transport: AvalancheMC (tabulated transport) in the drift gap,
// AvalancheMicroscopic from kZSwitch downwards (well / amplification).
// kZSwitch [cm] must sit ABOVE the strong-field region, where the field
// is still smooth. Read it off the cField plot. 0.02 is a placeholder.
// kDriftStepSize [cm] must be much smaller than (kZSwitch - top of well).
// ---------------------------------------------------------------------
inline constexpr bool kUseHybridDrift = true;
inline constexpr double kZSwitch = 0.01;
inline constexpr double kDriftStepSize = 2.e-3;
inline const std::string kGasTableFile = "/Users/alexandrosgardikioti/garfieldpp/urwell_v0_1/gas_files/settings_0/two_gas_mixtures/Ar_CO2_mixtures/Ar_CO2_80_20_NTP.gas";

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

// ---------------------------------------------------------------------
// GEANT4-hits electron source switch.
//   true  -> skip TrackHeed and the manual electron list; read
//            ionization clusters from a GEANT4 ROOT "hits" file (see
//            geant_hits.h / LoadGeantHitsByTrack) and avalanche those
//            instead. Takes priority over kUseTrackHeed if both are true.
//   false -> fall back to the existing kUseTrackHeed / manual behavior.
// ---------------------------------------------------------------------
inline constexpr bool kUseGeant4Hits = true;

// Path to the GEANT4 output ROOT file. Must contain a TTree named
// "hits" with the branches geant_hits.h expects: event_id, track_id,
// edep, x/y/z_pos_pre, x/y/z_pos_post, t_pre, t_post, and one of
// pdg_id / particle_pdg_id / particle_pdg / PDG.
inline const std::string kGeant4HitsFile = "/Users/alexandrosgardikioti/garfieldpp/urwell_v0_1/geant4root/output0.root";

// Which GEANT4 event to read; -1 reads every event in the file and
// concatenates all their clusters into one flat seed list.
inline constexpr int kGeant4EventId = -1;

// Max electron-ion pairs represented by a single ElectronSeed. Larger
// = fewer, "fatter" seeds (faster, coarser); must be >= 1.
inline constexpr int kGeant4ClusterSize = 1;

// x, y, z [cm], t [ns] for each manually seeded electron. Only used
// when kUseTrackHeed is false. Add/remove entries as needed, e.g. to
// scan across the gap or repeat the same point N times for statistics.
struct ElectronSeed { double x, y, z, t; int nElectrons = 1; };

inline const std::vector<ElectronSeed> kManualElectrons = {
  {kManualElectronX, kManualElectronY, kManualElectronZ0, kManualElectronT},
  {kManualElectronX, kManualElectronY, 0.10, kManualElectronT},
  {kManualElectronX, kManualElectronY, 0.03, kManualElectronT},
};

}  // namespace urwell
