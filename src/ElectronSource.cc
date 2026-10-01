#include "ElectronSource.h"

#include "Garfield/TrackHeed.hh"
#include "geant_hits.h"

#include <algorithm>
#include <iostream>
#include <vector>

using namespace Garfield;

namespace urwell {

static std::vector<ElectronSeed> GenerateManualSeeds() {
  std::vector<ElectronSeed> seeds;
  const int nSeeds = std::max(1, kManualElectronCount);
  seeds.reserve(static_cast<size_t>(nSeeds));

  for (int i = 0; i < nSeeds; ++i) {
    const double z = kManualElectronZ0 + i * kManualElectronDeltaZ;
    seeds.push_back({kManualElectronX, kManualElectronY, z, kManualElectronT});
  }

  return seeds;
}

static std::vector<ElectronSeed> GenerateGeant4Seeds(SimContext& ctx) {
  std::vector<ElectronSeed> seeds;

  // W-value must come from the already-Initialise()'d gas. SetupDetector()
  // runs before GenerateElectronSeeds() in main.cc, so this is safe here.
  const double wValue = ctx.gas->GetW();
  if (wValue <= 0.) {
    std::cerr << "GenerateGeant4Seeds: gas->GetW() returned " << wValue
              << " -- has the gas been Initialise()'d?\n";
    return seeds;
  }

  constexpr double kXOffsetMm = 100.;
  constexpr double kYOffsetMm = 150.;
  constexpr double kZOffsetMm = 1.55;

  const auto tracks = LoadGeantHitsByTrack(kGeant4HitsFile, kGeant4EventId,
                                            wValue, kGeant4ClusterSize,
                                            kXOffsetMm, kYOffsetMm, kZOffsetMm);

  for (const auto& track : tracks) {
    for (const auto& p : track.points) {
      seeds.push_back({p.x, p.y, p.z, p.t, p.n_electrons});
    }
  }

  std::cout << "GEANT4 hits: " << tracks.size() << " track(s) -> "
            << seeds.size() << " seed cluster(s) (<= " << kGeant4ClusterSize
            << " electrons each).\n";
  return seeds;
}

std::vector<ElectronSeed> GenerateElectronSeeds(SimContext& ctx) {
  std::vector<ElectronSeed> seeds;

  if (kUseGeant4Hits) {
    seeds = GenerateGeant4Seeds(ctx);
  } else if (kUseTrackHeed) {
    TrackHeed* track = new TrackHeed();
    track->SetParticle("alpha");
    track->SetKineticEnergy(1e6);
    track->SetSensor(ctx.sensor);

    double tx0 = 0., ty0 = 0., tz0 = 0.16, tt0 = 0.;
    double dx = 0., dy = 0., dz = -1.;
    track->NewTrack(tx0, ty0, tz0, tt0, dx, dy, dz);

    for (const auto& cluster : track->GetClusters()) {
      for (const auto& electron : cluster.electrons) {
        seeds.push_back({electron.x, electron.y, electron.z, electron.t});
      }
    }
    std::cout << "TrackHeed generated " << seeds.size() << " primary electrons.\n";
  } else {
    seeds = GenerateManualSeeds();
    std::cout << "Using " << seeds.size()
              << " manually seeded electron(s), TrackHeed skipped.\n";
  }

  return seeds;
}

}  // namespace urwell