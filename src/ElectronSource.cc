#include "ElectronSource.h"

#include "Garfield/TrackHeed.hh"

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

std::vector<ElectronSeed> GenerateElectronSeeds(SimContext& ctx) {
  std::vector<ElectronSeed> seeds;

  if (kUseTrackHeed) {
    TrackHeed* track = new TrackHeed();
    track->SetParticle("alpha");
    track->SetKineticEnergy(1e6);  // eV -- e.g. 5.5 MeV, adjust to your source
    track->SetSensor(ctx.sensor);

    // Starting point and direction of the track [cm], [ns]
    double tx0 = 0., ty0 = 0., tz0 = 0.16, tt0 = 0.;
    double dx = 0., dy = 0., dz = -1.;  // direction vector, need not be normalized
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
