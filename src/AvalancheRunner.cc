#include "AvalancheRunner.h"
#include "ProgressBar.h"

#include <algorithm>
#include <ctime>
#include <iostream>

using namespace Garfield;

namespace urwell {

void SetupAvalanche(SimContext& ctx) {
  const bool shouldStoreDriftLines = kComputePlots && kComputeDetectorPlots && kStoreDriftLines;

  AvalancheMicroscopic* avalanche = new AvalancheMicroscopic();
  avalanche->SetSensor(ctx.sensor);
  avalanche->EnableSignalCalculation();
  if (shouldStoreDriftLines) avalanche->EnableDriftLines();

  ViewDrift* viewDrift = nullptr;
  if (shouldStoreDriftLines) {
    viewDrift = new ViewDrift();
    avalanche->EnablePlotting(viewDrift);
  }

  AvalancheMC* driftIon = new AvalancheMC();
  driftIon->SetSensor(ctx.sensor);
  driftIon->SetDistanceSteps(kIonStepSize);
  driftIon->EnableSignalCalculation();
  if (shouldStoreDriftLines) {
    driftIon->EnableDriftLines();
    driftIon->EnablePlotting(viewDrift);  // reuse the SAME ViewDrift -> ions show on same plot
  }

  AvalancheMC* driftEl = nullptr;
  if (kUseHybridDrift) {
    driftEl = new AvalancheMC();
    driftEl->SetSensor(ctx.driftSensor);
    driftEl->SetDistanceSteps(kDriftStepSize);
    driftEl->EnableSignalCalculation();
    if (shouldStoreDriftLines) {
      driftEl->EnableDriftLines();
      driftEl->EnablePlotting(viewDrift);
    }
  }
  ctx.avalanche = avalanche;
  ctx.driftElectron = driftEl;
  ctx.driftIon = driftIon;
  ctx.viewDrift = viewDrift;
}

AvalancheResults RunAvalancheLoop(SimContext& ctx,
                                   const std::vector<ElectronSeed>& seeds) {
  AvalancheResults results;
  int nSkippedInvalidMedium = 0;

  auto processElectron = [&](double x, double y, double z, double t, int nElectrons) {
    for (int k = 0; k < nElectrons; ++k) {
      double xs = x, ys = y, zs = z, ts = t;

      if (kUseHybridDrift && z > kZSwitch) {
        const bool ok = ctx.driftElectron->DriftElectron(x, y, z, t);
        if (!ok || ctx.driftElectron->GetNumberOfElectronEndpoints() == 0) {
          ++results.nOutsideMap;   // start point not in a valid drift region
          continue;
        }
        double x0, y0, z0, t0, x1, y1, z1, t1;
        int st = 0;
        ctx.driftElectron->GetElectronEndpoint(0, x0, y0, z0, t0, x1, y1, z1, t1, st);
        if (z1 > kZSwitch + 2. * kDriftStepSize) { ++results.nLostInDrift; continue; }
        xs = x1; ys = y1; zs = z1; ts = t1;
      }

      try {
        ctx.avalanche->AvalancheElectron(xs, ys, zs, ts, 0.1, 0, 0, 0);
      } catch (const Garfield::Exception& e) {
        // Thrown by AvalancheMicroscopic::TransportElectrons when the
        // starting point isn't inside a valid medium (e.g. a GEANT4
        // cluster landing outside the gas region, or right on a mesh
        // boundary). Skip this electron rather than aborting the run.
        ++nSkippedInvalidMedium;
        continue;
      }

      ++results.nTried;
      bool entered = false;
      const size_t nEnd = ctx.avalanche->GetNumberOfElectronEndpoints();
      for (size_t j = 0; j < nEnd && !entered; ++j) {
        double x0, y0, z0, t0, e0, x1, y1, z1, t1, e1;
        int st = 0;
        ctx.avalanche->GetElectronEndpoint(j, x0, y0, z0, t0, e0,
                                           x1, y1, z1, t1, e1, st);
        if (z1 < kWellEntryZ) entered = true;
      }

      for (const auto& ion : ctx.avalanche->GetIons()) {
        const auto& p0 = ion.path.front();
        bool ok = ctx.driftIon->DriftIon(p0.x, p0.y, p0.z, p0.t);
        if (ok) ++results.nIonsDrifted; else ++results.nIonsFailed;
      }

      int ne = 0, ni = 0;
      ctx.avalanche->GetAvalancheSize(ne, ni);
      if (entered) {
        ++results.nEntered;
        results.gainsRaw.push_back(ne);
      }
      results.neTotal += ne;
      results.niTotal += ni;
      results.gains.push_back(ne);
    }
  };

  const size_t nSeeds = seeds.size();
  // Cap at ~200 printed updates regardless of nSeeds, so this doesn't
  // slow down large TrackHeed runs by flooding stdout every electron.
  const size_t progressEvery = std::max<size_t>(1, nSeeds / 200);
  const std::time_t tRunStart = std::time(nullptr);
  std::time_t tLastPrint = tRunStart;

  for (size_t i = 0; i < nSeeds; ++i) {
    const auto& e = seeds[i];
    processElectron(e.x, e.y, e.z, e.t, e.nElectrons);

    const bool countTrigger = (i + 1) % progressEvery == 0 || i + 1 == nSeeds;
    const bool firstElectron = (i == 0);  // instant feedback that the loop is alive
    const bool timeTrigger =
        std::difftime(std::time(nullptr), tLastPrint) >= kProgressIntervalSeconds;

    if (countTrigger || firstElectron || timeTrigger) {
      PrintProgressBar(i + 1, nSeeds, tRunStart);
      tLastPrint = std::time(nullptr);
    }
  }

  std::cout << "Total electrons: " << results.neTotal
            << "  ions: " << results.niTotal
            << "  drifted OK: " << results.nIonsDrifted
            << "  failed: " << results.nIonsFailed << "\n";
  if (kUseHybridDrift) {
    std::cout << "Hybrid drift: " << results.nOutsideMap
              << " outside-map, " << results.nLostInDrift
              << " lost in the drift gap (attachment/boundary).\n";
  }
  if (nSkippedInvalidMedium > 0) {
    std::cout << "Skipped " << nSkippedInvalidMedium
              << " electron(s): starting point not in a valid medium.\n";
  }

  int sizeLimit = ctx.avalanche->GetAvalancheSizeLimit();
  std::cout << "The current avalanche size limit is: " << sizeLimit << std::endl;

  return results;
}

}  // namespace urwell
