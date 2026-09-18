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

  ctx.avalanche = avalanche;
  ctx.driftIon = driftIon;
  ctx.viewDrift = viewDrift;
}

AvalancheResults RunAvalancheLoop(SimContext& ctx,
                                   const std::vector<ElectronSeed>& seeds) {
  AvalancheResults results;

  // ---- Helper: avalanche one electron + drift its resulting ions ----
  auto processElectron = [&](double x, double y, double z, double t) {
    ctx.avalanche->AvalancheElectron(x, y, z, t, 0.1, 0, 0, 0);

    for (const auto& ion : ctx.avalanche->GetIons()) {
      const auto& p0 = ion.path.front();
      bool ok = ctx.driftIon->DriftIon(p0.x, p0.y, p0.z, p0.t);
      if (ok) ++results.nIonsDrifted; else ++results.nIonsFailed;
    }

    int ne = 0, ni = 0;
    ctx.avalanche->GetAvalancheSize(ne, ni);
    results.neTotal += ne;
    results.niTotal += ni;
    results.gains.push_back(ne);
  };

  const size_t nSeeds = seeds.size();
  // Cap at ~200 printed updates regardless of nSeeds, so this doesn't
  // slow down large TrackHeed runs by flooding stdout every electron.
  const size_t progressEvery = std::max<size_t>(1, nSeeds / 200);
  const std::time_t tRunStart = std::time(nullptr);
  std::time_t tLastPrint = tRunStart;

  for (size_t i = 0; i < nSeeds; ++i) {
    const auto& e = seeds[i];
    processElectron(e.x, e.y, e.z, e.t);

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

  int sizeLimit = ctx.avalanche->GetAvalancheSizeLimit();
  std::cout << "The current avalanche size limit is: " << sizeLimit << std::endl;

  return results;
}

}  // namespace urwell
