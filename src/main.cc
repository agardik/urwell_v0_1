#include "TApplication.h"

#include "Config.h"
#include "Context.h"
#include "DetectorSetup.h"
#include "ElectronSource.h"
#include "AvalancheRunner.h"
#include "AmplificationDistribution.h"
#include "Plotting.h"

using namespace urwell;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  SimContext ctx;
  ctx.tStart = 0.;    // ns
  ctx.tStop = 200.;   // ns
  ctx.nSteps = 2000;  // prompt-signal (sensor time window) resolution

  // ---- Detector, gas, sensor, weighting field, pre-run diagnostics ----
  SetupDetector(ctx);

  // ---- Avalanche + ion-drift machinery ----
  SetupAvalanche(ctx);

  // ---- Primary electron seeds (TrackHeed alpha track, or manual list) ----
  const std::vector<ElectronSeed> seeds = GenerateElectronSeeds(ctx);

  // ---- Run the avalanche/ion-drift loop with a progress bar ----
  const AvalancheResults results = RunAvalancheLoop(ctx, seeds);

  // ---- Gain distribution across the simulated electron population ----
  if (kComputeAmplificationDistribution) {
    const GainHistogram gainHistogram = ComputeAmplificationDistribution(results);
    PrintAmplificationDistribution(gainHistogram);
  }

  // ---- Plots and summary diagnostics ----
  if (kComputePlots) {
    if (kComputeGainHistogram) PlotAmplificationDistribution(results);
    if (kComputeDetectorPlots) {
      PlotMesh(ctx);
      PlotSignal(ctx);
      PlotElectronIonSignals(ctx);
      PlotDelayedSignal(ctx);
      PlotIntegratedCharge(ctx);
      PlotFields(ctx);
      PrintDelayedSignalSummary(ctx);
    }
  }

  app.Run(kTRUE);
  return 0;
}
