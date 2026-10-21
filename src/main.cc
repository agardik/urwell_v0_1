#include "TApplication.h"

#include <iostream>

#include "Config.h"
#include "Context.h"
#include "DetectorSetup.h"
#include "ElectronSource.h"
#include "AvalancheRunner.h"
#include "AmplificationDistribution.h"
#include "Plotting.h"

using namespace urwell;

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--ion-simulation" || arg == "--with-ions" || arg == "--ions") {
      kEnableIonSimulation = true;
    } else if (arg == "--no-ion-simulation" || arg == "--without-ions" ||
               arg == "--no-ions" || arg == "--electron-only") {
      kEnableIonSimulation = false;
    } else if (arg.rfind("--ion-simulation=", 0) == 0) {
      const std::string value = arg.substr(std::string("--ion-simulation=").size());
      kEnableIonSimulation = (value == "1" || value == "true" || value == "yes" || value == "on");
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: urwell [--ion-simulation[=true|false]|--no-ion-simulation]\n"
                << "  --ion-simulation      enable ion drift and ion signal simulation (default)\n"
                << "  --no-ion-simulation   skip ion drift, leaving only electron avalanches\n"
                << "  --ion-simulation=false disable ion drift while keeping the electron avalanche path\n";
      return 0;
    }
  }

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
    PrintGainAndTransparency(results);
  }

  // ---- Plots and summary diagnostics ----
  if (kComputePlots) {
    if (kComputeGainHistogram) PlotAmplificationDistribution(results);
    if (kComputeGainHistogram) PlotRawGainAndTransparency(results);
    if (kComputeDetectorPlots) {
      //PlotMesh(ctx);
      PlotSignal(ctx);
      PlotElectronIonSignals(ctx);
      PlotDelayedSignal(ctx);
      PlotIntegratedCharge(ctx);
      //PlotFields(ctx);
      PrintDelayedSignalSummary(ctx);
    }
  }

  // Write all relevant simulation quantities to a ROOT file for later
  // analysis / plotting outside the GUI session.
  WriteRootOutput(results, ctx, "urwell_output.root");

  app.Run(kTRUE);
  return 0;
}
