#include "DetectorSetup.h"
#include "Config.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace Garfield;

namespace urwell {

void SetupDetector(SimContext& ctx) {
  // ---- 1. Real field (drift/avalanche) ----
  ComponentComsol* fmReal = new ComponentComsol();
  fmReal->Initialise("mesh_real.mphtxt", "dielectrics.dat",
                      "field_real.txt", "mum");

  // ---- 2. Gas ----
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 80., "co2", 20.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->Initialise(true);
  gas->EnableDrift();
  gas->LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas->EnablePenningTransfer();

  // ---- 3. Associate materials ----
  unsigned int nMat = fmReal->GetNumberOfMaterials();
  for (unsigned int i = 0; i < nMat; ++i) {
    std::cout << "Material " << i << ": eps_r = "
              << fmReal->GetPermittivity(i) << "\n";
  }
  fmReal->SetMedium(0, gas);
  fmReal->PrintRange();

  // ---- 4. Weighting field: SAME mesh as fmReal, different potential solve ----
  if (kUseDynamicWeightingField) {
    std::cout << "Using DYNAMIC (time-dependent) weighting field.\n";
    fmReal->SetDynamicWeightingPotential("field_weight_dynamic.txt", "anode");
  } else {
    std::cout << "Using STATIC weighting field (prompt signal only, no delayed component).\n";
    fmReal->SetWeightingField("field_weight.txt", "anode");
  }

  fmReal->EnablePeriodicityX();
  fmReal->EnablePeriodicityY();

  // ---- 5. Sensor ----
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fmReal);
  sensor->AddElectrode(fmReal, "anode");
  sensor->SetArea(-150e-3, -150e-3, -0.02, 150e-3, 150e-3, 2.0e-1);  // cm
  sensor->SetTimeWindow(ctx.tStart, (ctx.tStop - ctx.tStart) / ctx.nSteps,
                         ctx.nSteps);

  if (kUseDynamicWeightingField) {
    sensor->EnableDelayedSignal();  // only meaningful with a time-dependent weighting field

    // Time points at which the DELAYED signal is calculated. Deliberately
    // coarser than the nSteps prompt-signal binning -- see the
    // kNumDelayedSteps comment in Config.h for why this is safe here.
    std::vector<double> delayedTimes;
    const double dtDelayed = (ctx.tStop - ctx.tStart) / kNumDelayedSteps;
    for (int i = 0; i < kNumDelayedSteps; ++i) {
      delayedTimes.push_back(ctx.tStart + dtDelayed / 2 + i * dtDelayed);
    }
    sensor->SetDelayedSignalTimes(delayedTimes);
  }

  // ---- Sanity check: probe weighting field/potential BEFORE running ----
  {
    double wx, wy, wz, wpot;
    fmReal->WeightingField(0., 0., 0.01, wx, wy, wz, "anode");
    wpot = fmReal->WeightingPotential(0., 0., 0.01, "anode");
    std::cout << "Sanity check -- Wpot at (0,0,0.01 cm) = " << wpot
              << "  Wfield = (" << wx << ", " << wy << ", " << wz << ")\n";
    if (wpot == 0.) {
      std::cout << "WARNING: weighting potential is exactly zero here. "
                << "Check that field_weight.txt matches mesh_real.mphtxt "
                << "and that the electrode label 'anode' matches the COMSOL "
                << "boundary name used when exporting field_weight.txt.\n";
    }
  }

  // Direct check: is the delayed weighting potential actually non-zero
  // and evolving with time, independent of any plotting? Only meaningful
  // when a dynamic weighting field is actually loaded.
  if (kUseDynamicWeightingField) {
    for (double t : {0., 10., 50., 100., 150.}) {
      double dwp = fmReal->DelayedWeightingPotential(0., 0., 0.01, t, "anode");
      std::cout << "Delayed Wpot at (0,0,0.01), t=" << t << " ns: " << dwp << "\n";
    }
  }

  double ex, ey, ez, v;
  int status;
  Medium* m = nullptr;
  fmReal->ElectricField(0., 0., 0.002, ex, ey, ez, v, m, status);
  std::cout << "E field in well: (" << ex << ", " << ey << ", " << ez
            << ") V/cm,  |E| ~ " << std::sqrt(ex * ex + ey * ey + ez * ez)
            << " V/cm\n";

  ctx.fmReal = fmReal;
  ctx.gas = gas;
  ctx.sensor = sensor;
}

}  // namespace urwell
