#include "Plotting.h"
#include "Config.h"

#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"

#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1I.h"

#include <algorithm>
#include <iostream>

using namespace Garfield;

namespace urwell {

void PlotAmplificationDistribution(const AvalancheResults& results) {
  if (results.gains.empty()) return;

  const auto minmax = std::minmax_element(results.gains.begin(), results.gains.end());
  const int minGain = *minmax.first;
  const int maxGain = *minmax.second;
  const int nBins = std::max(10, std::min(60, maxGain - minGain + 1));

  TH1I* hGain = new TH1I("hAmplification", "Electron amplification distribution;Gain;Entries",
                        nBins, minGain - 0.5, maxGain + 0.5);
  for (const int gain : results.gains) {
    hGain->Fill(gain);
  }

  TCanvas* cGain = new TCanvas("cGain", "Electron amplification distribution", 800, 600);
  hGain->SetLineColor(kBlue + 1);
  hGain->SetFillColorAlpha(kAzure - 3, 0.4);
  hGain->Draw("HIST");
  cGain->Update();
}

void PlotMesh(SimContext& ctx) {
  if (!kShowMeshPlot) return;

  ViewFEMesh* viewMesh = new ViewFEMesh();
  viewMesh->SetComponent(ctx.fmReal);
  viewMesh->SetPlane(1, 0, 0, 0, 0, 0);  // normal (1,0,0) -> cut at x=0, shows y-z
  viewMesh->SetArea(-100e-4, -155.1e-4, 100e-4, 155.1e-4);  // y0,z0,y1,z1 (cm)
  viewMesh->SetFillMesh(true);
  viewMesh->EnableAxes();
  viewMesh->SetXaxisTitle("y [cm]");
  viewMesh->SetYaxisTitle("z [cm]");
  if (ctx.viewDrift) viewMesh->SetViewDrift(ctx.viewDrift);

  TCanvas* cMesh = new TCanvas("cMesh", "Mesh + avalanche (y-z, x=0)", 800, 800);
  viewMesh->SetCanvas(cMesh);

  viewMesh->SetColor(0, kCyan);      // gas
  viewMesh->SetColor(1, kYellow);    // kapton
  viewMesh->SetColor(2, kGray + 2);  // DLC
  viewMesh->SetColor(3, kGreen + 2); // prepreg
  viewMesh->SetColor(4, kOrange + 7);// conductor (copper/cathode/anode)

  viewMesh->Plot(true);
  cMesh->Update();
}

void PlotSignal(SimContext& ctx) {
  ViewSignal* viewSignal = new ViewSignal();
  viewSignal->SetSensor(ctx.sensor);
  TCanvas* cSignal = new TCanvas("cSignal", "Induced signal on anode", 800, 600);
  viewSignal->SetCanvas(cSignal);
  // Non-empty option strings for optPrompt/optDelayed are what makes
  // ViewSignal actually draw (and legend) those components alongside
  // the total. "t"/"p"/"d" are just placeholder ROOT draw-option
  // strings -- any non-empty string works. With a static field there is
  // no delayed component to draw, so optDelayed is left off in that case.
  if (kUseDynamicWeightingField) {
    viewSignal->PlotSignal("anode", "t", "p", "d");
  } else {
    viewSignal->PlotSignal("anode", "t", "p");
  }
  cSignal->Update();
}

void PlotElectronIonSignals(SimContext& ctx) {
  const int nSteps = ctx.nSteps;
  TGraph* grElectron = new TGraph(nSteps);
  TGraph* grIon = new TGraph(nSteps);
  for (int i = 0; i < nSteps; ++i) {
    const double t = ctx.tStart + (i + 0.5) * (ctx.tStop - ctx.tStart) / nSteps;
    grElectron->SetPoint(i, t, ctx.sensor->GetElectronSignal("anode", i));
    grIon->SetPoint(i, t, ctx.sensor->GetIonSignal("anode", i));
  }

  TCanvas* cSignalElectron = new TCanvas("cSignalElectron", "Electron-induced signal", 800, 600);
  grElectron->SetTitle("Electron-induced signal;Time [ns];Signal [fC/ns]");
  grElectron->SetLineColor(kOrange + 7);
  grElectron->SetLineWidth(2);
  grElectron->GetXaxis()->SetRangeUser(0, 2);
  grElectron->Draw("AL");
  cSignalElectron->Update();

  TCanvas* cSignalIon = new TCanvas("cSignalIon", "Ion-induced signal", 800, 600);
  grIon->SetTitle("Ion-induced signal;Time [ns];Signal [fC/ns]");
  grIon->SetLineColor(kBlue + 1);
  grIon->SetLineWidth(2);
  grIon->Draw("AL");
  cSignalIon->Update();
}

void PlotDelayedSignal(SimContext& ctx) {
  if (!kUseDynamicWeightingField) return;

  // Built directly from Sensor::GetDelayedElectronSignal/GetDelayedIonSignal,
  // the same way the electron/ion graphs above are built, so this is
  // guaranteed to show the exact same data ViewSignal's "delayed" curve
  // is drawing.
  // NOTE: these accessors are indexed against the PROMPT (nSteps)
  // binning, not kNumDelayedSteps -- Garfield internally maps the
  // coarser delayedTimes grid back onto that binning for you.
  const int nSteps = ctx.nSteps;
  TGraph* grDelayed = new TGraph(nSteps);
  for (int i = 0; i < nSteps; ++i) {
    const double t = ctx.tStart + (i + 0.5) * (ctx.tStop - ctx.tStart) / nSteps;
    const double de = ctx.sensor->GetDelayedElectronSignal("anode", i);
    const double di = ctx.sensor->GetDelayedIonSignal("anode", i);
    grDelayed->SetPoint(i, t, de + di);
  }

  TCanvas* cSignalDelayed = new TCanvas("cSignalDelayed", "Delayed (resistive) signal", 800, 600);
  grDelayed->SetTitle("Delayed signal;Time [ns];Signal [fC/ns]");
  grDelayed->SetLineColor(kGreen + 2);
  grDelayed->SetLineWidth(2);
  grDelayed->Draw("AL");
  cSignalDelayed->Update();
}

double PlotIntegratedCharge(SimContext& ctx) {
  // sensor->GetSignal(...) returns the induced *current* in each time
  // bin [fC/ns]; multiplying by the bin width and summing gives the
  // running total induced charge [fC]. This is a simple rectangle-rule
  // integral, which is fine here since the bins are uniform.
  // Note: with kUseDynamicWeightingField == true, GetSignal already
  // includes the delayed contribution, so this total is prompt+delayed
  // in that case, and prompt-only when using the static field.
  const int nSteps = ctx.nSteps;
  const double dt = (ctx.tStop - ctx.tStart) / nSteps;
  TGraph* grCharge = new TGraph(nSteps);
  double runningCharge = 0.;
  for (int i = 0; i < nSteps; ++i) {
    const double t = ctx.tStart + (i + 0.5) * dt;
    runningCharge += ctx.sensor->GetSignal("anode", i) * dt;
    grCharge->SetPoint(i, t, runningCharge);
  }
  const double totalCharge = runningCharge;  // charge at the end of the time window
  std::cout << "Total integrated charge on anode: " << totalCharge << " fC\n";

  TCanvas* cCharge = new TCanvas("cCharge", "Integrated charge on anode", 800, 600);
  grCharge->SetTitle("Integrated charge on anode;Time [ns];Charge [fC]");
  grCharge->SetLineColor(kBlack);
  grCharge->SetLineWidth(2);
  grCharge->Draw("AL");
  cCharge->Update();

  return totalCharge;
}

void PlotFields(SimContext& ctx) {
  ViewField* viewField = new ViewField();
  viewField->SetComponent(ctx.fmReal);
  viewField->SetPlane(0, -1, 0, 0, 0, 0);
  viewField->SetArea(-140e-4, -200.1e-5, 140e-4, 800.1e-4);

  TCanvas* cField = new TCanvas("cField", "Electric field |E| (z-x, y=0)", 800, 800);
  viewField->SetCanvas(cField);
  viewField->PlotContour("v");
  cField->Update();

  ViewField* viewWeight = new ViewField();
  viewWeight->SetComponent(ctx.fmReal);
  viewWeight->SetPlane(0, -1, 0, 0, 0, 0);
  viewWeight->SetArea(-140e-4, -200.1e-5, 140e-4, 800.1e-4);

  TCanvas* cWeight = new TCanvas("cWeight", "Weighting potential (z-x, y=0)", 800, 800);
  viewWeight->SetCanvas(cWeight);
  viewWeight->PlotContourWeightingField("anode", "v");
  cWeight->Update();
}

void PrintDelayedSignalSummary(SimContext& ctx) {
  if (!kUseDynamicWeightingField) return;

  const int nSteps = ctx.nSteps;
  double delayedTotal = 0.;
  int nNonZeroDelayed = 0;
  for (int i = 0; i < nSteps; ++i) {
    const double de = ctx.sensor->GetDelayedElectronSignal("anode", i);
    const double di = ctx.sensor->GetDelayedIonSignal("anode", i);
    delayedTotal += (de + di);
    if (de != 0. || di != 0.) ++nNonZeroDelayed;
  }
  std::cout << "Delayed signal: nonzero bins = " << nNonZeroDelayed
            << " / " << nSteps
            << ", summed delayed contribution = " << delayedTotal << "\n";
}

}  // namespace urwell
