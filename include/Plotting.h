#pragma once

#include "Context.h"

namespace urwell {

// 2D y-z cross section: mesh + avalanche overlaid. No-op unless
// kShowMeshPlot is set.
void PlotMesh(SimContext& ctx);

// Induced signal on the anode: total (+ prompt + delayed if dynamic).
void PlotSignal(SimContext& ctx);

// Electron-only and ion-only signal, split from the same sensor, on
// their own canvases.
void PlotElectronIonSignals(SimContext& ctx);

// Delayed (resistive) signal on its own canvas. No-op unless
// kUseDynamicWeightingField is set.
void PlotDelayedSignal(SimContext& ctx);

// Integrated (cumulative) charge on the anode vs. time. Prints the
// final total and returns it [fC].
double PlotIntegratedCharge(SimContext& ctx);

// Electric field and weighting potential contour plots (z-x, y=0).
void PlotFields(SimContext& ctx);

// Prints the nonzero-bin / summed-contribution delayed-signal
// diagnostic. No-op unless kUseDynamicWeightingField is set.
void PrintDelayedSignalSummary(SimContext& ctx);

}  // namespace urwell
