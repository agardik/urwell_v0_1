#pragma once

#include <vector>

#include "AvalancheRunner.h"

namespace urwell {

struct GainHistogram {
  std::vector<int> gains;
  std::vector<double> binEdges;
  std::vector<int> counts;
  double mean = 0.;
  double sigma = 0.;
  int nBins = 0;
  int minGain = 0;
  int maxGain = 0;
};

// Compute a per-electron amplification histogram from the AvalancheResults
// collected in RunAvalancheLoop(). The gain assigned to each seed is the
// final number of electrons reported by AvalancheMicroscopic::GetAvalancheSize.
GainHistogram ComputeAmplificationDistribution(const AvalancheResults& results,
                                              int nBins = 50);

// Print a compact summary of the gain histogram and a few representative bins.
void PrintAmplificationDistribution(const GainHistogram& histogram);

}  // namespace urwell
