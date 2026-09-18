#include "AmplificationDistribution.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace urwell {

GainHistogram ComputeAmplificationDistribution(const AvalancheResults& results,
                                              int nBins) {
  GainHistogram histogram;
  histogram.gains = results.gains;

  if (histogram.gains.empty()) return histogram;

  if (nBins <= 0) nBins = 1;
  histogram.nBins = nBins;
  histogram.minGain = *std::min_element(histogram.gains.begin(), histogram.gains.end());
  histogram.maxGain = *std::max_element(histogram.gains.begin(), histogram.gains.end());

  histogram.binEdges.resize(nBins + 1);
  histogram.counts.assign(nBins, 0);

  const int gainRange = std::max(1, histogram.maxGain - histogram.minGain);
  const double binWidth = static_cast<double>(gainRange) / nBins;

  for (int i = 0; i <= nBins; ++i) {
    histogram.binEdges[i] = histogram.minGain + static_cast<double>(i) * binWidth;
  }

  double sum = 0.;
  double sumSq = 0.;
  for (const int gain : histogram.gains) {
    sum += gain;
    sumSq += static_cast<double>(gain) * gain;

    int idx = 0;
    if (gainRange > 0) {
      idx = static_cast<int>((gain - histogram.minGain) / binWidth);
      if (idx >= nBins) idx = nBins - 1;
      if (idx < 0) idx = 0;
    }
    ++histogram.counts[idx];
  }

  const double mean = sum / histogram.gains.size();
  histogram.mean = mean;
  histogram.sigma = std::sqrt(std::max(0., sumSq / histogram.gains.size() - mean * mean));

  return histogram;
}

void PrintAmplificationDistribution(const GainHistogram& histogram) {
  if (histogram.gains.empty()) {
    std::cout << "No electron amplification samples were recorded." << std::endl;
    return;
  }

  std::cout << "Electron amplification distribution:\n"
            << "  entries: " << histogram.gains.size() << "\n"
            << "  mean gain: " << histogram.mean << "\n"
            << "  sigma: " << histogram.sigma << "\n"
            << "  min gain: " << histogram.minGain << "\n"
            << "  max gain: " << histogram.maxGain << std::endl;

  const int printBins = std::min<int>(10, histogram.nBins);
  for (int i = 0; i < printBins; ++i) {
    const double left = histogram.binEdges[i];
    const double right = histogram.binEdges[i + 1];
    std::cout << "  [" << left << ", " << right << ") : " << histogram.counts[i]
              << std::endl;
  }
}

}  // namespace urwell
