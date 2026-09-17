#include "ProgressBar.h"

#include <iostream>
#include <iomanip>
#include <string>

namespace urwell {

void PrintProgressBar(size_t done, size_t total, std::time_t tStart) {
  const double frac = total > 0 ? double(done) / double(total) : 1.;
  const double elapsed = std::difftime(std::time(nullptr), tStart);
  const double eta = frac > 0. ? elapsed / frac - elapsed : 0.;

  const int barWidth = 40;
  const int filled = static_cast<int>(frac * barWidth);

  // std::fixed/std::setprecision below are STICKY manipulators -- they
  // permanently change std::cout's formatting for every print that comes
  // after this function returns, not just this call. Save the current
  // state and restore it before returning so nothing downstream (e.g.
  // "Total integrated charge on anode: " << totalCharge) gets silently
  // reformatted as a 0-decimal fixed-point number.
  const std::ios_base::fmtflags oldFlags = std::cout.flags();
  const std::streamsize oldPrecision = std::cout.precision();

  // NOTE: prints a normal newline-terminated line rather than an
  // in-place "\r" overwrite. Garfield/TrackHeed/Magboltz print a lot of
  // their own diagnostic messages during the avalanche loop, and those
  // were constantly clobbering an in-place progress bar before you ever
  // saw it. A plain scrolling line survives that interleaving -- you'll
  // see it appear in your terminal/log every progressEvery electrons.
  std::cout << ">>> PROGRESS [" << std::string(filled, '=')
            << std::string(barWidth - filled, ' ') << "] "
            << std::fixed << std::setprecision(1) << (frac * 100.) << "%  "
            << "electron " << done << "/" << total
            << std::setprecision(0)
            << "   elapsed=" << elapsed << "s"
            << "   ETA=" << eta << "s" << std::endl;

  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

}  // namespace urwell
