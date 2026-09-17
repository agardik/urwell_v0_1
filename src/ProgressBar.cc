#include "ProgressBar.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

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

  // Print an in-place-updating single-line progress bar. This uses '\r'
  // to return the cursor to the start of the line and flushes the output
  // so the bar updates in-place instead of spamming the terminal with
  // many newline-terminated lines. Track the previous printed length and
  // pad with spaces when the new line is shorter to ensure leftover
  // characters are cleared.
  std::ostringstream oss;
  oss << ">>> PROGRESS [" << std::string(filled, '=')
      << std::string(barWidth - filled, ' ') << "] "
      << std::fixed << std::setprecision(1) << (frac * 100.) << "%  "
      << "electron " << done << "/" << total
      << std::setprecision(0)
      << "   elapsed=" << elapsed << "s"
      << "   ETA=" << eta << "s";

  const std::string out = oss.str();
  static std::size_t prev_len = 0;
  const std::size_t pad = prev_len > out.size() ? prev_len - out.size() : 0;

  // Emit carriage return, the content, pad to clear previous longer output,
  // and flush. If we're complete, terminate the line with a newline.
  std::cout << '\r' << out << std::string(pad, ' ');
  if (done == total) {
    std::cout << std::endl;
    prev_len = 0;
  } else {
    std::cout << std::flush;
    prev_len = out.size();
  }

  std::cout.flags(oldFlags);
  std::cout.precision(oldPrecision);
}

}  // namespace urwell
