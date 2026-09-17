#pragma once

#include <cstddef>
#include <ctime>

namespace urwell {

// Simple in-place-looking progress bar with elapsed time / ETA. Call
// once per completed item; the caller is responsible for rate-limiting
// calls (see kProgressIntervalSeconds / progressEvery in
// AvalancheRunner.cc) so this doesn't flood stdout for million-electron
// runs.
void PrintProgressBar(size_t done, size_t total, std::time_t tStart);

}  // namespace urwell
