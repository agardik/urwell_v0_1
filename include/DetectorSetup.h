#pragma once

#include "Context.h"

namespace urwell {

// Builds the real-field component, gas, and sensor; loads the (static
// or dynamic) weighting field; and runs the pre-avalanche sanity-check
// diagnostics (weighting potential probe, delayed-Wpot probe, E-field
// probe). Populates ctx.fmReal, ctx.gas, ctx.sensor.
void SetupDetector(SimContext& ctx);

}  // namespace urwell
