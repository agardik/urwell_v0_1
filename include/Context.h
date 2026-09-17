#pragma once

#include "Garfield/ComponentComsol.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"

namespace urwell {

// Everything the various stages of the pipeline (detector setup,
// electron source, avalanche loop, plotting) need to share. Built up
// piece by piece as main.cc calls into each module -- nothing here owns
// its memory (matches the original file's raw `new` + TApplication
// lifetime), so there's no dtor to write.
struct SimContext {
  Garfield::ComponentComsol* fmReal = nullptr;
  Garfield::MediumMagboltz* gas = nullptr;
  Garfield::Sensor* sensor = nullptr;

  Garfield::AvalancheMicroscopic* avalanche = nullptr;
  Garfield::AvalancheMC* driftIon = nullptr;
  Garfield::ViewDrift* viewDrift = nullptr;  // only set if kStoreDriftLines

  double tStart = 0.;   // ns
  double tStop = 200.;  // ns
  int nSteps = 2000;    // prompt-signal (sensor time window) resolution
};

}  // namespace urwell
