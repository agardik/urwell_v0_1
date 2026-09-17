# u-RWELL Garfield++ Simulation

Simulation of avalanche multiplication, ion back-drift, and induced signal formation
in a u-RWELL detector, using a COMSOL-imported electric field map and the
Garfield++ toolkit.

## What it does

1. **Loads a real electric field** from a COMSOL export (`ComponentComsol`)
   defined by a mesh file, a dielectric properties file, and a field solution.
2. **Defines the gas** as Ar/CO₂ (80/20) at 293.15 K and 760 Torr using
   `MediumMagboltz`, with an externally supplied ion mobility table.
3. **Assigns the gas medium** to the appropriate mesh region and prints the
   permittivity of each material in the mesh for a sanity check.
4. **Loads a weighting field** (same mesh, different potential solve) for the
   "anode" electrode, used later to compute induced signals.
5. **Sets up a `Sensor`** with periodic boundary conditions in x and y, a
   defined active volume, and a 200 ns time window split into 2000 steps.
6. **Produces diagnostic plots**, each in its own `TCanvas`:
   - Detector mesh (y–z cross section at x = 0) with overlaid drift lines.
   - Total induced signal on the anode vs. time.
   - Electron-only and ion-only components of the induced signal.
   - Electric potential/field map (z–x cross section at y = 0).
   - Weighting potential map for the anode (same cross section).

## Requirements

- [ROOT](https://root.cern/) (tested with a version providing `TApplication`,
  `TCanvas`, `TGraph`, `TAxis`)
- [Garfield++](https://garfieldpp.web.cern.ch/garfieldpp/) built with:
  - `ComponentComsol`
  - `MediumMagboltz` (requires Magboltz)
  - `TrackHeed` (requires Heed++)
  - `AvalancheMicroscopic`, `AvalancheMC`
  - `ViewField`, `ViewDrift`, `ViewSignal`, `ViewFEMesh`
- A C++17-capable compiler
- CMake (recommended, if building alongside a Garfield++ CMake project)

## Input files

Place the following files in the working directory before running (not
included in this repository — see notes below):

| File | Purpose |
|---|---|
| `mesh_real.mphtxt` | COMSOL mesh export |
| `dielectrics.dat` | Material/dielectric property definitions matched to mesh regions |
| `field_real.txt` | COMSOL electrostatic field solution (operating voltages) |
| `field_weight.txt` | COMSOL field solution with the anode set to unit potential, all other electrodes grounded (weighting field) |
| `IonMobility_Ar+_Ar.txt` | Ar⁺ ion mobility table for the drift gas |

**Important:** `field_real.txt` and `field_weight.txt` must be solved on the
*same* mesh (`mesh_real.mphtxt`), and the electrode name used when exporting
the weighting field from COMSOL must exactly match the string `"anode"` used
in the code (`SetWeightingField`, `AddElectrode`, `WeightingField`, etc.). The
code includes a built-in sanity check at startup that prints a warning if the
weighting potential evaluates to zero at a sample point — if you see that
warning, check the mesh/electrode-name consistency before trusting any of the
signal output.

## Building

Example build against an existing Garfield++ installation (adjust include/
library paths as needed for your setup, e.g. via `garfield-config` or the
Garfield++ CMake package):

```bash
g++ -std=c++17 urwell_sim.cc -o urwell_sim \
    $(garfield-config --cflags --libs) \
    $(root-config --cflags --libs)
```

Or, if this file lives inside a Garfield++-style CMake project, add it as an
executable target in `CMakeLists.txt` and build with the usual:

```bash
mkdir build && cd build
cmake ..
make
```

## Running

```bash
./urwell_sim
```

The program opens a ROOT `TApplication` and displays five interactive
canvases:

- `cMesh` — detector cross section with avalanche/ion drift lines
- `cSignal` — total induced signal on the anode
- `cSignalElectron` — electron-induced signal component
- `cSignalIon` — ion-induced signal component
- `cField` — electric potential/field map
- `cWeight` — weighting potential map for the anode

Console output includes per-material permittivity, the weighting-field sanity
check, and totals for electrons produced, ions produced, and ions
successfully vs. unsuccessfully drifted.


## License

Add a license of your choice (e.g. MIT, GPLv3) here.
