#ifndef GEANT_HITS_H
#define GEANT_HITS_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

// One ionization point produced from a GEANT4 step.
// Positions are in cm (Garfield convention), time in ns.
struct IonizationPoint {
  double x;          // cm
  double y;          // cm
  double z;          // cm
  double t;          // ns
  int    n_electrons; // number of electron-ion pairs in this cluster
};

// Ionization clusters grouped by originating GEANT4 track.
struct TrackIonization {
  int track_id;
  int pdg_id;   // PDG particle ID (e.g. 1000020040 = alpha, 1000030060 = Li-6)
  double phi;   // azimuthal angle of initial step direction [rad]
  double theta; // polar angle of initial step direction from z-axis [rad]
  std::vector<IonizationPoint> points;
};

namespace geant_hits_detail {

// Max length (including null terminator) of a G4AnalysisManager string
// ntuple column, as stored in the ROOT backend's fixed-size char branch.
// 256 is comfortably larger than any Geant4 particle/process name.
constexpr int kStringBranchBufSize = 256;

// Try a list of tree-name candidates (to support different ntuple naming
// conventions across projects, e.g. "hits" vs "Hits").
inline TTree* FindHitsTree(TFile* f, const std::vector<std::string>& treeNameCandidates) {
  for (const auto& name : treeNameCandidates) {
    if (auto* t = dynamic_cast<TTree*>(f->Get(name.c_str()))) return t;
  }
  return nullptr;
}

// Try a list of scalar branch-name candidates; bind the first one found to addr.
// Returns the branch name actually used, or "" if none of the candidates exist.
// NOT for string branches -- see BindCharBranch below.
inline std::string BindBranch(TTree* tree, const std::vector<std::string>& candidates, void* addr) {
  for (const auto& name : candidates) {
    if (tree->FindBranch(name.c_str())) {
      tree->SetBranchAddress(name.c_str(), addr);
      return name;
    }
  }
  return "";
}

// G4AnalysisManager string ntuple columns are stored as a fixed-size Char_t
// branch (ROOT leaflist type "/C"), not as std::string objects, so they need
// SetBranchAddress bound to a char buffer.
inline std::string BindCharBranch(TTree* tree, const std::vector<std::string>& candidates, char* buf) {
  for (const auto& name : candidates) {
    if (tree->FindBranch(name.c_str())) {
      tree->SetBranchAddress(name.c_str(), buf);
      return name;
    }
  }
  return "";
}

// Map a Geant4 particle name to a PDG (or PDG-like nuclear) code.
// Covers elementary/light-composite particles by name, and falls back to
// parsing generic ion names of the form "<ElementSymbol><MassNumber>"
// (e.g. "Li6", "C12", "B10") into the standard nuclear code
// 1000000000 + Z*10000 + A*10 (matches e.g. 1000020040 for alpha).
// Returns 0 if the name is not recognized.
inline int ParticleNameToPdgId(const std::string& name) {
  static const std::map<std::string, int> kKnown = {
    {"gamma", 22},
    {"e-", 11}, {"e+", -11},
    {"mu-", 13}, {"mu+", -13},
    {"neutron", 2112}, {"anti_neutron", -2112},
    {"proton", 2212}, {"anti_proton", -2212},
    {"pi+", 211}, {"pi-", -211}, {"pi0", 111},
    {"deuteron", 1000010020},
    {"triton",   1000010030},
    {"He3",      1000020030},
    {"alpha",    1000020040},
  };
  auto it = kKnown.find(name);
  if (it != kKnown.end()) return it->second;

  // Generic ion convention: leading element symbol + trailing mass number.
  size_t i = 0;
  while (i < name.size() && std::isalpha(static_cast<unsigned char>(name[i]))) ++i;
  if (i == 0 || i == name.size()) return 0;  // no letters, or no digits after them

  const std::string symbol = name.substr(0, i);
  const std::string digits = name.substr(i);
  if (!std::all_of(digits.begin(), digits.end(),
                    [](unsigned char c) { return std::isdigit(c); })) {
    return 0;
  }
  const int A = std::stoi(digits);

  // Element symbol -> Z. Extend this table if your simulation produces
  // heavier recoils/fragments than what's listed here.
  static const std::map<std::string, int> kElementZ = {
    {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7},
    {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13},
    {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}, {"K", 19},
    {"Ca", 20}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26}, {"Ni", 28},
  };
  auto zit = kElementZ.find(symbol);
  if (zit == kElementZ.end()) return 0;

  return 1000000000 + zit->second * 10000 + A * 10;
}

}  // namespace geant_hits_detail

// Read a GEANT4 hits ROOT file and return ionization points grouped by track.
//
// Parameters:
//   filename    - path to the ROOT file containing the hits TTree
//   event_id    - select only entries with this event_id; pass -1 for all events
//   w_value_eV  - W value of the gas in eV (use gas.GetW() from Garfield)
//   cluster_size - max number of electron-ion pairs per returned point
//   x_offset_mm, y_offset_mm, z_offset_mm
//               - coordinates (in the tree's native mm units) of the point
//                 you want to become the new origin (0,0,0) in the returned
//                 ionization points. Subtracted from every pre/post position
//                 BEFORE the mm -> cm conversion. Default 0 (no shift), so
//                 existing calls are unaffected.
//
// Coordinate convention in the tree: positions in mm, time in ns.
// Positions are converted to cm (Garfield convention) after the offset
// is applied.
//
// Supports two ntuple naming conventions (tree name and branch names):
//   - the plain "hits" tree (event_id, track_id, particle_pdg_id, edep,
//     x_pos_pre/post, y_pos_pre/post, z_pos_pre/post, t_pre, t_post)
//   - the NeutronSource "Hits" tree (fEvent, fTrackID, Edep,
//     fXpre/fXpost, fYpre/fYpost, fZpre/fZpost, fTpre, fTpost)
//
// PDG id resolution order: a numeric PDG-like branch (pdg_id /
// particle_pdg_id / particle_pdg / PDG / fPdgID) if present, otherwise
// derived from a string particle-name branch (fParticleName), otherwise 0.
//
// For each step with N = int(edep * 1e6 / W) electron-ion pairs:
//   - M clusters are placed uniformly along the line segment from
//     (x_pos_pre, y_pos_pre, z_pos_pre) to (x_pos_post, y_pos_post, z_pos_post)
//     (after the offset above is subtracted)
//   - M = ceil(N / cluster_size)
//   - The time at each cluster is linearly interpolated between t_pre and t_post
//     assuming constant particle speed along the segment
//   - Each IonizationPoint carries up to cluster_size pairs
//
// Steps with edep <= 0 or that round to zero pairs are skipped.
inline std::vector<TrackIonization> LoadGeantHitsByTrack(const std::string& filename,
                                                          int event_id,
                                                          double w_value_eV,
                                                          int cluster_size = 1,
                                                          double x_offset_mm = 0.,
                                                          double y_offset_mm = 0.,
                                                          double z_offset_mm = 0.) {
  using namespace geant_hits_detail;
  std::vector<TrackIonization> tracks;

  if (cluster_size < 1) cluster_size = 1;

  TFile* f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "LoadGeantHits: cannot open file '" << filename << "'\n";
    return tracks;
  }

  TTree* tree = FindHitsTree(f, {"hits", "Hits"});
  if (!tree) {
    std::cerr << "LoadGeantHits: no tree named 'hits' or 'Hits' found in '" << filename << "'\n";
    f->Close();
    return tracks;
  }

  Int_t    ev_id;
  Int_t    track_id;
  Int_t    pdg_id = 0;
  Double_t edep;
  Double_t x_pre,  y_pre,  z_pre,  t_pre;
  Double_t x_post, y_post, z_post, t_post;

  // Mandatory branches: try both naming conventions.
  const bool ok =
      !BindBranch(tree, {"event_id", "fEvent"}, &ev_id).empty() &&
      !BindBranch(tree, {"track_id", "fTrackID"}, &track_id).empty() &&
      !BindBranch(tree, {"edep", "Edep"}, &edep).empty() &&
      !BindBranch(tree, {"x_pos_pre", "fXpre"}, &x_pre).empty() &&
      !BindBranch(tree, {"y_pos_pre", "fYpre"}, &y_pre).empty() &&
      !BindBranch(tree, {"z_pos_pre", "fZpre"}, &z_pre).empty() &&
      !BindBranch(tree, {"t_pre", "fTpre"}, &t_pre).empty() &&
      !BindBranch(tree, {"x_pos_post", "fXpost"}, &x_post).empty() &&
      !BindBranch(tree, {"y_pos_post", "fYpost"}, &y_post).empty() &&
      !BindBranch(tree, {"z_pos_post", "fZpost"}, &z_post).empty() &&
      !BindBranch(tree, {"t_post", "fTpost"}, &t_post).empty();

  if (!ok) {
    std::cerr << "LoadGeantHits: required branch missing in '" << filename
              << "' (checked both 'hits' and 'Hits' naming conventions)\n";
    f->Close();
    return tracks;
  }

  // Prefer a numeric PDG branch if one exists.
  const std::string pdg_branch_name =
      BindBranch(tree, {"pdg_id", "particle_pdg_id", "particle_pdg", "PDG", "fPdgID"}, &pdg_id);

  // Otherwise fall back to deriving pdg_id from a string particle-name branch.
  // G4AnalysisManager stores string columns as a fixed-size Char_t branch,
  // so bind to a plain char buffer (not std::string*).
  char particle_name_buf[kStringBranchBufSize];
  particle_name_buf[0] = '\0';
  std::string particle_name_branch;
  if (pdg_branch_name.empty()) {
    particle_name_branch = BindCharBranch(tree, {"fParticleName", "particle_name"}, particle_name_buf);
  }

  static bool pdg_branch_reported = false;
  static bool pdg_missing_warned = false;
  if (!pdg_branch_name.empty()) {
    if (!pdg_branch_reported) {
      std::cout << "LoadGeantHits: using numeric PDG branch '" << pdg_branch_name << "'\n";
      pdg_branch_reported = true;
    }
  } else if (!particle_name_branch.empty()) {
    if (!pdg_branch_reported) {
      std::cout << "LoadGeantHits: no numeric PDG branch found; deriving pdg_id from '"
                << particle_name_branch << "'\n";
      pdg_branch_reported = true;
    }
  } else if (!pdg_missing_warned) {
    std::cerr << "LoadGeantHits: no PDG branch or particle-name branch found; "
                 "pdg_id will default to 0\n";
    pdg_missing_warned = true;
  }

  if (x_offset_mm != 0. || y_offset_mm != 0. || z_offset_mm != 0.) {
    std::cout << "LoadGeantHits: applying origin offset (mm): ("
              << x_offset_mm << ", " << y_offset_mm << ", " << z_offset_mm << ")\n";
  }

  const Long64_t nEntries = tree->GetEntries();
  int total_pairs = 0;
  std::map<int, std::vector<IonizationPoint>> grouped;
  std::map<int, int>                          track_pdg;
  std::map<int, std::pair<double, double>>    track_angles; // phi, theta

  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    if (event_id != -1 && ev_id != event_id) continue;
    if (edep <= 0.) continue;

    const int n = static_cast<int>(std::lround(edep * 1.e6 / w_value_eV));
    if (n < 1) continue;

    // Apply the origin shift in the tree's native mm units.
    const double x_pre_s  = x_pre  - x_offset_mm;
    const double y_pre_s  = y_pre  - y_offset_mm;
    const double z_pre_s  = z_pre  - z_offset_mm;
    const double x_post_s = x_post - x_offset_mm;
    const double y_post_s = y_post - y_offset_mm;
    const double z_post_s = z_post - z_offset_mm;

    // Record PDG id and initial direction from the first step of each track.
    if (track_pdg.find(track_id) == track_pdg.end()) {
      int resolved_pdg = static_cast<int>(pdg_id);
      if (pdg_branch_name.empty() && !particle_name_branch.empty()) {
        resolved_pdg = ParticleNameToPdgId(std::string(particle_name_buf));
      }
      track_pdg[track_id] = resolved_pdg;

      const double dx = x_post_s - x_pre_s;
      const double dy = y_post_s - y_pre_s;
      const double dz = z_post_s - z_pre_s;
      const double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
      track_angles[track_id] = {
        std::atan2(dy, dx),
        (r > 0.) ? std::acos(dz / r) : 0.
      };
    }

    // Line segment in cm (mm -> cm: divide by 10), after the offset shift.
    const double ax = 0.1 * x_pre_s,  ay = 0.1 * y_pre_s,  az = 0.1 * z_pre_s;
    const double bx = 0.1 * x_post_s, by = 0.1 * y_post_s, bz = 0.1 * z_post_s;

    // Group pairs into clusters to reduce drift calls.
    const int n_clusters = (n + cluster_size - 1) / cluster_size;
    for (int k = 0; k < n_clusters; ++k) {
      const int first_pair = k * cluster_size;
      const int remaining = n - first_pair;
      const int n_in_cluster = remaining < cluster_size ? remaining : cluster_size;
      const double lam = (k + 0.5) / n_clusters;  // uniform in (0,1)

      IonizationPoint p;
      p.x = ax + lam * (bx - ax);
      p.y = ay + lam * (by - ay);
      p.z = az + lam * (bz - az);
      p.t = t_pre + lam * (t_post - t_pre);  // ns, keep the absolute simulation time here
      p.n_electrons = n_in_cluster;

      grouped[track_id].push_back(p);
    }

    total_pairs += n;
  }

  f->Close();

  size_t total_clusters = 0;
  tracks.reserve(grouped.size());
  for (auto& item : grouped) {
    total_clusters += item.second.size();

    // GEANT4 stores global simulation time, not detector-relative drift time.
    // Rebase each track so its earliest hit starts at t = 0 ns, which keeps the
    // drift histogram window near the physical detector response without
    // inflating the bin width or shifting the whole event to unrealistic times.
    double min_t = std::numeric_limits<double>::infinity();
    for (const auto& p : item.second) {
      if (p.t < min_t) min_t = p.t;
    }
    if (std::isfinite(min_t)) {
      for (auto& p : item.second) {
        p.t -= min_t;
      }
    }

    TrackIonization tr;
    tr.track_id = item.first;
    tr.pdg_id   = track_pdg.count(item.first)    ? track_pdg[item.first]            : 0;
    tr.phi      = track_angles.count(item.first) ? track_angles[item.first].first   : 0.;
    tr.theta    = track_angles.count(item.first) ? track_angles[item.first].second  : 0.;
    tr.points   = std::move(item.second);
    tracks.push_back(std::move(tr));
  }

  std::cout << "LoadGeantHits: loaded " << total_pairs << " electron-ion pairs"
            << " in " << total_clusters << " clusters"
            << " across " << tracks.size() << " tracks"
            << " from " << nEntries << " tree entries"
            << " (event_id=" << (event_id == -1 ? "all" : std::to_string(event_id)) << ")\n";

  return tracks;
}

// Return all unique event IDs present in the hits tree, in ascending order.
inline std::vector<int> GetGeantEventIds(const std::string& filename) {
  using namespace geant_hits_detail;
  std::vector<int> ids;

  TFile* f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "GetGeantEventIds: cannot open file '" << filename << "'\n";
    return ids;
  }

  TTree* tree = FindHitsTree(f, {"hits", "Hits"});
  if (!tree) {
    std::cerr << "GetGeantEventIds: no tree named 'hits' or 'Hits' found in '" << filename << "'\n";
    f->Close();
    return ids;
  }

  Int_t ev_id;
  if (BindBranch(tree, {"event_id", "fEvent"}, &ev_id).empty()) {
    std::cerr << "GetGeantEventIds: no event_id/fEvent branch found in '" << filename << "'\n";
    f->Close();
    return ids;
  }

  std::set<int> seen;
  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    seen.insert(static_cast<int>(ev_id));
  }

  f->Close();
  ids.assign(seen.begin(), seen.end()); // already sorted (std::set)
  std::cout << "GetGeantEventIds: found " << ids.size() << " events in '" << filename << "'\n";
  return ids;
}

// Backward-compatible flat view: concatenate all per-track clusters.
inline std::vector<IonizationPoint> LoadGeantHits(const std::string& filename,
                                                   int event_id,
                                                   double w_value_eV,
                                                   int cluster_size = 1,
                                                   double x_offset_mm = 0.,
                                                   double y_offset_mm = 0.,
                                                   double z_offset_mm = 0.) {
  std::vector<IonizationPoint> points;
  const auto tracks = LoadGeantHitsByTrack(filename, event_id, w_value_eV, cluster_size,
                                           x_offset_mm, y_offset_mm, z_offset_mm);
  for (const auto& tr : tracks) {
    points.insert(points.end(), tr.points.begin(), tr.points.end());
  }
  return points;
}

#endif // GEANT_HITS_H