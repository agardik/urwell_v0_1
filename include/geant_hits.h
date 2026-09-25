#ifndef GEANT_HITS_H
#define GEANT_HITS_H

#include <cmath>
#include <iostream>
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

// Read a GEANT4 hits ROOT file and return ionization points grouped by track.
//
// Parameters:
//   filename    - path to the ROOT file containing a TTree named "hits"
//   event_id    - select only entries with this event_id; pass -1 for all events
//   w_value_eV  - W value of the gas in eV (use gas.GetW() from Garfield)
//   cluster_size - max number of electron-ion pairs per returned point
//
// Coordinate convention in the tree: positions in mm, time in ns.
// Positions are converted to cm (Garfield convention).
//
// For each step with N = int(edep * 1e6 / W) electron-ion pairs:
//   - M clusters are placed uniformly along the line segment from
//     (x_pos_pre, y_pos_pre, z_pos_pre) to (x_pos_post, y_pos_post, z_pos_post)
//   - M = ceil(N / cluster_size)
//   - The time at each cluster is linearly interpolated between t_pre and t_post
//     assuming constant particle speed along the segment
//   - Each IonizationPoint carries up to cluster_size pairs
//
// Steps with edep <= 0 or that round to zero pairs are skipped.
inline std::vector<TrackIonization> LoadGeantHitsByTrack(const std::string& filename,
                                                          int event_id,
                                                          double w_value_eV,
                                                          int cluster_size = 1) {
  std::vector<TrackIonization> tracks;

  if (cluster_size < 1) cluster_size = 1;

  TFile* f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "LoadGeantHits: cannot open file '" << filename << "'\n";
    return tracks;
  }

  TTree* tree = dynamic_cast<TTree*>(f->Get("hits"));
  if (!tree) {
    std::cerr << "LoadGeantHits: TTree 'hits' not found in '" << filename << "'\n";
    f->Close();
    return tracks;
  }

  Int_t    ev_id;
  Int_t    track_id;
  Int_t    pdg_id = 0;
  Double_t edep;
  Double_t x_pre,  y_pre,  z_pre,  t_pre;
  Double_t x_post, y_post, z_post, t_post;

  tree->SetBranchAddress("event_id",   &ev_id);
  tree->SetBranchAddress("track_id",   &track_id);
  tree->SetBranchAddress("edep",       &edep);
  tree->SetBranchAddress("x_pos_pre",  &x_pre);
  tree->SetBranchAddress("y_pos_pre",  &y_pre);
  tree->SetBranchAddress("z_pos_pre",  &z_pre);
  tree->SetBranchAddress("t_pre",      &t_pre);
  tree->SetBranchAddress("x_pos_post", &x_post);
  tree->SetBranchAddress("y_pos_post", &y_post);
  tree->SetBranchAddress("z_pos_post", &z_post);
  tree->SetBranchAddress("t_post",     &t_post);
  const char* pdg_branch_name = nullptr;
  if (tree->FindBranch("pdg_id")) {
    pdg_branch_name = "pdg_id";
  } else if (tree->FindBranch("particle_pdg_id")) {
    pdg_branch_name = "particle_pdg_id";
  } else if (tree->FindBranch("particle_pdg")) {
    pdg_branch_name = "particle_pdg";
  } else if (tree->FindBranch("PDG")) {
    pdg_branch_name = "PDG";
  }

  static bool pdg_branch_reported = false;
  static bool pdg_missing_warned = false;
  if (pdg_branch_name) {
    tree->SetBranchAddress(pdg_branch_name, &pdg_id);
    if (!pdg_branch_reported) {
      std::cout << "LoadGeantHits: using PDG branch '" << pdg_branch_name << "'\n";
      pdg_branch_reported = true;
    }
  } else if (!pdg_missing_warned) {
    std::cerr << "LoadGeantHits: no PDG branch found (tried pdg_id, particle_pdg_id, particle_pdg, PDG); pdg_id will default to 0\n";
    pdg_missing_warned = true;
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

    // Record PDG id and initial direction from the first step of each track.
    if (track_pdg.find(track_id) == track_pdg.end()) {
      track_pdg[track_id] = static_cast<int>(pdg_id);
      const double dx = x_post - x_pre;
      const double dy = y_post - y_pre;
      const double dz = z_post - z_pre;
      const double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
      track_angles[track_id] = {
        std::atan2(dy, dx),
        (r > 0.) ? std::acos(dz / r) : 0.
      };
    }

    // Line segment in cm (mm -> cm: divide by 10)
    const double ax = 0.1 * x_pre,  ay = 0.1 * y_pre,  az = 0.1 * z_pre;
    const double bx = 0.1 * x_post, by = 0.1 * y_post, bz = 0.1 * z_post;

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
      p.t = t_pre + lam * (t_post - t_pre);  // ns, constant-speed interpolation
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
  std::vector<int> ids;

  TFile* f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "GetGeantEventIds: cannot open file '" << filename << "'\n";
    return ids;
  }

  TTree* tree = dynamic_cast<TTree*>(f->Get("hits"));
  if (!tree) {
    std::cerr << "GetGeantEventIds: TTree 'hits' not found in '" << filename << "'\n";
    f->Close();
    return ids;
  }

  Int_t ev_id;
  tree->SetBranchAddress("event_id", &ev_id);

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
                                                   int cluster_size = 1) {
  std::vector<IonizationPoint> points;
  const auto tracks = LoadGeantHitsByTrack(filename, event_id, w_value_eV, cluster_size);
  for (const auto& tr : tracks) {
    points.insert(points.end(), tr.points.begin(), tr.points.end());
  }
  return points;
}

#endif // GEANT_HITS_H