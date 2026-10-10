// Polya fit of the per-electron gain distribution.
// Usage:  root -l 'polya_fit.C("output_data/urwell_output.root")'
//         root -l 'polya_fit.C("output_data/urwell_output.root", "gainRawTree", "", 60)'
//
// Polya (gamma) distribution, with mean gain Gbar and shape theta:
//   P(G) = 1/Gbar * (theta+1)^(theta+1)/Gamma(theta+1) * (G/Gbar)^theta * exp(-(theta+1) G/Gbar)
// theta = 0 gives the exponential (Furry) distribution; larger theta = narrower peak.
// Fit function = N * binWidth * P(G), so the parameters are:
//   [0] = N (number of entries), [1] = Gbar (mean gain), [2] = theta

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"
#include <algorithm>
#include <iostream>
#include <string>

double polya(double *x, double *p) {
  const double G = x[0], N = p[0], Gbar = p[1], th = p[2];
  if (G <= 0. || Gbar <= 0. || th <= -1.) return 0.;
  const double a = th + 1.;
  const double logP = a * TMath::Log(a) - TMath::LnGamma(a) + th * TMath::Log(G / Gbar)
                      - a * G / Gbar - TMath::Log(Gbar);
  return N * TMath::Exp(logP);
}

void polya_fit(const char *fileName = "/Users/alexandrosgardikioti/garfieldpp/urwell_v0_1/output_data/urwell_output.root",
               const char *treeName = "gainRawTree",
               const char *branchName = "",   // empty -> first branch of the tree
               int nBins = 1000) {   // target number of bins; width is rounded up to an integer >= 1
  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fileName << "\n"; return; }

  TTree *t = dynamic_cast<TTree *>(f->Get(treeName));
  if (!t) { std::cerr << "Tree " << treeName << " not found\n"; return; }

  std::string br = branchName;
  if (br.empty()) {
    TObjArray *bl = t->GetListOfBranches();
    if (!bl || bl->GetEntries() == 0) { std::cerr << "Tree has no branches\n"; return; }
    br = bl->At(0)->GetName();
  }
  std::cout << "Using " << treeName << "::" << br << " (" << t->GetEntries() << " entries)\n";

  const double gmin = t->GetMinimum(br.c_str());
  const double gmax = t->GetMaximum(br.c_str());
  // Gains are integers: use an integer bin width, with bin edges on half-integers,
  // otherwise bins contain alternating numbers of integer values (zigzag).
  const int w = std::max(1, (int)TMath::Ceil((gmax - gmin + 1.) / nBins));
  const int nb = (int)TMath::Ceil((gmax - gmin + 1.) / w);
  const double lo = gmin - 0.5;
  TH1D *h = new TH1D("hGainFit", ";Gain (electrons per primary electron);Entries",
                     nb, lo, lo + nb * w);
  t->Draw((br + ">>hGainFit").c_str(), "", "goff");
  if (h->GetEntries() == 0) { std::cerr << "No entries\n"; return; }

  const double binW = h->GetBinWidth(1);
  const double mean = h->GetMean();

  TF1 *fit = new TF1("fPolya", polya, lo, lo + nb * w, 3);
  fit->SetParNames("N", "Gbar", "theta");
  fit->SetParameters(h->GetEntries() * binW, mean, 1.);
  fit->SetParLimits(2, -0.99, 50.);
  fit->FixParameter(0, h->GetEntries() * binW);  // normalisation fixed by the entries
  fit->SetLineColor(kRed + 1);
  fit->SetLineWidth(3);

  // L = likelihood, I = integrate function over each bin, E = Minos errors.
  h->Fit(fit, "RLEI0");

  // Pulls: (data - fit)/sqrt(fit), fit integrated over the bin.
  TH1D *hp = (TH1D *)h->Clone("hPull");
  hp->Reset();
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double m = fit->Integral(h->GetBinLowEdge(i), h->GetBinLowEdge(i + 1)) / h->GetBinWidth(i);
    if (m > 0.) hp->SetBinContent(i, (h->GetBinContent(i) - m) / TMath::Sqrt(m));
  }

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("cPolya", "Gain distribution with Polya fit", 900, 800);
  TPad *p1 = new TPad("p1", "", 0, 0.3, 1, 1);
  p1->SetBottomMargin(0.02); p1->Draw();
  TPad *p2 = new TPad("p2", "", 0, 0, 1, 0.3);
  p2->SetTopMargin(0.03); p2->SetBottomMargin(0.3); p2->Draw();
  p1->cd();
  h->GetXaxis()->SetLabelSize(0);
  h->SetLineColor(kBlue + 1);
  h->SetFillColorAlpha(kAzure - 3, 0.35);
  h->Draw("HIST");
  fit->Draw("SAME");

  const double Gbar = fit->GetParameter(1), eG = fit->GetParError(1);
  const double th = fit->GetParameter(2), eth = fit->GetParError(2);

  TLegend *leg = new TLegend(0.55, 0.65, 0.89, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(h, "Simulation", "f");
  leg->AddEntry(fit, "Polya fit", "l");
  leg->Draw();

  TLatex lt; lt.SetNDC(); lt.SetTextSize(0.035);
  lt.DrawLatex(0.58, 0.60, Form("#bar{G} = %.1f #pm %.1f", Gbar, eG));
  lt.DrawLatex(0.58, 0.55, Form("#theta = %.2f #pm %.2f", th, eth));
  lt.DrawLatex(0.58, 0.50, Form("Sample mean = %.1f, RMS = %.1f", mean, h->GetRMS()));

  p2->cd();
  hp->SetTitle(";Gain;Pull");
  hp->SetFillColor(kGray + 1);
  hp->GetXaxis()->SetTitleSize(0.1); hp->GetXaxis()->SetLabelSize(0.09);
  hp->GetYaxis()->SetTitleSize(0.1); hp->GetYaxis()->SetLabelSize(0.09);
  hp->GetYaxis()->SetTitleOffset(0.45);
  hp->Draw("HIST");
  TLine *l0 = new TLine(hp->GetXaxis()->GetXmin(), 0, hp->GetXaxis()->GetXmax(), 0);
  l0->SetLineColor(kRed + 1); l0->Draw();
  std::cout << "Fit: Gbar = " << Gbar << " +- " << eG << ", theta = " << th << " +- " << eth << "\n";
  c->Update();
  c->SaveAs("gain_polya_fit.png");
  c->SaveAs("gain_polya_fit.pdf");
}