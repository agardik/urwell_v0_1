// Polya fit of the per-electron gain distribution.
// Usage:  root -l 'polya_fit.C("output_data/urwell_output.root")'
//         root -l 'polya_fit.C("output_data/urwell_output.root", "gainRawTree", "", 60)'
//         root -l 'polya_fit.C("output_data/urwell_output.root", "gainRawTree", "", 200, 0., 10)'
//
// Polya (gamma) distribution, with mean gain Gbar and shape theta:
//   P(G) = 1/Gbar * (theta+1)^(theta+1)/Gamma(theta+1) * (G/Gbar)^theta * exp(-(theta+1) G/Gbar)
// theta = 0 gives the exponential (Furry) distribution; larger theta = narrower peak.
// Fit function = N * binWidth * P(G), so the parameters are:
//   [0] = N (number of entries), [1] = Gbar (mean gain), [2] = theta
//
// This version uses ROOT's built-in gamma distribution, TMath::GammaDist:
//   TMath::GammaDist(x, gamma, mu, beta)
// with gamma = theta+1, mu = 0, beta = Gbar/(theta+1), which is exactly P(G).
//
// Bins with fewer than minCounts entries are excluded from the fit (they are given zero
// error, which ROOT's fitter skips). Two canvases are saved: linear and log-y views.

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

// ROOT built-in gamma distribution, parameterised as the Polya distribution.
double polya_root(double *x, double *p) {
  const double G = x[0], N = p[0], Gbar = p[1], th = p[2];
  if (G <= 0. || Gbar <= 0. || th <= -1.) return 0.;

  const double gamma = th + 1.;
  const double beta  = Gbar / gamma;

  return N * TMath::GammaDist(G, gamma, 0., beta);
}

void polya_fit(const char *fileName = "/Users/alexandrosgardikioti/garfieldpp/urwell_v0_1/output_data/urwell_output.root",
               const char *treeName = "gainRawTree",
               const char *branchName = "",   // empty -> first branch of the tree
               int nBins = 200,       // target number of bins; width is rounded up to an integer >= 1
               double gainMin = 0.,   // ignore gains below this value in the fit (0 = use everything)
               double minCounts = 10.) {  // ignore bins with fewer than this many entries in the fit
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
  std::cout << "Gain range [" << gmin << ", " << gmax << "], requested nBins = " << nBins
            << " -> bin width = " << w << ", actual bins = " << nb << "\n";
  TH1D *h = new TH1D("hGainFit", ";Gain (electrons per primary electron);Entries",
                     nb, lo, lo + nb * w);
  t->Draw((br + ">>hGainFit").c_str(), "", "goff");
  if (h->GetEntries() == 0) { std::cerr << "No entries\n"; return; }

  const double binW = h->GetBinWidth(1);
  const double mean = h->GetMean();

  const double hi = lo + nb * w;

  // Fit range starts at the lower edge of the bin containing gainMin, so no bin is cut in half.
  const bool cut = gainMin > lo;
  const int firstBin = cut ? h->FindBin(gainMin) : 1;
  const double fitLo = h->GetBinLowEdge(firstBin);
  if (cut) {
    const double nIn = h->Integral(firstBin, nb);
    std::cout << "Fit range starts at G = " << fitLo << " (gainMin = " << gainMin << "); "
              << nIn << " of " << h->GetEntries() << " entries used\n";
    if (nIn <= 0.) { std::cerr << "No entries above gainMin\n"; return; }
  }
  const double nUsed = cut ? h->Integral(firstBin, nb) : h->GetEntries();

  // Clone the histogram for the fit. Bins with fewer than minCounts entries get their error
  // set to zero, which tells ROOT's fitter to skip them. The display histogram h is untouched.
  TH1D *hFit = (TH1D *)h->Clone("hGainFitForFit");
  hFit->SetDirectory(0);
  int nExcluded = 0, nKept = 0;
  for (int i = 1; i <= hFit->GetNbinsX(); ++i) {
    if (hFit->GetBinContent(i) < minCounts) {
      hFit->SetBinError(i, 0.);
      ++nExcluded;
    } else {
      ++nKept;
    }
  }
  std::cout << "Fit uses " << nKept << " bins; " << nExcluded
            << " bins with < " << minCounts << " entries are excluded\n";
  if (nKept < 4) {
    std::cerr << "Too few bins with >= " << minCounts << " entries to constrain a 3-parameter fit\n";
    return;
  }
  const bool anyExcluded = nExcluded > 0;

  TF1 *fit = new TF1("fPolya", polya_root, fitLo, hi, 3);
  fit->SetParNames("N", "Gbar", "theta");
  fit->SetParameters(nUsed * binW, mean, 1.);
  fit->SetParLimits(2, -0.99, 50.);
  // N can only be fixed to the total entries when *all* bins in the fit range contribute.
  // If any bin is excluded (either by gainMin or by minCounts), the Polya is no longer
  // normalised over the fitted data, so N is left free.
  if (cut || anyExcluded) {
    fit->SetParLimits(0, 0., 10. * nUsed * binW);
  } else {
    fit->FixParameter(0, nUsed * binW);
  }
  fit->SetLineColor(kRed + 1);
  fit->SetLineWidth(3);

  // L = likelihood, I = integrate function over each bin, E = Minos errors.
  hFit->Fit(fit, "REIM0");

  // Pulls: (data - fit)/sqrt(fit), fit integrated over the bin. Computed for all bins,
  // including excluded ones, so you can see where the fit is failing outside its active set.
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

  TLine *lc = nullptr;
  if (cut) {
    // Dashed extrapolation of the fit into the ignored region, plus a line at the cut.
    TF1 *ext = new TF1("fPolyaExt", polya_root, lo, fitLo, 3);
    ext->SetParameters(fit->GetParameters());
    ext->SetLineColor(kRed + 1); ext->SetLineStyle(2); ext->SetLineWidth(2);
    ext->Draw("SAME");
    lc = new TLine(fitLo, 0, fitLo, h->GetMaximum());
    lc->SetLineColor(kGray + 2); lc->SetLineStyle(7); lc->Draw();
  }

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
  // With the "L" option, GetChisquare() is the Baker-Cousins likelihood chi2 (2*FCN).
  const int ndf = fit->GetNDF();
  const double chi2 = fit->GetChisquare();
  lt.DrawLatex(0.58, 0.45, Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2, ndf, ndf > 0 ? chi2 / ndf : 0.));
  // p-value: probability of a chi2 at least this large if the Polya model were correct.
  const double pval = ndf > 0 ? TMath::Prob(chi2, ndf) : 0.;
  lt.DrawLatex(0.58, 0.40, Form("p-value = %.2e", pval));
  if (cut) lt.DrawLatex(0.58, 0.35, Form("Fit range: G #geq %.0f", fitLo));
  if (anyExcluded) lt.DrawLatex(0.58, cut ? 0.30 : 0.35,
                                Form("Excluded: bins with < %.0f entries", minCounts));
  std::cout << "chi2/ndf = " << chi2 << " / " << ndf << ", p-value = " << pval << "\n";

  p2->cd();
  hp->SetTitle(";Gain;Pull");
  hp->SetFillColor(kGray + 1);
  hp->GetXaxis()->SetTitleSize(0.1); hp->GetXaxis()->SetLabelSize(0.09);
  hp->GetYaxis()->SetTitleSize(0.1); hp->GetYaxis()->SetLabelSize(0.09);
  hp->GetYaxis()->SetTitleOffset(0.45);
  hp->Draw("HIST");
  TLine *l0 = new TLine(hp->GetXaxis()->GetXmin(), 0, hp->GetXaxis()->GetXmax(), 0);
  l0->SetLineColor(kRed + 1); l0->Draw();
  if (cut) {
    TLine *lcp = new TLine(fitLo, hp->GetMinimum(), fitLo, hp->GetMaximum());
    lcp->SetLineColor(kGray + 2); lcp->SetLineStyle(7); lcp->Draw();
  }

  // ----- Linear view -----
  c->Update();
  c->SaveAs("gain_polya_fit.png");
  c->SaveAs("gain_polya_fit.pdf");

  // ----- Logarithmic view -----
  p1->cd();
  p1->SetLogy();
  // Log scale needs a positive lower limit; empty bins (content 0) simply won't be drawn.
  const double ymax = h->GetMaximum();
  h->SetMinimum(0.5);
  h->SetMaximum(ymax * 3.);
  // Extend the vertical cut line to the new top of the axis.
  if (lc) lc->SetY2(ymax * 3.);
  c->Update();
  c->SaveAs("gain_polya_fit_log.png");
  c->SaveAs("gain_polya_fit_log.pdf");

  std::cout << "Fit: Gbar = " << Gbar << " +- " << eG
            << ", theta = " << th << " +- " << eth << "\n";
}