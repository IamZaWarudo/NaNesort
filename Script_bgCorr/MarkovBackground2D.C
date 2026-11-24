#include <GCommands.h>
#include <GCanvas.h>
#include <GH2D.h>
#include <TH2.h>
#include <TSpectrum2.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <vector>

void MarkovBackground2D(int iterations = 10, int smoothingWindow = 3)
{
    auto* hist = dynamic_cast<TH2*>(GrabHist());

    if (!hist) {
        std::cerr << "Select the original 2D histogram first.\n";
        return;
    }

    const int nx = hist->GetNbinsX();
    const int ny = hist->GetNbinsY();

    if (iterations < 1 ||
        iterations > (nx - 1) / 2 ||
        iterations > (ny - 1) / 2) {
        std::cerr << "Invalid clipping-window extent.\n";
        return;
    }

    if (smoothingWindow < 1) {
        std::cerr << "The Markov smoothing window must be positive.\n";
        return;
    }

    std::vector<std::vector<double>> counts(nx, std::vector<double>(ny));
    std::vector<double*> rows(nx);

    for (int x = 0; x < nx; ++x) {
        rows[x] = counts[x].data();

        for (int y = 0; y < ny; ++y) {
            const double value = hist->GetBinContent(x + 1, y + 1);

            if (!std::isfinite(value) || value < 0.0) {
                std::cerr << "Input must contain finite, nonnegative counts.\n";
                return;
            }

            counts[x][y] = value;
        }
    }

    TSpectrum2 estimator;

    // SNIP replaces the copied counts with the background estimate.
    if (const char* error = estimator.Background(
            rows.data(), nx, ny, iterations, iterations,
            TSpectrum2::kBackIncreasingWindow,
            TSpectrum2::kBackSuccessiveFiltering)) {
        std::cerr << "Background estimation failed: " << error << '\n';
        return;
    }

    // Smooth the background estimate in place before subtraction.
    if (const char* error =
            estimator.SmoothMarkov(rows.data(), nx, ny, smoothingWindow)) {
        std::cerr << "Markov smoothing failed: " << error << '\n';
        return;
    }

    for (const auto& row : counts) {
        for (double value : row) {
            if (!std::isfinite(value) || value < 0.0) {
                std::cerr << "Smoothed background contains invalid bin values.\n";
                return;
            }
        }
    }

    static int run = 0;
    const TString name = Form("markov_background2d_%d", ++run);

    auto* background = new GH2D(hist);
    background->SetName((name + "_background").Data());
    background->SetTitle("SNIP background after Markov smoothing");
    background->SetDirectory(nullptr);
    background->Reset();
    background->SetStats(false);
    background->GetXaxis()->SetRange(0, 0);
    background->GetYaxis()->SetRange(0, 0);

    auto* net = new GH2D(hist);
    net->SetName((name + "_net").Data());
    net->SetTitle("Gamma-gamma matrix after subtracting smoothed SNIP background");
    net->SetDirectory(nullptr);
    net->Reset();
    net->SetStats(false);
    net->GetXaxis()->SetRange(0, 0);
    net->GetYaxis()->SetRange(0, 0);

    for (int x = 1; x <= nx; ++x) {
        for (int y = 1; y <= ny; ++y) {
            const double estimate = counts[x - 1][y - 1];

            background->SetBinContent(x, y, estimate);
            net->SetBinContent(x, y, hist->GetBinContent(x, y) - estimate);

            // Carry over the measured bin error only.
            net->SetBinError(x, y, hist->GetBinError(x, y));
        }
    }

    auto* cBackground = new GCanvas(
        (name + "_background_canvas").Data(),
        "SNIP background after Markov smoothing", 1000, 800);
    cBackground->cd();
    cBackground->SetLogz();
    background->Draw("COLZ");
    cBackground->Update();

    auto* cNet = new GCanvas(
        (name + "_net_canvas").Data(),
        "Background-subtracted matrix", 1000, 800);
    cNet->cd();
    net->Draw("COLZ");
    cNet->Update();
}
