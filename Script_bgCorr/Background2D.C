#include <GCommands.h>
#include <GCanvas.h>
#include <GH2D.h>
#include <TSpectrum2.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <vector>

void Background2D(int iterations = 10)
{
    auto* hist = dynamic_cast<TH2*>(GrabHist());

    if (!hist) {
        std::cerr << "Select a 2D histogram first.\n";
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

    // Copy counts into the array required by TSpectrum2.
    std::vector<std::vector<double>> counts(
        nx, std::vector<double>(ny));
    std::vector<double*> rows(nx);

    for (int x = 0; x < nx; ++x) {
        rows[x] = counts[x].data();

        for (int y = 0; y < ny; ++y) {
            double value = hist->GetBinContent(x + 1, y + 1);

            if (!std::isfinite(value) || value < 0) {
                std::cerr << "Input must contain finite, "
                             "nonnegative counts.\n";
                return;
            }

            counts[x][y] = value;
        }
    }

    // The input array is replaced by the estimated background.
    TSpectrum2 estimator;
    const char* error = estimator.Background(
        rows.data(), nx, ny,
        iterations, iterations,
        TSpectrum2::kBackIncreasingWindow,
        TSpectrum2::kBackSuccessiveFiltering);

    if (error) {
        std::cerr << error << '\n';
        return;
    }

    static int run = 0;
    TString name = Form("background2d_%d", ++run);

    auto* background = new GH2D(hist);
    background->SetName(name + "_background");
    background->SetTitle("Estimated 2D background");
    background->SetDirectory(nullptr);
    background->Reset();
    background->SetStats(false);

    auto* net = new GH2D(hist);
    net->SetName(name + "_net");
    net->SetTitle("Background-subtracted gamma-gamma matrix");
    net->SetDirectory(nullptr);
    net->Reset();
    net->SetStats(false);

    for (int x = 1; x <= nx; ++x) {
        for (int y = 1; y <= ny; ++y) {
            double b = counts[x - 1][y - 1];

            background->SetBinContent(x, y, b);
            net->SetBinContent(
                x, y, hist->GetBinContent(x, y) - b);

            // Input uncertainty only; background uncertainty is excluded.
            net->SetBinError(x, y, hist->GetBinError(x, y));
        }
    }

    background->GetXaxis()->SetRange(0, 0);
    background->GetYaxis()->SetRange(0, 0);
    net->GetXaxis()->SetRange(0, 0);
    net->GetYaxis()->SetRange(0, 0);

    auto* cBackground = new GCanvas(
        (name + "_background_canvas").Data(),
        "Estimated background", 1000, 800);

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
