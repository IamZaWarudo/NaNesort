#include <GCommands.h>
#include <GCanvas.h>
#include <GH1D.h>
#include <TSpectrum.h>
#include <TString.h>
#include <iostream>

void Background1D(int iterations = 20)
{
    auto* hist = dynamic_cast<TH1D*>(GrabHist());

    if (!hist) {
        std::cerr << "Select a TH1D or GH1D spectrum first.\n";
        return;
    }

    if (iterations < 1 ||
        iterations > (hist->GetNbinsX() - 1) / 2) {
        std::cerr << "Invalid number of iterations.\n";
        return;
    }

    static int run = 0;
    TString name = Form("background_%d", ++run);

    auto* original = new GH1D(*hist);
    original->SetName(name + "_original");
    original->SetDirectory(nullptr);
    original->SetStats(false);
    original->GetXaxis()->SetRange(0, 0);

    TSpectrum estimator;
    auto* background = estimator.Background(
        original, iterations,
        "BackIncreasingWindow BackOrder2 nosmoothing");

    if (!background) {
        std::cerr << "Background estimation failed.\n";
        delete original;
        return;
    }

    background->SetName(name + "_estimate");
    background->SetDirectory(nullptr);
    background->SetStats(false);

    auto* net = new GH1D(*original);
    net->SetName(name + "_net");
    net->SetDirectory(nullptr);
    net->SetStats(false);
    net->SetTitle("Background-subtracted spectrum");

    for (int i = 1; i <= original->GetNbinsX(); ++i) {
        net->SetBinContent(
            i, original->GetBinContent(i)
             - background->GetBinContent(i));
    }

    auto* cOriginal = new GCanvas(
        (name + "_original_canvas").Data(),
        "Original and background", 1000, 600);

    cOriginal->cd();
    original->SetTitle("continuem check");
    original->SetLineColor(kBlack);
    original->Draw("HIST");

    background->SetLineColor(kRed);
    background->SetLineWidth(2);
    background->Draw("HIST SAME");

    cOriginal->Modified();
    cOriginal->Update();

    auto* cNet = new GCanvas(
        (name + "_net_canvas").Data(),
        "Background-subtracted spectrum", 1000, 600);

    cNet->cd();
    net->SetLineColor(kBlue);
    net->Draw("HIST");

    cNet->Modified();
    cNet->Update();
}
