#include <TH1.h>
#include <TCutG.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TKey.h>
#include <TVirtualPad.h>

void pidprint(){
  
  auto gates = TFile::Open("pid.cuts");

  if(!gates || gates->IsZombie() || !gPad)
    return;

  TIter next(gates->GetListOfKeys());

  while(auto key = (TKey*)next()) {
    auto cut = dynamic_cast<TCutG*>(key->ReadObj());
    if(!cut) continue;

    cut->SetLineColor(kRed);
    cut->SetLineWidth(2);
    cut->Draw("L SAME");
  
  
  int n = cut->GetN() - 1;
  if(n <= 0) continue;

  double x = 0;
  double y = 0;
  for(int i = 0; i < n; i++) {
    x += cut->GetX() [i];
    y += cut->GetY() [i];
  }
    x = x / n;
    y = y / n;
    
  auto label = new TLatex(x, y, cut->GetName());
  label->SetTextColor(kRed);
  label->SetTextSize(0.035);
  label->SetTextAlign(22);
  label->Draw();
  }

  gPad->Modified();
  gPad->Update();
}

