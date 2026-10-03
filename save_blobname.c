#include <TH1.h>
#include <TCutG.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>

void save_blobname(TH1* hist) {


  const char* names[10] = {
  
  "S1", //Cut0
  "A1", //Cut1
  "A2", //Cut2
  "A3", //Cut3
  "B1", //Cut4
  "B2", //Cut5
  "B3", //Cut6
  "S2", //Cut7
  "S3", //Cut8
  "S4" //Cut9
  
  };

  TFile file("pid.cuts", "RECREATE");
  file.cd();

  for (int i = 0; i < 10; i++) {
    auto cut = (TCutG*)hist->GetListOfFunctions()->FindObject(Form("cut%d", i));
    cut->SetName(names[i]);
    cut->Write();
  }

  file.Close();
}







