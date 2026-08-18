#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>

#include <TSystem.h>
#include <TOFCorrector.h>
#include <evtLoop.h>
#include <ddasLoop.h>
#include <ddasHit.h>
#include <GChannel.h>
#include <GHistogramer.h>
#include <GBCS.h>
#include <globals.h>
#include <Unpacker.h>
#include <utils.h>

void ProcessEvent(Unpacker& Event, GBCS& bcs, const std::vector<ddasHit>& event);

int main(int argc, char** argv) {
  if(argc < 2) { printf("usage: MakeTOF run-<RUN>-*.evt [...]\n"); return 1; }

  std::vector<std::string> inputFiles;
  for(int i = 1; i < argc; ++i) {
    std::filesystem::path argument(argv[i]);
    if(argument.extension().string() != ".evt") {
      printf("Unsupported input file: %s\n", argv[i]);
      return 1;
    }
    inputFiles.push_back(argument.string());
  }
  std::sort(inputFiles.begin(), inputFiles.end());

  std::string PATH = programPath();
  GChannel::ReadDetmap(Form("%s/../cals/detmapfix.tsv", PATH.c_str()));

  std::filesystem::path p(inputFiles.front());
  std::string stem = p.stem().string();

  int run = -1;
  if(std::sscanf(stem.c_str(), "run-%d-%*d", &run) != 1) {
    printf("Invalid EVT filename: %s\n", stem.c_str());
    return 1;
  }

  gSystem->mkdir("hist", true);
  std::string histFile = Form("hist/hist%04d.root", run);
  GHistogramer::Get().SetOutFile(histFile);

  evtLoop  reader(inputFiles, 500000, true);
  ddasLoop converter(reader, 200, 1);
  reader.Start();
  converter.Start();

  Unpacker Event;
  std::vector<ddasHit> event;
  GBCS bcs;

  while(!converter.Finished() || !converter.Empty()) {
    if(converter.TryPop(event)) {
      ProcessEvent(Event, bcs, event);
      event.clear();
    }
  }

  reader.Stop();
  converter.Stop();
  GHistogramer::Get().Close();

  printf("Run %d: wrote %s, generating TOF correction...\n", run, histFile.c_str());
  TOFCorrector tc;
  tc.MakeCorrectionFile(histFile);   // writes tof/tof%04d.tof

  return 0;
}

void ProcessEvent(Unpacker& Event, GBCS& bcs, const std::vector<ddasHit>& event) {
  Event.Reset();
  Event.Unpack(event);
  bcs.Reset();
  bcs.Fill(Event);

  bool hasPin = Event.fPin1.HasHit() && Event.fPin2.HasHit();
  if(!hasPin) return;

  double rawtof  = Event.fI2SPin1.fCharge;
  double runtime = Event.fPin1.fTimestamp / 1.e8;

  GHistogramer::Get().Fill("tof_time", 3600, 0, 7200,  runtime,
                                       4000, 0, 64000, rawtof);
}
