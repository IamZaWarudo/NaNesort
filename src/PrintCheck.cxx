// Add stuff`


#include <cstdio>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <filesystem>
#include <cmath>

#include <TSystem.h>
#include <memory>
#include <TOFCorrector.h>
#include <evtLoop.h>
#include <ddasLoop.h>
#include <ddasHit.h>
#include <GChannel.h>
#include <GHistogramer.h>
#include <GBCS.h>
#include <globals.h>
#include <PIDGates.h>

#include <TTree.h>
#include <TFile.h>

// for when I tree the TTree
static int isotopeCode(const std::string& iso) {
  if(iso == "Na") return 0;
  if(iso == "Ne") return 1;
  if(iso == "F")  return 2;
  return -1;    // for the things that are not in the TCut gate
}

struct EventRecord {
  int       type;        // 0 = implant, 1 = decay
  long long id;
  int       run;
  int       isotope;     // -1 if unidentified
  int       xStrip, yStrip;
  double    timestamp;
  double    energy;
  double    frontE, backE;
  int       frontMult, backMult;
  std::vector<double> gammas;
  std::vector<int>    gammaId;
  std::vector<double> abGammas;
  std::vector<int>    abId;
};

static TTree*      gTree         = nullptr;
static TFile*      gTreeFile     = nullptr;
static EventRecord gRec;
static int         gRun          = -1;
static long long   gEventCounter = 0;

static void FillTree(GBCS& bcs, int isotope);


void ProcessEvent(GBCS& bcs, const std::vector<ddasHit>& event, const TOFCorrector* tofCorrector);

int main(int argc, char** argv) {
  if(argc < 2) { printf("usage: PrintCheck file.evt [file.evt ...]\n"); return 1; }

   std::vector<std::string> inputFiles;
   std::string correctionFile;

   for(int i = 1; i < argc; ++i) {
     std::filesystem::path argument(argv[i]);
     const std::string extension = argument.extension().string();

     if(extension == ".evt") {
       inputFiles.push_back(argument.string());
     } else if(extension == ".tof") {
       if(!correctionFile.empty()) {
         printf("Only one .tof file may be supplied\n");
         return 1;
       }
       correctionFile = argument.string();
     } else {
       printf("Unsupported input file: %s\n", argv[i]);
       return 1;
     }
  }

  if(inputFiles.empty()) { printf("No EVT files supplied\n"); return 1; }
  std::sort(inputFiles.begin(), inputFiles.end());

  std::string homedir = std::getenv("HOME");
  GChannel::ReadDetmap(Form("%s/Sandbox/NaNesort/cals/detmapfix.tsv", homedir.c_str()));
  PIDGates::Get().Load(Form("%s/Sandbox/NaNesort/cals/pid_cuts.root", homedir.c_str()));

  std::filesystem::path p(inputFiles.front());
  std::string stem = p.stem().string();

  int run = -1;
  if(std::sscanf(stem.c_str(), "run-%d-%*d", &run) != 1) {
    printf("Invalid EVT filename: %s\n", stem.c_str());
    return 1;
  }

  gSystem->mkdir("hist", true);
  gSystem->mkdir("tofcorr", true);
  std::string ofile;
  if(correctionFile.empty())
    ofile = Form("hist/hist%04d.root", run);
  else
    ofile = Form("tofcorr/tofcorr%04d.root", run);

  GHistogramer::Get().SetOutFile(ofile);

  if(!correctionFile.empty()) {
    gRun = run;
    gSystem->mkdir("tree", true);
    gTreeFile = new TFile(Form("tree/tree%04d.root", run), "recreate");

    gTree = new TTree("events", "implant/decay events");
    gTree->Branch("type",      &gRec.type);
    gTree->Branch("id",        &gRec.id);
    gTree->Branch("run",       &gRec.run);
    gTree->Branch("isotope",   &gRec.isotope);
    gTree->Branch("xStrip",    &gRec.xStrip);
    gTree->Branch("yStrip",    &gRec.yStrip);
    gTree->Branch("timestamp", &gRec.timestamp);
    gTree->Branch("energy",    &gRec.energy);
    gTree->Branch("frontE",    &gRec.frontE);
    gTree->Branch("backE",     &gRec.backE);
    gTree->Branch("frontMult", &gRec.frontMult);
    gTree->Branch("backMult",  &gRec.backMult);
    gTree->Branch("gammas",    &gRec.gammas);
    gTree->Branch("gammaId",   &gRec.gammaId);
    gTree->Branch("abGammas", &gRec.abGammas);
    gTree->Branch("abId",     &gRec.abId);

    printf("Writing tree: tree/tree%04d.root\n", run);
  } else {
    printf("No .tof supplied - histograms only, no tree written.\n");
  }

  std::unique_ptr<TOFCorrector> tofCorrector;
    if(!correctionFile.empty()) {
      printf("Loading TOF correction: %s\n", correctionFile.c_str());
      tofCorrector = std::make_unique<TOFCorrector>(correctionFile);
    }

  evtLoop  reader(inputFiles, 500000, true);
  ddasLoop converter(reader, 200, 1);

  reader.Start();
  converter.Start();

  int64_t  lastPos    = 0;
  uint64_t lastBlocks = 0;
  uint64_t lastHits   = 0;
  auto lastTime  = std::chrono::steady_clock::now();
  auto lastPrint = lastTime;

  printf(HIDE_CURSOR);
  fflush(stdout);

  GBCS bcs;
  std::vector<ddasHit> event;

  while(!converter.Finished() || !converter.Empty()) {
    if(converter.TryPop(event)) {
      ProcessEvent(bcs, event, tofCorrector.get());
      event.clear();
    }

    auto now = std::chrono::steady_clock::now();
    if(now - lastPrint > std::chrono::milliseconds(500)) {
      auto e = reader.GetStats();
      auto d = converter.GetStats();

      constexpr int barWidth = 30;
      double percent = e.Percent();
      int filled = std::clamp(static_cast<int>((percent / 100.0) * barWidth), 0, barWidth);
      std::string bar(filled, '#');
      bar.append(barWidth - filled, '-');

      double dt = std::chrono::duration<double>(now - lastTime).count();
      double mbps = 0, blockRate = 0, hitRate = 0;
      if(dt > 0) {
        mbps      = (e.filePos    - lastPos)    / dt / 1024.0 / 1024.0;
        blockRate = (e.blocksRead - lastBlocks) / dt;
        hitRate   = (d.hitsBuilt  - lastHits)   / dt;
      }

      printf(CLEAR_LINE "[%s] %6.2f%%  file %llu/%llu  %.1f/%.1f MB  %7.1f MB/s\n",
          bar.c_str(), percent,
          (unsigned long long)e.currentFile, (unsigned long long)e.totalFiles,
          e.filePos / 1024.0 / 1024.0, e.fileSize / 1024.0 / 1024.0, mbps);
      printf(CLEAR_LINE "blocks=%llu (%6.0f/s)  hits=%llu (%6.0f/s)",
          (unsigned long long)e.blocksRead, blockRate,
          (unsigned long long)d.hitsBuilt,  hitRate);
      fflush(stdout);
      printf(CURSOR_UP);
      fflush(stdout);

      lastPos = e.filePos; lastBlocks = e.blocksRead; lastHits = d.hitsBuilt;
      lastTime = now; lastPrint = now;
    }
  }

  printf(CURSOR_DOWN "\n" SHOW_CURSOR);

  reader.Stop();
  converter.Stop();
  if(gTreeFile) {
    gTreeFile->cd();
    gTree->Write();
    gTreeFile->Close();}
  GHistogramer::Get().Close();
  return 0;
}

void ProcessEvent(GBCS& bcs, const std::vector<ddasHit>& event, const TOFCorrector* tofCorrector) {
  bcs.Reset();
  bcs.Unpack(event);
  int isotope = -1;                    

// all channel
  for(const auto& hit : event) {
    GHistogramer::Get().Fill("All_Channel/ecal", 16000, 0, 64000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());
    GHistogramer::Get().Fill("All_Channel/raw",  16000, 0, 64000, hit.GetCharge(),
                                                   300,   0, 300,   hit.GetId());}
  for(const auto& hit : bcs.fImplant) {
    GHistogramer::Get().Fill("All_Channel/Imp_ecal", 16000, 0, 64000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());}

  for(const auto& hit : bcs.fDecay) {
    GHistogramer::Get().Fill("All_Channel/Dec_ecal", 16000, 0, 64000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());}
  for(const auto& hit : bcs.fVeto) {
    GHistogramer::Get().Fill("All_Channel/Veto_ecal", 16000, 0, 64000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());}


// for TOF and PID only Pin1 and I2S used, because the other detectors failed and are unreliable
  if(bcs.fPin1.fTime > 10 && bcs.fI2S.fTime > 10) {
    const double runtime = bcs.fPin1.fTime / 1.e8;
    const double rawTOF  = bcs.fI2S.fCharge;
    const double dE      = bcs.fPin1.fEcal;

// TOF
    GHistogramer::Get().Fill("TOF/tof_raw", 3600, 0, 7200,  runtime,
                                            4000, 0, 64000, rawTOF);
// PID
    GHistogramer::Get().Fill("PID/pid_raw", 3600, 0, 64000, rawTOF,
                                            4000, 0, 16000, dE);
    double correctedTOF = rawTOF;
    if(tofCorrector) {
      correctedTOF = tofCorrector->Correct(rawTOF, runtime);

      GHistogramer::Get().Fill("TOF/tof_cr", 3600, 0, 7200,  runtime,
                                             4000, 0, 64000, correctedTOF);

      GHistogramer::Get().Fill("PID/pid_cr", 3600, 0, 64000, correctedTOF,
                                             4000, 0, 16000, dE);

      if(bcs.fPixel.ImplantValid()) {
      GHistogramer::Get().Fill("PID/pid_cr_implant", 3600, 0, 64000, correctedTOF,
                                                 4000, 0, 16000, dE);}
        }

   //Position with Pid Gate
      std::string iso = PIDGates::Get().Identify(correctedTOF, dE);
      isotope = isotopeCode(iso);
      if(!iso.empty() && bcs.fPixel.ImplantValid()) {
        GHistogramer::Get().Fill(Form("dssd/implant_%s", iso.c_str()),
                                       40, 0, 40, bcs.fPixel.ImplantPixel[0],
                                       40, 0, 40, bcs.fPixel.ImplantPixel[1]);

    GHistogramer::Get().Fill(Form("dssd/esum/implant_esum_%s", iso.c_str()),
                                       40, 0, 40, bcs.fPixel.ImplantPixel[0],
                                       40, 0, 40, bcs.fPixel.ImplantPixel[1],
                                       bcs.fPixel.ImplantEnergy);     
    }
  }


//Position
  if(bcs.fPixel.ImplantValid()) {
    GHistogramer::Get().Fill("dssd/implant_map", 40, 0, 40, bcs.fPixel.ImplantPixel[0],
                                                 40, 0, 40, bcs.fPixel.ImplantPixel[1]);

    GHistogramer::Get().Fill("dssd/esum/implant_esum", 40, 0, 40, bcs.fPixel.ImplantPixel[0],
                                                  40, 0, 40, bcs.fPixel.ImplantPixel[1],
                                                  bcs.fPixel.ImplantEnergy);}

  if(bcs.fPixel.DecayValid()) {
    GHistogramer::Get().Fill("dssd/decay_map", 40, 0, 40, bcs.fPixel.DecayPixel[0],
                                               40, 0, 40, bcs.fPixel.DecayPixel[1]);

    GHistogramer::Get().Fill("dssd/esum/decay_esum", 40, 0, 40, bcs.fPixel.DecayPixel[0],
                                                  40, 0, 40, bcs.fPixel.DecayPixel[1],
                                                  bcs.fPixel.DecayEnergy);}


// Gamma Spectra
  for(const auto& hit : bcs.fDecay) {
  switch(hit.GetId()) {
    case 192 ... 207:
      GHistogramer::Get().Fill("Gamma/LaBr_spec",16000,0,64000,hit.GetEcal());
      break;
    case 208 ... 271:
      GHistogramer::Get().Fill("Gamma/Clover_spec",16000,0,64000,hit.GetEcal()); 
      break;
    default:
      break;
    }

  for(const auto& h : bcs.fCrystal.hits)
    GHistogramer::Get().Fill("gamma/crystal", 8000, 0, 8000, h.fEcal,
                                                64, 0, 64,   h.fId);

  for(const auto& ab : bcs.fCrystal.addbackHits)
    GHistogramer::Get().Fill("gamma/addback", 8000, 0, 8000, ab.fEcal,
                                                16, 0, 16,   ab.fId);
}

//TTree
  FillTree(bcs, isotope);
  

}


//FIlline Tree

static void FillTree(GBCS& bcs, int isotope) {
  if(!gTree) return;                       // no .tof -> no tree

  const long long eventId = gEventCounter++;

  // ---- Implant ----
  if(bcs.fPixel.ImplantValid()) {
    gRec.type      = 0;
    gRec.id        = eventId;
    gRec.run       = gRun;
    gRec.isotope   = isotope;
    gRec.xStrip    = bcs.fPixel.ImplantFrontStrip;
    gRec.yStrip    = bcs.fPixel.ImplantBackStrip;
    gRec.timestamp = bcs.fPixel.ImplantTime;
    gRec.energy    = bcs.fPixel.ImplantEnergy;
    gRec.frontE    = bcs.fPixel.ImplantFrontE;
    gRec.backE     = bcs.fPixel.ImplantBackE;
    gRec.frontMult = bcs.fPixel.ImplantFrontMult;
    gRec.backMult  = bcs.fPixel.ImplantBackMult;
    gRec.gammas.clear();
    gRec.gammaId.clear();
    gRec.abGammas.clear();
    gRec.abId.clear();
    gTree->Fill();
  }

  // ---- Decay ----
  if(bcs.fPixel.DecayValid()) {
    gRec.type      = 1;
    gRec.id        = eventId;
    gRec.run       = gRun;
    gRec.isotope   = -1;                   // decays are not PID-identified
    gRec.xStrip    = bcs.fPixel.DecayFrontStrip;
    gRec.yStrip    = bcs.fPixel.DecayBackStrip;
    gRec.timestamp = bcs.fPixel.DecayTime;
    gRec.energy    = bcs.fPixel.DecayEnergy;
    gRec.frontE    = bcs.fPixel.DecayFrontE;
    gRec.backE     = bcs.fPixel.DecayBackE;
    gRec.frontMult = bcs.fPixel.DecayFrontMult;
    gRec.backMult  = bcs.fPixel.DecayBackMult;

    gRec.gammas.clear();
    gRec.gammaId.clear();
    for(const auto& h : bcs.fCrystal.hits) {
      gRec.gammas.push_back(h.fEcal);
      gRec.gammaId.push_back(h.fId);       // crystal 0-63
    }

    gRec.abGammas.clear();
    gRec.abId.clear();
    for(const auto& ab : bcs.fCrystal.addbackHits) {
      gRec.abGammas.push_back(ab.fEcal);
      gRec.abId.push_back(ab.fId);         // clover 0-15
    } 

    gTree->Fill();
  }
}
