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


void ProcessEvent(GBCS& bcs, const std::vector<ddasHit>& event, const TOFCorrector* tofCorrector);

// Progress bar eta
static std::string FormatDuration(double seconds) {
  if(!std::isfinite(seconds) || seconds < 0) return "--:--";
  long long total = static_cast<long long>(seconds + 0.5);
  long long h = total / 3600, m = (total % 3600) / 60, s = total % 60;
  char buf[32];
  if(h > 0) std::snprintf(buf, sizeof(buf), "%lld:%02lld:%02lld", h, m, s);
  else      std::snprintf(buf, sizeof(buf), "%02lld:%02lld", m, s);
  return buf;
}


// Main function
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

  std::filesystem::path p(inputFiles.front());
  std::string stem = p.stem().string();

  int run = -1;
  if(std::sscanf(stem.c_str(), "run-%d-%*d", &run) != 1) {
    printf("Invalid EVT filename: %s\n", stem.c_str());
    return 1;
  }

  gSystem->mkdir("diagn", true);
  std::string ofile;
  if(correctionFile.empty())
    ofile = Form("diagn/diag%04d.root", run);
  else
    ofile = Form("diagn/diagCorr%04d.root", run);

  GHistogramer::Get().SetOutFile(ofile);

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
  const auto startTime = lastTime;
  double avgMBps = 0.0;

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

      if(mbps > 0) avgMBps = (avgMBps > 0) ? 0.7 * avgMBps + 0.3 * mbps : mbps;

      double elapsed  = std::chrono::duration<double>(now - startTime).count();
      double remainMB = (e.fileSize - e.filePos) / 1024.0 / 1024.0;
      double etaSec   = (avgMBps > 1e-6) ? remainMB / avgMBps : -1.0;
      if(etaSec < 0 && e.filePos > 0 && elapsed > 0)
        etaSec = remainMB / ((e.filePos / 1024.0 / 1024.0) / elapsed);
      std::string etaStr = (e.filePos >= e.fileSize)
                             ? std::string("draining")
                             : FormatDuration(etaSec);

      printf(CLEAR_LINE "[%s] %6.2f%%  file %llu/%llu  %.1f/%.1f MB  %7.1f MB/s\n",
          bar.c_str(), percent,
          (unsigned long long)e.currentFile, (unsigned long long)e.totalFiles,
          e.filePos / 1024.0 / 1024.0, e.fileSize / 1024.0 / 1024.0, mbps);
      printf(CLEAR_LINE "blocks=%llu (%6.0f/s)  hits=%llu (%6.0f/s)  elapsed %s  ETA %s",
          (unsigned long long)e.blocksRead, blockRate,
          (unsigned long long)d.hitsBuilt,  hitRate,
          FormatDuration(elapsed).c_str(),
          etaStr.c_str());      
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
  GHistogramer::Get().Close();
  return 0;
}

void ProcessEvent(GBCS& bcs, const std::vector<ddasHit>& event, const TOFCorrector* tofCorrector) {
  bcs.Reset();
  bcs.Unpack(event);

// all channel
  for(const auto& hit : event) {
    GHistogramer::Get().Fill("All_Channel/ecal", 10000, 0, 32000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());
    GHistogramer::Get().Fill("All_Channel/raw",  10000, 0, 32000, hit.GetCharge(),
                                                   300,   0, 300,   hit.GetId());
    if(bcs.fPin1.fTime > 0 && bcs.fPin2.fTime > 0){
    GHistogramer::Get().Fill("PIN/Pin1vs2", 4000,0,16000, bcs.fPin1.fEcal,
                                             4000,0,16000, bcs.fPin2.fEcal);}
    }
 
  for(const auto& hit : bcs.fImplant) {
    GHistogramer::Get().Fill("All_Channel/Implant_ecal", 10000, 0, 32000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());
    double implantfhbhdt =  std::abs(bcs.FrontHighGain.fTime - bcs.BackHighGain.fTime);
    GHistogramer::Get().Fill("dt/implantfhbhdt",  1000, 0,500, implantfhbhdt);
 
    GHistogramer::Get().Fill("PIN/Pin1vs2_veto", 10000,0,16000, bcs.fPin1.fEcal,
                                             10000,0,16000, bcs.fPin2.fEcal);
  }

  for(const auto& hit : bcs.fDecay) {
    GHistogramer::Get().Fill("All_Channel/Decay_ecal", 16000, 0, 64000, hit.GetEcal(),
                                                   300,   0, 300,   hit.GetId());
    double decayfhbhdt =  std::abs(bcs.FrontHighGain.fTime - bcs.BackHighGain.fTime);
    GHistogramer::Get().Fill("dt/decayfhbhdt",  1000, 0,500, decayfhbhdt);

  }

  for(const auto& hit : bcs.fPDCheck) {
  double dtfl = bcs.FrontLowGain.fTime - bcs.fPin1.fTime;
  double dtbl = bcs.BackLowGain.fTime - bcs.fPin1.fTime;
  double dtfh = bcs.FrontHighGain.fTime - bcs.fPin1.fTime;
  double dtbh = bcs.BackHighGain.fTime - bcs.fPin1.fTime;
  double dtfhbh = std::abs(bcs.FrontHighGain.fTime - bcs.BackHighGain.fTime);
  double dtflbl = std::abs(bcs.FrontLowGain.fTime - bcs.BackLowGain.fTime);
  double dtflfh = std::abs(bcs.FrontHighGain.fTime - bcs.FrontLowGain.fTime);
  double dtblbh = std::abs(bcs.BackHighGain.fTime - bcs.BackLowGain.fTime);
  double dtpin1n2 = bcs.fPin2.fTime - bcs.fPin1.fTime;

    GHistogramer::Get().Fill("dt/Pin-fl",  1000, 0,500, dtfl);
    GHistogramer::Get().Fill("dt/Pin-bl",  1000, 0,500, dtbl);
    GHistogramer::Get().Fill("dt/Pin-fh",  1000, 0,500, dtfh);
    GHistogramer::Get().Fill("dt/Pin-bh",  1000, 0,500, dtbh);
    GHistogramer::Get().Fill("dt/dtfhbh",  1000, 0,500, dtfhbh);
    GHistogramer::Get().Fill("dt/dtflbl",  1000, 0,500, dtflbl);
    GHistogramer::Get().Fill("dt/dtflfh",  1000, 0,500, dtflfh);
    GHistogramer::Get().Fill("dt/dtblbh",  1000, 0,500, dtblbh);
    GHistogramer::Get().Fill("dt/dtpin1n2", 1000,0,500, dtpin1n2);
}

  for(const auto& hit : bcs.fPSCheck) {
  double dtsssdlow = bcs.fSSSDLow.fTime - bcs.fPin1.fTime;
  double dtsssdhigh = bcs.fSSSDHigh.fTime - bcs.fPin1.fTime;
  double dtsssdlowhigh = std::abs(bcs.fSSSDHigh.fTime - bcs.fSSSDLow.fTime);


  double dt_dssd_sssd_low = bcs.fSSSDLow.fTime - bcs.FrontLowGain.fTime; 
    GHistogramer::Get().Fill("dt/Pin-sssdlow",  1000, 0,500, dtsssdlow);
    GHistogramer::Get().Fill("dt/Pin-sssdhigh",  1000, 0,500, dtsssdhigh);
    GHistogramer::Get().Fill("dt/dtsssdLowHigh",  1000, 0,500, dtsssdlowhigh);
  }

// PID & TOF
  if(bcs.fPin1.fTime > 10 && bcs.fPin2.fTime > 10 && bcs.fI2S.fTime > 10) {
    const double runtime = bcs.fPin1.fTime / 1.e8;
    const double rawTOF  = bcs.fI2S.fCharge;
    const double dE      = bcs.fPin1.fEcal;
    const double runtime2 = bcs.fPin2.fTime / 1.e8;
    const double dE2      = bcs.fPin2.fEcal;

      GHistogramer::Get().Fill("TOF/tof2_cr", 3600, 0, 7200,  runtime2,
                                             4000, 0, 64000, rawTOF);

      GHistogramer::Get().Fill("PID/pid2_cr", 3600, 0, 64000, rawTOF,
                                             4000, 0, 16000, dE2);

    double correctedTOF = rawTOF;
    if(tofCorrector) {
      correctedTOF = tofCorrector->Correct(rawTOF, runtime);

      GHistogramer::Get().Fill("TOF/tof1_cr", 3600, 0, 7200,  runtime,
                                             4000, 0, 64000, correctedTOF);

      GHistogramer::Get().Fill("PID/pid1_cr", 2000, 0, 64000, correctedTOF,
                                             3000, 0, 16000, dE);
        }
}

  for(const auto& hit : bcs.fImplant) {
    if(bcs.fPin1.fTime > 10 && bcs.fI2S.fTime > 10) {
    const double runtime = bcs.fPin1.fTime / 1.e8;
    const double rawTOF  = bcs.fI2S.fCharge;
    const double dE      = bcs.fPin1.fEcal;

    double correctedTOF = rawTOF;
    if(tofCorrector) {
      GHistogramer::Get().Fill("TOF/tof_with_veto", 3600, 0, 7200,  runtime,
                                             4000, 0, 64000, correctedTOF);

      correctedTOF = tofCorrector->Correct(rawTOF, runtime);
      GHistogramer::Get().Fill("PID/pid_with_veto", 2000, 0, 64000, correctedTOF,
                                             3000, 0, 16000, dE);}
        }
  }


//Position
  if(bcs.fPixel.Implant()) {
    GHistogramer::Get().Fill("dssd/implant_map", 40, 0, 40, bcs.fPixel.ImPixel[0],
                                                 40, 0, 40, bcs.fPixel.ImPixel[1]);
}

  if(bcs.fPixel.Decay()) {
    GHistogramer::Get().Fill("dssd/decay_map", 40, 0, 40, bcs.fPixel.DePixel[0],
                                               40, 0, 40, bcs.fPixel.DePixel[1]);
}


/*
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
       }
}
*/

}
