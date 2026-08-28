#include <cstdio>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <filesystem>

#include <TSystem.h>
#include <TFile.h>
#include <TH2.h>
#include <TOFCorrector.h>
#include <evtLoop.h>
#include <ddasLoop.h>
#include <ddasHit.h>
#include <GChannel.h>
#include <GHistogramer.h>
#include <GBCS.h>
#include <Unpacker.h>
#include <FillHistogram.h>
#include <GCorrelator.h>
#include <globals.h>
#include <utils.h>


void ProcessEvent(Unpacker& unpacker, GBCS& bcs, GCorrelator& corr, std::vector<ddasHit>& event);


// ---- event type tally ----
enum { kNEventTypes = 6 };

static const char* kTypeName[kNEventTypes] = {
  "UNUSED", "IMPLANT", "DECAY", "LIGHTION", "VETO", "UNKNOWN"
};

static long gTypeCount[kNEventTypes] = {0};
static long gBadType = 0;   // EventType() out of range
static long gTotalEvents = 0;

static void PrintEventTypeSummary() {
  const long total = gTotalEvents ? gTotalEvents : 1;

  printf("\n--- Event type summary ---\n");
  for(int i = 0; i < kNEventTypes; ++i) {
    printf("  %-9s %d : %12ld  (%6.2f%%)\n",
           kTypeName[i], i, gTypeCount[i],
           100.0 * gTypeCount[i] / total);
  }
  if(gBadType) {
    printf("  %-9s   : %12ld  (%6.2f%%)\n",
           "BADTYPE", gBadType, 100.0 * gBadType / total);
  }
  printf("  %-11s : %12ld\n", "TOTAL", gTotalEvents);
}




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


int main(int argc, char** argv) {
  if(argc < 2) { printf("usage: sort file.evt [file.evt ...]\n"); return 1; }

  std::vector<std::string> inputFiles;

  for(int i = 1; i < argc; ++i) {
    std::filesystem::path argument(argv[i]);
    const std::string extension = argument.extension().string();

    if(extension == ".evt") {
      inputFiles.push_back(argument.string());
    } else {
      printf("Unsupported input file: %s\n", argv[i]);
      return 1;
    }
  }

  if(inputFiles.empty()) { printf("No EVT files supplied\n"); return 1; }
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

  gSystem->mkdir("maketof", true);
  std::string ofile = Form("maketof/Spline%04d.root", run);
  LoadTOF(run);
  GHistogramer::Get().SetOutFile(ofile);
  printf("Output histograms: %s\n", ofile.c_str());

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

  Unpacker unpacker;
  GBCS bcs;
  GCorrelator corr;
  std::vector<ddasHit> event;

  while(!converter.Finished() || !converter.Empty()) {
    if(converter.TryPop(event)) {
      ProcessEvent(unpacker, bcs, corr, event);
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

  PrintEventTypeSummary(); /// temporary

  corr.Flush();
  corr.Print();

  GHistogramer::Get().Close();

  // ---- build the TOF correction from the histogram we just wrote ----
  TDirectory* current = gDirectory;

  TFile* spline = TFile::Open(ofile.c_str());
  if(!spline || spline->IsZombie()) {
    printf(RED "Could not reopen %s - no TOF correction made\n" RESET_COLOR, ofile.c_str());
    return 1;
  }

  TH2* tof_time = dynamic_cast<TH2*>(spline->Get("tof"));
  if(!tof_time) {
    printf(RED "No 'tof' histogram in %s - no TOF correction made\n" RESET_COLOR, ofile.c_str());
    spline->Close();
    current->cd();
    return 1;
  }

  TOFCorrector corrector;
  corrector.MakeCorrectionFile(tof_time, run);

  spline->Close();
  current->cd();
  return 0;
}



void ProcessEvent(Unpacker& unpacker, GBCS& bcs, GCorrelator& corr, std::vector<ddasHit>& event) {

  unpacker.Unpack(bcs, event);   // Unpack() calls bcs.Reset() & bcs.classification()


if(!bcs.fPin1.HasHit() && !bcs.fPin2.HasHit()) return;

  double tof = bcs.fI2SPin1.fCharge;
  double runtime = bcs.fPin1.fTimestamp / 1.e8 ;
    
    GHistogramer::Get().Fill("tof",500,0,5000, runtime,
                                               1500,0,25000, tof);


}
