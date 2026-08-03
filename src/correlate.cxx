#include <cstdio>
#include <deque>
#include <string>
#include <chrono>
#include <algorithm>
#include <filesystem>

#include <TChain.h>
#include <TString.h>
#include <TSystem.h>

#include <GHistogramer.h>
#include <globals.h>


static constexpr double fTicksPerMs  = 100000.0;            
static constexpr double fWindowMs    = 3000.0;     // 1 second window (yiyi defense)
static constexpr double fWindowTicks = fWindowMs * fTicksPerMs;

static constexpr int fXMin = 24, fXMax = 40;                
static constexpr int fYMin = 24, fYMax = 40;
static constexpr int fNeighbour = 1;

static bool ActiveRegion(int x, int y) {
  return x >= fXMin && x < fXMax &&
         y >= fYMin && y < fYMax;
}

struct Implant {
  double t;
  double energy;
  int isotope;
  bool   dead = false;
};

static std::deque<Implant> gGrid[40][40];

static void Prune(std::deque<Implant>& list, double now) {
  while(!list.empty() && now - list.front().t > fWindowTicks)
    list.pop_front();
}


int main (int argc, char** argv) {
  if(argc < 2) { printf("usage: correlate tree####.root\n"); return 1; }

  TChain chain("events");
  for(int i = 1; i < argc; ++i) chain.Add(argv[i]);

  int    type;
  int    xStrip; 
  int    yStrip;
  double timestamp; 
  double energy;
  int    isotope;
  std::vector<double>* abGammas = nullptr;

  chain.SetBranchAddress("type", &type);
  chain.SetBranchAddress("xStrip", &xStrip);
  chain.SetBranchAddress("yStrip", &yStrip);
  chain.SetBranchAddress("timestamp", &timestamp);
  chain.SetBranchAddress("energy", &energy);
  chain.SetBranchAddress("abGammas", &abGammas);
  chain.SetBranchAddress("isotope", &isotope);

  std::filesystem::path p(argv[1]);
  int run = -1;
  std::sscanf(p.stem().string().c_str(), "tree%d", &run);

  gSystem->mkdir("corl", true);
  GHistogramer::Get().SetOutFile(Form("corl/corl%04d.root", run));


  const Long64_t nEntries = chain.GetEntries();
  printf("correlating %lld entries   region x[%d,%d) y[%d,%d)   window %.0f ms   %dx%d search\n",
                         nEntries,    fXMin,fXMax,   fYMin,fYMax,   fWindowMs, 2*fNeighbour+1, 2*fNeighbour+1);

  long long nImp = 0, nDec = 0, nPairs = 0, nRejected = 0;

  auto     lastTime  = std::chrono::steady_clock::now();
  Long64_t lastEntry = 0;

  printf(HIDE_CURSOR);
  fflush(stdout);


  for(Long64_t i = 0; i < nEntries; ++i) {
    chain.GetEntry(i);

    if(ActiveRegion(xStrip, yStrip)) {

      if(type == 0) {                                  // ---- implant ----
        auto& list = gGrid[xStrip][yStrip];
        Prune(list, timestamp);
        list.push_back(Implant{timestamp, energy, isotope});
        ++nImp;

        GHistogramer::Get().Fill("corr/implant_map", 40, 0, 40, xStrip,
                                                     40, 0, 40, yStrip);
      }
      else {                                           // ---- decay ----
       ++nDec;

        GHistogramer::Get().Fill("corr/decay_map", 40, 0, 40, xStrip,
                                                   40, 0, 40, yStrip);
        GHistogramer::Get().Fill("corr/decay_energy", 4000, 0, 32000, energy);

        const int ixLo = std::max(0,  xStrip - fNeighbour);
        const int ixHi = std::min(39, xStrip + fNeighbour);
        const int iyLo = std::max(0,  yStrip - fNeighbour);
        const int iyHi = std::min(39, yStrip + fNeighbour);

        // prune the whole cluster first, so nothing shifts while we hold pointers
        for(int ix = ixLo; ix <= ixHi; ++ix)
          for(int iy = iyLo; iy <= iyHi; ++iy)
            Prune(gGrid[ix][iy], timestamp);

        // collect every live implant in the window
        std::vector<Implant*> candidates;
        for(int ix = ixLo; ix <= ixHi; ++ix) {
          for(int iy = iyLo; iy <= iyHi; ++iy) {
            for(auto& im : gGrid[ix][iy]) {
              if(im.dead) continue;
              const double dt = timestamp - im.t;
              if(dt > 0 && dt <= fWindowTicks)
                candidates.push_back(&im);
            }
          }
        }

        const int nCand = candidates.size();
       // GHistogramer::Get().Fill("time/n_candidates", 10, 0, 10, nCand);

        if(nCand == 1) {                               // unambiguous - accept
          const double dt = timestamp - candidates[0]->t;
          GHistogramer::Get().Fill("time/dt_vs_count", 3000, 0, fWindowMs, dt / fTicksPerMs);
          GHistogramer::Get().Fill("time/dt_vs_denergy", 3000, 0, fWindowMs, dt / fTicksPerMs,
                                                         4000, 0, 32000, energy);   // decay energy
          if(abGammas) { for(double e : *abGammas) {
          GHistogramer::Get().Fill("gamma/coinc_addback", 8000, 0, 8000, e);
          GHistogramer::Get().Fill("gamma/dt_vs_gamma", 1000, 0, fWindowMs, dt / fTicksPerMs,
                                                          8000, 0, 8000,     e);}
          }

          const int Isotope = candidates[0]->isotope;
           if(Isotope != -1) {
             if(Isotope == 0) {
              GHistogramer::Get().Fill("decay_time/Na_dt_vs_count",3000,0, fWindowMs, dt / fTicksPerMs);
            }
            if(Isotope == 1) {
              GHistogramer::Get().Fill("decay_time/Ne_dt_vs_count",3000,0, fWindowMs, dt / fTicksPerMs);
            }
            if(Isotope == 2) {
              GHistogramer::Get().Fill("decay_time/F_dt_vs_count",3000,0, fWindowMs, dt / fTicksPerMs);
            }
          }    
         ++nPairs;
        }
        else if(nCand > 1) {                           // back-to-back implants - reject
          for(auto* im : candidates) im->dead = true;
          ++nRejected;
        }
      }  
    }

 if(i % 100000 == 0) {
      auto now = std::chrono::steady_clock::now();
      double elapsed = std::chrono::duration<double>(now - lastTime).count();
      double rate = elapsed > 0 ? (i - lastEntry) / elapsed : 0;

      constexpr int barWidth = 30;
      double percent = 100.0 * i / nEntries;
      int filled = std::clamp(static_cast<int>((percent / 100.0) * barWidth), 0, barWidth);
      std::string bar(filled, '#');
      bar.append(barWidth - filled, '-');

      printf(CLEAR_LINE "[%s] %6.2f%%  %.2f M/s  imp=%lld dec=%lld pairs=%lld",
             bar.c_str(), percent, rate / 1e6, nImp, nDec, nPairs);
      fflush(stdout);

      lastTime = now; lastEntry = i;
    }
  }

  printf(CLEAR_LINE "[%s] 100.00%%\n" SHOW_CURSOR, std::string(30, '#').c_str());
  printf("implants %lld   decays %lld   pairs %lld  rejected %lld\n", nImp, nDec, nPairs, nRejected);

  GHistogramer::Get().Close();
  printf("wrote corr/corr%04d.root\n", run);
  return 0;
}



