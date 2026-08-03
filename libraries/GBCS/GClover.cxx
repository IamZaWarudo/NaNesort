#include <GClover.h>

GClover::GClover()  { Reset(); }
GClover::~GClover() {}

void GClover::Reset() {
  hits.clear();
  addbackHits.clear();
}

void GClover::Get(int crystal, const ddasHit& hit) {
  if(crystal < 0 || crystal > 63) return;
  if(hit.GetEcal() <= 0.0)        return;

  GCloverHit h;
  h.fId     = crystal;
  h.fCharge = hit.GetCharge();
  h.fEcal   = hit.GetEcal();
  h.fTime   = hit.GetTime();
  hits.push_back(h);
}

void GClover::BuildAddback() {
  addbackHits.clear();

  double eSum[16] = {0};
  double qSum[16] = {0};
  double maxE[16] = {0};      // local only - picks whose time to keep
  double tMax[16] = {0};
  bool   fired[16] = {false};

  for(const auto& h : hits) {
    const int c = h.fId / 4;
    if(c < 0 || c >= 16) continue;

    eSum[c] += h.fEcal;
    qSum[c] += h.fCharge;
    fired[c] = true;

    if(h.fEcal > maxE[c]) { maxE[c] = h.fEcal; tMax[c] = h.fTime; }
  }

  for(int c = 0; c < 16; ++c) {
    if(!fired[c]) continue;
    GCloverHit ab;
    ab.fId     = c;
    ab.fCharge = qSum[c];
    ab.fEcal   = eSum[c];
    ab.fTime   = tMax[c];
    addbackHits.push_back(ab);
  }
}
