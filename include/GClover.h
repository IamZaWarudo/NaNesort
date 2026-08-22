#ifndef GCLOVER_H
#define GCLOVER_H

#include <ddasHit.h>
#include <vector>


class GClover {
 public:
  GClover();
  ~GClover();

  std::vector<ddasHit> fCloverHits;

  void Reset();
  void Unpack(const ddasHit& hit);
  bool HasHit() const;
 
  void AddBack();

  int    fId;
  double fCharge;
  double fEcal;
  double fTimestamp;

 ClassDef(GClover,1);
};


#endif
