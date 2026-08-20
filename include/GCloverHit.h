#ifndef GCLOVERHIT_H
#define GCLOVERHIT_H

#include <ddasHit.h>
#include <vector>


class GCloverHit {
 public:
  GCloverHit();
  ~GCloverHit();

  void Reset();
  void Unpack(const ddasHit& hit, int id);

  //bool HasHit() const;
 
  int    fId;   //should be 0-63?
  double fCharge;
  double fEcal;
  double fTimestamp;

 ClassDef(GCloverHit,1);
};


#endif
