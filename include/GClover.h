#ifndef GCLOVER_H
#define GCLOVER_H

#include <ddasHit.h>
#include <vector>


class GClover {
 public:
  GClover();
  ~GClover();

  void Reset();
  void Unpack(const ddasHit& hit, int id);

  bool HasHit() const;
 

  int    fId;
  double fCharge;
  double fEcal;
  double fTimestamp;

 ClassDef(GClover,1);
};


#endif
