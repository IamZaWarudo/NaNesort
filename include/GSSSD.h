#ifndef GSSSD_H
#define GSSSD_H

#include <ddasHit.h>


class GSSSD {
  public:
  GSSSD();
  ~GSSSD();
 
  void Unpack(const ddasHit& hit);
  void Reset();

//  private:

  bool HasHit();
  bool Hit = false;

  double fEcal;
  double Timestamp;
  int strip;

  ClassDef(GSSSD,1);
};


#endif
