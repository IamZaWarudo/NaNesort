#ifndef GLABR_H
#define GLABR_H

#include <ddasHit.h>


class GLaBr {
 public:
  GLaBr();
  ~GLaBr();
 
  void Unpack(const ddasHit& hit); 
  void Reset();

 private: 

  double fEcal;
  double fCharge;
  double fTimestamp;


  ClassDef(GLaBr,1);
};


#endif
