#ifndef GLABR_H
#define GLABR_H

#include <ddasHit.h>


class GLaBr {
 public:
  GLaBr();
  ~GLaBr();
 
  void Get(int Crystal, const ddasHit& hit); 
  void Reset();

 private: 

  double fCharge;
  double fEcal;
  double fTime;


  ClassDef(GLaBr,1);
};


#endif
