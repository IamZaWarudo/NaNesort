#ifndef __UNPACKER_H__
#define __UNPACKER_H__

#include <vector>

#include <ddasHit.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>
#include <GCloverHit.h>

class Unpacker {
 public:
  Unpacker();
  ~Unpacker() = default;

  void Reset();
  void Unpack (const std::vector<ddasHit>& event);

// all the detectors 
//  GPIN fI2N;  // this only intermitently worked, so largely useless 
  GPIN fI2SPin1;   
  GPIN fI2SPin2;
  GPIN fI2TAC; // I2S-I2N
  GPIN fPin1;
  GPIN fPin2;
//  GPIN fPin3;

  GDSSD fDSSDLow;
  GDSSD fDSSDHigh;

  GSSSD fSSSDLow;  // Veto Particles
  GSSSD fSSSDHigh; // gamma spectrum

  GLaBr fLaBr;
  GClover fCrystal;

  std::vector<GCloverHit> fCloverHits;

};


#endif
