#ifndef __GBCS_H__
#define __GBCS_H__

#include <vector>
#include <Unpacker.h>
#include <ddasHit.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>
#include <GCloverHit.h>

/****************
Implant -----> 1
*****************
Decay -------> 2
*****************
Veto ----> 3
*****************
Undentified -> 4
*****************/


class GBCS {
 public:
  GBCS();
  ~GBCS() = default;

  void Reset();


  //GDSSD DSSDLow;

  double Timestamp() const { return DSSD.fTimestamp; }
  double TOF() const { return I2SPin1.fEcal>0 ? I2SPin1.fEcal : -1; } 
  int    X()   const { return DSSD.Xpos; }
  int    Y()   const { return DSSD.Ypos; }

  void Print(int type=-1) const;
  int EventType() const { return fEventType; }


//private:
  GDSSD DSSD;  // this DSSD is High Gain since the low gain seemed to be acting weird
  GPIN Pin1;
  GPIN Pin2;
  GPIN I2SPin1;
  GPIN I2SPin2;

  GSSSD SSSDLow;
  GSSSD SSSDHigh;

  std::vector<GCloverHit> fCloverHits;

  int fEventType;

  void Fill(const Unpacker& Event);

};




#endif
