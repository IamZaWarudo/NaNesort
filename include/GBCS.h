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
  GDSSD DSSD;  // this DSSD is High Gain since the low gain seemed to be acting weird

  GPIN Pin1;
  GPIN Pin2;
  GPIN I2SPin1;
  GPIN I2SPin2;

  GSSSD SSSDLow;
  GSSSD SSSDHigh;

  int EventType();

  void Fill(const Unpacker& Event);

};




#endif
