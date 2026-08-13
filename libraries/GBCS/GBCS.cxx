

#include <ddasHit.h>
#include <Unpacker.h>
#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

GBCS::GBCS() { }

void GBCS::Reset() {

 DSSDLow.Reset();
 DSSDHigh.Reset();

 Pin1.Reset();
 Pin2.Reset();

 SSSDLow.Reset();
 SSSDHigh.Reset();
  
}


void GBCS::Fill(const Unpacker& Event) {
  

 DSSDLow = Event.fDSSDLow;
 DSSDHigh = Event.fDSSDHigh;

 Pin1 = Event.fPin1;
 Pin2 = Event.fPin2;
 I2SPin1 = Event.fI2SPin1;
 I2SPin2 = Event.fI2SPin2;


 SSSDLow = Event.fSSSDLow;
 SSSDHigh = Event.fSSSDHigh;

}


/****************************************
 * Correlation flags:                   *
 * ==================================== *
 * 4  - implant good position    -------------*
 * 8  - implant bad position            *
 * 12 - decay good position      -------------*
 * 16 - decay bad position              *
 * 20 - light ion good position  -------------*
 * 24 - light ion bad position
 * 28 - decay before ion (BAD!)         *
 * 32 - good correlation                *
 * 44 - no implants to correlate        *
 * 56 - Clover only     *
 * 96 - unknown event                   *
 * ==================================== *
 * 99 - Reset correlation aray          *
 ****************************************/


/*
Implant - 1
decay   - 2
LightIon - 3
Undentified - 4
*/

int GBCS::EventType(){

 int condition = 123456789;

//signal flags
 bool DSSDLowHasGoodPosition = false;
 bool DSSDHighHasGoodPosition = false;
 
 bool Pin1HasHit = false;
 bool Pin2HasHit = false;

 bool HasPinImplant = false;
 bool HasPinLightIon = false;

 bool Pin1HasImplant = false;
 bool Pin1HasLightIon = false;
 bool Pin2HasImplant = false;
 bool Pin2HasLightIon = false;

 bool SSSDLowHasHit = false;
 bool SSSDHighHasHit = false;

 bool HasTOF = false;

 bool HasClover = false;
 bool HasLaBr = false;


//event identification
 bool IsImplant = false;        // 1
 bool IsDecay = false;          // 2
 bool IsLightIon = false;       // 3  
 bool IsUnidentified = false;   // 4


// add ifs here
 
  if(DSSDLow.HasGoodPosition()){
   DSSDLowHasGoodPosition = true;
  }
  
  if(DSSDHigh.HasGoodPosition()){
   DSSDHighHasGoodPosition = true;
  }

  if(Pin1.HasHit()){
    Pin1HasHit = true;
  }
  
  if(Pin2.HasHit()){
    Pin2HasHit = true;
  }


  if(SSSDLow.HasHit()){
    SSSDLowHasHit = true;
  }
  
  if(SSSDHigh.HasHit()){
    SSSDHighHasHit = true;
  }


if(Pin1HasHit && Pin2HasHit && DSSDLowHasGoodPosition && !SSSDLowHasHit && !SSSDHighHasHit){
 IsImplant = true;
 condition = 1;
} else
if(!Pin1HasHit && !Pin2HasHit && !DSSDLowHasGoodPosition && DSSDHighHasGoodPosition && !SSSDLowHasHit && !SSSDHighHasHit){
 IsDecay = true;
 condition = 2;
} else
if(Pin1HasHit && Pin2HasHit && !DSSDLowHasGoodPosition && DSSDHighHasGoodPosition && !SSSDLowHasHit && SSSDHighHasHit){
 IsLightIon = true;
 condition = 3;
}

if(!IsImplant && !IsDecay && !IsLightIon) { 
 IsUnidentified = true;
 condition = 4;
}


  return condition;
}
