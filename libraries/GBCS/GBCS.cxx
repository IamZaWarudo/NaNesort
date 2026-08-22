

#include <ddasHit.h>
#include <Unpacker.h>
#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

#include <globals.h>

GBCS::GBCS() { }

void GBCS::Reset() {

 // DSSDLow.Reset();  // Not using this - unreliable
 fDSSD.Reset();

 fPin1.Reset();
 fPin2.Reset();
  
 fI2SPin1.Reset();
 fI2SPin2.Reset();

 fSSSDLow.Reset();
 fSSSDHigh.Reset();

}





int GBCS::EventType() const{

int condition = 4;  // Unidentified

//  fEventType =  condition;
 
 return condition;
}




/*
void GBCS::Print(int type) const {
  switch(fEventType) {
    case 1: //implant
      printf(DRED);
      break;
    case 2: //decay
      printf(DGREEN);
      break;
    case 3: //veto
      printf(DYELLOW);
      break;
    case 4: //unknown
    default:
      printf(CYAN);
      break;
  }
  
  if(type<0 || type==fEventType) {
    printf("BCS [%02i][%02i] @ %.1f\n",X(),Y(),Timestamp());
    printf("EvnetType = %i\n",fEventType);
    printf("\tPin1: %f\n",Pin1.fEcal);
    printf("\tTOF:  %f\n",TOF());
    printf("\tDSSD MULTI [%i][%i]\n",DSSD.fxMult,DSSD.fyMult);
    printf("\tDSSD charege %.1f\n",DSSD.fxE + DSSD.fyE);
    printf("\tSSSDLow MULTI [%lu]\n",SSSDLow.fStrips.size());
    printf("\tSSSDHigh MULTI [%lu]\n",SSSDHigh.fStrips.size());
  }
  printf(RESET_COLOR);
  return;
}   */
