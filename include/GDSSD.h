#ifndef GDSSD_H
#define GDSSD_H

#include <ddasHit.h>

class GDSSD {
  public:
  GDSSD();
  ~GDSSD();
  
  void UnpackFront(const ddasHit& hit);
  void UnpackBack(const ddasHit& hit);
  
  void Reset();

  bool HasGoodPosition();

//  private:

  int     fId;
  double  fEcal;      
  double  frontTimestamp;
  double  backTimestamp;
  int     xpos;       
  int     ypos;       
  int     Xpos;       
  int     Ypos;       
  bool    GoodX;      
  bool    GoodY;      
  
  double  fxEsum;     
  double  fxE;        
  double  fxMaxE;     
  double  fxMaxStrip; 
  double  fyEsum;     
  double  fyE;        
  double  fyMaxE;     
  double  fyMaxStrip; 


 ClassDef(GDSSD,1);
};


#endif
