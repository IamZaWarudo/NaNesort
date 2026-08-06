#ifndef __GPIXEL_H__
#define __GPIXEL_H__

#include <array>
#include <vector>
#include <ddasHit.h>

class GPixel {
 public:
  GPixel();
  ~GPixel();

  void Reset();

  void BuildImplantPixel (const std::vector<ddasHit>& hits);
  void BuildDecayPixel   (const std::vector<ddasHit>& hits);

  bool Implant() const { return ImPixel[0] >= 0.0 && ImPixel[1] >= 0.0; }
  bool Decay()   const { return DePixel[0] >= 0.0 && DePixel[1] >= 0.0; }

  std::array<double,2> ImPixel = {-1.0, -1.0};
  std::array<double,2> DePixel = {-1.0, -1.0};

  double ImTime = -1.0;
  double ImEnergy = 0.0;
  int ImfMaxStrip = -1;
  int ImbMaxStrip = -1;
  int ImfMult = 0;
  int ImbMult = 0;

  double DeTime = -1.0;
  double DeEnergy = 0.0;
  int DefMaxStrip = -1;
  int DebMaxStrip = -1;
  int DefMult = 0;
  int DebMult = 0;


 ClassDef(GPixel, 1);
};

#endif
