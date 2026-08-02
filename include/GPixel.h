#ifndef __GPIXEL_H__
#define __GPIXEL_H__

#include <vector>
#include <ddasHit.h>

class GPixel {
 public:
  GPixel();
  ~GPixel();

  void Reset();

  std::array<double,2> ImplantPixel  {-1,-1};
  std::array<double,2> DecayPixel    {-1,-1};

  double ImplantTime   {-1};
  double DecayTime     {-1};
  double ImplantEnergy {0};
  double DecayEnergy   {0};

  void BuildImplantPixel (const std::vector<ddasHit>& hit);
  void BuildDecayPixel (const std::vector<ddasHit>& hit);

  bool ImplantValid() const { return ImplantPixel[0] >= 0.0 && ImplantPixel[1] >= 0.0; }
  bool DecayValid()   const { return DecayPixel[0]   >= 0.0 && DecayPixel[1]   >= 0.0; }


 ClassDef(GPixel, 1);
};

#endif
