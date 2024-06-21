#ifndef CT_IDEALGASSURF_H
#define CT_IDEALGASSURF_H

#include "IdealGasPhase.h"
#include "SurfPhase.h"

namespace Cantera
{

class IdealGasSurfPhase: public SurfPhase
{
 public:
  explicit IdealGasSurfPhase(const string& inputFile="", const string& id=""):
    SurfPhase(inputFile, id) 
  {}

  bool isCompressible() const override {
    return true;
  }

  string type() const override {
    return "ideal-gas-surface";
  }

  string phaseOfMatter() const override {
    return "gas";
  }
};

}

#endif
