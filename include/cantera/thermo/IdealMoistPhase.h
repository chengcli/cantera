#ifndef CT_IDEALMOISTPHASE_H
#define CT_IDEALMOISTPHASE_H

#include "IdealGasPhase.h"

namespace Cantera
{

class IdealMoistPhase : public IdealGasPhase {
 public:
  explicit IdealMoistPhase(const std::string& infile = "", const std::string& id=""):
    IdealGasPhase(infile, id) {}

  string type() const override {
    return "ideal-moist"; 
  }

  string phaseOfMatter() const override {
    return "gas"; 
  }

  void setTemperature(double temp) override {
    IdealGasPhase::setTemperature(temp);
    m_pressure_sets_temperature = false;
  }

  void setDensity(double dens) override {
    IdealGasPhase::setDensity(dens);
    m_pressure_sets_temperature = true;
  }

  double entropy_mole() const override;

  double pressure() const override;

  void setPressure(double p) override;

  void getChemPotentials(double* mu) const override {
    getStandardChemPotentials(mu);
  }

  size_t nVapor() const {
    return m_nvapor;
  }

  void setNumVapors(size_t nvapors) {
    m_nvapor = nvapors;
  }

 protected:
  bool  m_pressure_sets_temperature = false;
  size_t m_nvapor = 0;
};

}

#endif
