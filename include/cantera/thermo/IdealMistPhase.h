#ifndef CT_IDEALMISTPHASE_H
#define CT_IDEALMISTPHASE_H

#include "IdealGasPhase.h"

namespace Cantera
{

class IdealMistPhase : public IdealGasPhase {
 public:
  explicit IdealMistPhase(const std::string& infile = "", const std::string& id=""):
    IdealGasPhase(infile, id) {}

  string type() const override {
    return "ideal-mist"; 
  }

  string phaseOfMatter() const override {
    return "condensed"; 
  }

  /**
   * Molar entropy. Units: J/kmol/K.
   * For an ideal mist mixture,
   * @f[
   * \hat s(T, P) = \sum_k X_k \hat s^0_k(T)
   * @f]
   * The reference-state pure-species entropies @f$ \hat s^0_k(T) @f$ are
   * computed by the species thermodynamic property manager.
   * @see MultiSpeciesThermo
   */
  double entropy_mole() const override;

  /**
   * Pressure. Units: Pa.
   * For an ideal gas mixture,
   * @f[ P = n \hat R T. @f]
   */
  double pressure() const override {
      return 0.;
  }

  double isothermalCompressibility() const override { 
    return 0.; 
  }

  double thermalExpansionCoeff() const override { 
    return 0.; 
  }

  void getChemPotentials(double* mu) const override {
    getStandardChemPotentials(mu);
  }

  void getPartialMolarEntropies(double* sbar) const override;
};

}

#endif
