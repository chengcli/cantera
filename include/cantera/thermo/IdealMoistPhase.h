#ifndef CT_IDEALMOISTPHASE_H
#define CT_IDEALMOISTPHASE_H

#include "IdealGasPhase.h"

namespace Cantera
{

class Kinetics;
class ReactionRate;

class IdealMoistPhase : public IdealGasPhase {
 public:
  explicit IdealMoistPhase(const std::string& infile = "", const std::string& id=""):
    IdealGasPhase(infile, id) {}

  string type() const override {
    return "ideal-moist"; 
  }

  void setTemperature(double temp) override {
    IdealGasPhase::setTemperature(temp);
    m_pressure_sets_temperature = false;
  }

  void setDensity(double dens) override {
    IdealGasPhase::setDensity(dens);
    m_pressure_sets_temperature = true;
  }

  /*!
   * @f[
   * \hat s(T, P) = \sum_k X_k \hat s^0_k(T) - \hat R \ln \frac{P}{P^0}.
   *                - \sum_{k \elem G} X_k \ln X_k
   * @f]
   */
  double entropy_mole() const override;

  /*!
   * @f[ 
   * \hat c_v = \hat c_p - g * \hat R. 
   * @f]
   */
  double cv_mole() const override {
    return cp_mole() - _g_ov_mu() * meanMolecularWeight() * GasConstant;
  }

  double pressure() const override;

  void setPressure(double p) override;

  bool addSpecies(shared_ptr<Species> spec) override;

  void setState_DP(double rho, double p) override {
    setDensity(rho);
    setPressure(p);
  }

  void getChemPotentials(double* mu) const override {
    getStandardChemPotentials(mu);
  }

  void getCv_R(double* cvr) const override;

  void getIntEnergy_RT_ref(double* urt) const override;

  size_t nGas() const {
    return nSpecies() - m_ncloud;
  }

  void updateFromKinetics(Kinetics& kin) override;

 protected:
  /*!
   * @f[
   * \frac{\sum_{i \elem G} X_i}{\bar{\mu}}
   * @f]
   */
  double _g_ov_mu() const;

  /*!
   * @f[
   * \frac{\sum_{i \elem G} X_i \ln X_i}{\sum_{i \elem G} X_i}
   * @f]
   */
  double _sum_xlogx_g() const;

  bool  m_pressure_sets_temperature = false;

  //! Number of cloud species
  size_t m_ncloud = 0;

  //! Index of vapor species forming clouds
  std::vector<std::vector<int>> m_vapor_index;

  //! Vector of reaction rates for evaluating cloud entropy
  std::vector<std::shared_ptr<ReactionRate>> m_rate;

  //! Vector of molar volumes for each species in the solution
  /**
   * Species molar volumes (@f$ m^3 kmol^-1 @f$) at the current mixture state.
   */
  mutable vector<double> m_speciesMolarVolume;

  void updateThermo() const override;
};

}

#endif
