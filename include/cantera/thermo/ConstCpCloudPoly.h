#pragma once

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/speciesThermoTypes.h"

namespace Cantera
{

class Kinetics;
class ReactionRate;

/**
 * This thermodynamic property manager class uses vapor phase as
 * the reference state for calculating the entropy of a condensed phase.
 */
class ConstCpCloudPoly: public SpeciesThermoInterpType
{
 public:
  ConstCpCloudPoly();

  ConstCpCloudPoly(double tlow, double thigh, double pref, 
                   const double* coeffs);

  void setParameters(double t0, double h0, double cp0);
  void setReaction(std::string name, Kinetics& kin);

  int reportType() const override {
    return CONSTANT_CP_CLOUD;
  }

  void updateProperties(const double* tt, double* cp_R, double* h_RT,
                        double* s_R) const override {
    updatePropertiesTemp(*tt, cp_R, h_RT, s_R);
  }

  void updatePropertiesTemp(const double temp, double* cp_R, double* h_RT,
                            double* s_R) const override;

  size_t nCoeffs() const override { return 3; }

  void reportParameters(size_t& n, int& type, double& tlow, double& thigh,
                        double& pref, double* const coeffs) const override;

  void getParameters(AnyMap& thermo) const override;

 protected:
  //! Base temperature
  double m_t0 = 298.15;

  //! Dimensionless value of the heat capacity
  double m_cp0_R = 0.0;

  //! dimensionless value of the enthalpy at t0
  double m_h0_R = 0.0;

  std::shared_ptr<ReactionRate> m_rate;

  std::vector<std::shared_ptr<SpeciesThermoInterpType>> m_vapor_thermo;
};

}
