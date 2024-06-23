#include "cantera/thermo/IdealMoistPhase.h"
#include "cantera/base/utilities.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

double IdealMoistPhase::entropy_mole() const
{
  return GasConstant * mean_X(entropy_R_ref());
}

void IdealMoistPhase::setPressure(double p) {
  double g_ov_mu = std::accumulate(ym().begin(), ym().begin() + m_nvapor, 0.0);
  if (m_pressure_sets_temperature) {
    setTemperature(p / (g_ov_mu * GasConstant * density()));
  } else {
    setDensity(p / (GasConstant * temperature() * g_ov_mu));
  }
}

double IdealMoistPhase::pressure() const {
  double g_ov_mu = std::accumulate(ym().begin(), ym().begin() + m_nvapor, 0.0);
  return GasConstant * temperature() * density() * g_ov_mu;
}

}
