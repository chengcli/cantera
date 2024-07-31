#include "cantera/thermo/IdealMoistPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

double IdealMoistPhase::_g_ov_mu() const {
  return std::accumulate(ym().begin(), ym().begin() + nGas(), 0.0);
}

double IdealMoistPhase::_sum_xlogx_g() const {
  double sumxlogx = 0;
  double sumx = 0;
  for (size_t k = 0; k < nGas(); k++) {
    sumxlogx += ym()[k] * std::log(std::max(ym()[k], SmallNumber));
    sumx += ym()[k];
  }
  return sumxlogx / sumx + std::log(meanMolecularWeight());
}

void IdealMoistPhase::setPressure(double p) {
  if (p <= 0) {
    throw CanteraError("IdealGasPhase::setState_DP",
                       "pressure must be positive");
  }

  if (m_pressure_sets_temperature) {
    setTemperature(p / (_g_ov_mu() * GasConstant * density()));
  } else {
    setDensity(p / (GasConstant * temperature() * _g_ov_mu()));
  }
}

double IdealMoistPhase::pressure() const {
  return GasConstant * temperature() * density() * _g_ov_mu();
}

bool IdealMoistPhase::addSpecies(shared_ptr<Species> spec) {
  bool added = IdealGasPhase::addSpecies(spec);
  if (spec->input.hasKey("equation-of-state")) {
    auto& eos = spec->input["equation-of-state"].getMapWhere("model", "constant-volume");
    double mv;
    if (eos.hasKey("density")) {
      mv = molecularWeight(m_kk-1) / eos.convert("density", "kg/m^3");
    } else if (eos.hasKey("molar-density")) {
      mv = 1.0 / eos.convert("molar-density", "kmol/m^3");
    } else if (eos.hasKey("molar-volume")) {
      mv = eos.convert("molar-volume", "m^3/kmol");
    } else {
      throw CanteraError("IdealMoistPhase::addSpecies",
          "equation-of-state entry for species '{}' is missing "
          "'density', 'molar-volume', or 'molar-density' "
          "specification", spec->name);
    }
    m_speciesMolarVolume.push_back(mv);
    m_ncloud++;
  } else {
    m_speciesMolarVolume.push_back(0.);
  }
}

void IdealMoistPhase::getCv_R(double* cvr) const
{
  const vector<double>& _cpr = cp_R_ref();
  // gas
  for (size_t k = 0; k < nGas(); k++) {
    cvr[k] = _cpr[k] - 1.0;
  }
    
  // clouds
  copy(_cpr.begin() + nGas(), _cpr.end(), cvr + nGas());
}

void IdealMoistPhase::getIntEnergy_RT_ref(double* urt) const {
  const vector<double>& _h = enthalpy_RT_ref();

  // gas
  for (size_t k = 0; k < nGas(); k++) {
    urt[k] = _h[k] - 1.0;
  }

  // clouds
  copy(_h.begin() + nGas(), _h.end(), urt + nGas());
}

}
