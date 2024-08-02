#include "cantera/thermo/IdealMoistPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/utilities.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/Reaction.h"

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

double IdealMoistPhase::entropy_mole() const {
  double xg = _g_ov_mu() * meanMolecularWeight();
  return GasConstant * (mean_X(entropy_R_ref()) - _sum_xlogx_g()
           - xg * std::log(pressure() / (xg * refPressure())));
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

void IdealMoistPhase::updateFromKinetics(Kinetics& kin) 
{
  Composition reactants;
  m_vapor_index.resize(m_ncloud);

  for (size_t i = 0; i < m_ncloud; i++) {
    auto name = speciesName(m_kk - m_ncloud + i);
    for (size_t j = 0; j < kin.nReactions(); j++) {
      auto rxn = kin.reaction(j);
      if (rxn->products.size() == 1 && rxn->products.begin()->first == name) 
      {
        m_rate.push_back(rxn->rate());
        reactants = rxn->reactants;
        break;
      }
    }

    if (reactants.empty()) {
      throw CanteraError("ConstCpCloudPoly::setReaction",
                         "No reaction found for species " + name);
    }

    for (auto& [name, _] : reactants) {
      m_vapor_index[i].push_back(speciesIndex(name));
    }
  }
}

void IdealMoistPhase::updateThermo() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tnow) {
        m_spthermo.update(tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        cached.state1 = tnow;

        // revise cloud entropy
        for (size_t i = 0; i < m_ncloud; i++) {
          int j = m_kk - m_ncloud + i;
          int order = m_vapor_index[i].size();
          auto svp = m_rate[i]->eval(tnow) * pow(GasConstant * tnow, order);
          auto pref = m_spthermo.getSpeciesThermo(j)->refPressure();
          m_s0_R[j] = - log(svp / pref) + m_h0_RT[j];

          /*std::cout << "svp = " << svp << std::endl;
          std::cout << "pref = " << pref << std::endl;*/
          for (auto k : m_vapor_index[i]) {
            m_s0_R[j] += m_s0_R[k] - m_h0_RT[k];
            /*std::cout << "k = " << k << ", "
                      << "sv_R = " << m_s0_R[k] << ", "
                      << "hv_RT = " << m_h0_RT[k] << std::endl;*/
          }
          //std::cout << std::endl;
        }

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
}

}
