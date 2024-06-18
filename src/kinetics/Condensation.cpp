// C/C++
#include <numeric>
#include <algorithm>

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Condensation.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

// This works for these two types of reactions:
// (1) A -> B (order 1)
// (2) A + B -> C (order 2)
inline double saturation_function(double delta, double react, double prod, size_t order)
{
  if (delta < 0) {
    return -std::min(prod, -delta);
  }

  if (order == 1) {
    return delta;
  } else if (order == 2) {
    return (react - sqrt(react * react - 4 * delta)) / 2.;
  }

  return 0.;
}

void Condensation::resizeReactions()
{
  Kinetics::resizeReactions();

  for (auto& rates : m_interfaceRates) {
    rates->resize(nTotalSpecies(), nReactions(), nPhases());
  }
  m_satf.resize(nReactions());
}

void Condensation::getActivityConcentrations(double* const conc)
{
  _update_rates_C();
  copy(m_actConc.begin(), m_actConc.end(), conc);
}

void Condensation::getFwdRateConstants(double* kfwd)
{
  updateROP();
  copy(m_rfn.begin(), m_rfn.end(), kfwd);
}

void Condensation::resizeSpecies()
{
  size_t kOld = m_kk;
  Kinetics::resizeSpecies();
  if (m_kk != kOld && nReactions()) {
      throw CanteraError("InterfaceKinetics::resizeSpecies", "Cannot add"
          " species to InterfaceKinetics after reactions have been added.");
  }

  m_actConc.resize(m_kk);
}

bool Condensation::addReaction(shared_ptr<Reaction> r_base, bool resize)
{
  size_t i = nReactions();
  bool added = Kinetics::addReaction(r_base, resize);
  if (!added) {
      return false;
  }

  // Set index of rate to number of reaction within kinetics
  shared_ptr<ReactionRate> rate = r_base->rate();
  rate->setRateIndex(nReactions() - 1);
  rate->setContext(*r_base, *this);

  string rtype = rate->subType();
  if (rtype == "") {
      rtype = rate->type();
  }

  // If necessary, add new interface MultiRate evaluator
  if (m_interfaceTypes.find(rtype) == m_interfaceTypes.end()) {
    m_interfaceTypes[rtype] = m_interfaceRates.size();
    m_interfaceRates.push_back(rate->newMultiRate());
    m_interfaceRates.back()->resize(m_kk, nReactions(), nPhases());
  }

  // Add reaction rate to evaluator
  size_t index = m_interfaceTypes[rtype];
  m_interfaceRates[index]->add(nReactions() - 1, *rate);

  return true;
}

void Condensation::updateROP() {
  // evaluate rate constants and equilibrium constants at temperature and phi
  // (electric potential)
  _update_rates_T();
  // get updated activities (rates updated below)
  _update_rates_C();

  if (m_ROP_ok) {
    return;
  }

  std::fill(m_ropf.begin(), m_ropf.end(), 1.0);

  for (size_t i = 0; i < m_actConc.size(); i++) {
    std::cout << "actConc[" << i << "] = " << m_actConc[i] << std::endl;
  }

  // multiply ropf by the activity concentration reaction orders to obtain
  // the forward rates of progress.
  m_reactantStoich.multiply(m_actConc.data(), m_ropf.data());

  for (size_t j = 0; j != nReactions(); ++j) {
    std::cout << "ropf[" << j << "] = " << m_ropf[j] << std::endl;
  }

  // m_rfn -> saturation activity concentration
  // m_ropf -> current activity concentration

  for (int j = 0; j < nReactions(); j++) {
    // inactive reactions
    if (m_rfn[j] < 0.0) {
      m_satf[j] = 0.0;
      continue;
    }

    // calculate saturation function
    auto& R = m_reactions[j];
    double react = std::accumulate(R->reactants.begin(), R->reactants.end(), 0.0,
                    [&](double sum, const std::pair<std::string, double>& r) {
                      return sum + m_actConc[kineticsSpeciesIndex(r.first)] / r.second;
                    });
    double prod = std::accumulate(R->products.begin(), R->products.end(), 
                    std::numeric_limits<double>::infinity(),
                    [&](double min, const std::pair<std::string, double>& p) {
                      return std::min(min, m_actConc[kineticsSpeciesIndex(p.first)] / p.second);
                    });
    m_satf[j] = saturation_function(m_ropf[j] - m_rfn[j], react, prod, R->reactants.size());
  }

  for (size_t j = 0; j != nReactions(); ++j) {
    m_ropf[j] = std::max(0., m_satf[j]);
    m_ropr[j] = std::max(0., -m_satf[j]);
    m_ropnet[j] = m_ropf[j] - m_ropr[j];
  }
}

void Condensation::_update_rates_T()
{
  // Go find the temperature from the surface
  double T = thermo(0).temperature();

  if (T != m_temp) {
    m_temp = T;
    m_ROP_ok = false;
  }

  // loop over interface MultiRate evaluators for each reaction type
  for (auto& rates : m_interfaceRates) {
    bool changed = rates->update(thermo(0), *this);
    if (changed) {
      rates->getRateConstants(m_rfn.data());
      m_ROP_ok = false;
    }
  }
}

void Condensation::_update_rates_C()
{
  // dry air has activity concentration of 0.0
  m_actConc[0] = 0.0;

  for (size_t n = 1; n < nPhases(); n++) {
    const auto& tp = thermo(n);
    /*
     * We call the getActivityConcentrations function of each ThermoPhase
     * class that makes up this kinetics object to obtain the generalized
     * concentrations for species within that class. This is collected in
     * the vector m_conc. m_start[] are integer indices for that vector
     * denoting the start of the species for each phase.
     */
    tp.getActivityConcentrations(m_actConc.data() + m_start[n]);
  }
  m_ROP_ok = false;
}


}
