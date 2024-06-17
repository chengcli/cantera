#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Condensation.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

void Condensation::resizeReactions()
{
  Kinetics::resizeReactions();

  for (auto& rates : m_interfaceRates) {
    rates->resize(nTotalSpecies(), nReactions(), nPhases());
  }
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

  // multiply ropf by the activity concentration reaction orders to obtain
  // the forward rates of progress.
  m_reactantStoich.multiply(m_actConc.data(), m_ropf.data());

  // products
  m_revProductStoich.multiply(m_actConc.data(), m_ropr.data());

  for (size_t j = 0; j != nReactions(); ++j) {
    m_ropnet[j] = m_ropf[j];
    m_ropr[j] = 0.;
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
  for (size_t n = 0; n < nPhases(); n++) {
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
