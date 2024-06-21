// C/C++
#include <numeric>
#include <algorithm>

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Condensation.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

// Only works for A -> B (order 1)
inline double saturation_function1(double ss, double react, double prod)
{
  if (ss < 0) {
    return -std::min(prod, -ss);
  }
  return ss;
}

// Only works for A + B -> C (order 2)
inline double saturation_function2(double ss, double react, double prod)
{
  if (ss < 0.) {
    return -std::min(prod, (- react + sqrt(react * react - 4 * ss)) / 2.);
  }
  return (react - sqrt(react * react - 4 * ss)) / 2.;
}

void Condensation::resizeReactions()
{
  Kinetics::resizeReactions();
  m_rbuf0.resize(nReactions());
  m_rbuf1.resize(nReactions());

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
  auto jac = netRatesOfProgress_ddCi();

  if (m_ROP_ok) {
    return;
  }

  std::fill(m_rbuf0.begin(), m_rbuf0.end(), 1.0);

  /*for (size_t i = 0; i < m_actConc.size(); i++) {
    std::cout << "actConc[" << i << "] = " << m_actConc[i] << std::endl;
  }*/

  // multiply ropf by the activity concentration reaction orders to obtain
  // the forward rates of progress.
  m_reactantStoich.multiply(m_actConc.data(), m_rbuf0.data());

  /*for (size_t j = 0; j != nReactions(); ++j) {
    std::cout << "rbuf[" << j << "] = " << m_rbuf0[j] << std::endl;
  }*/

  // m_rfn -> saturation activity 
  // m_rbuf -> current activity 

  Eigen::VectorXd b(nReactions());
  Eigen::SparseMatrix<double> stoich(m_stoichMatrix);

  for (int j = 0; j < nReactions(); j++) {
    // inactive reactions
    if (m_rfn[j] < 0.0) {
      b(j) = 0.0;
      for (int i = 0; i < nTotalSpecies(); i++) {
        stoich.coeffRef(i,j) = 0.0;
      }
      continue;
    }

    // calculate saturation function
    double ss = m_rbuf0[j] - m_rfn[j];
    size_t order = m_reactions[j]->reactants.size();
    auto& R = m_reactions[j];

    if (order == 1) { // order 1
      size_t iy = kineticsSpeciesIndex(R->products.begin()->first);
      double y = m_actConc[iy];
      b(j) = saturation_function1(ss, 0. /* dummy */, y);
    } else { // order 2
      size_t ix1 = kineticsSpeciesIndex(R->reactants.begin()->first);
      size_t ix2 = kineticsSpeciesIndex(next(R->reactants.begin())->first);
      size_t iy = kineticsSpeciesIndex(R->products.begin()->first);
      double x1 = m_actConc[ix1];
      double x2 = m_actConc[ix2];
      double y = m_actConc[iy];
      b(j) = saturation_function2(ss, x1 + x2, y);
    }
  }

  //auto A = jac * m_stoichMatrix;
  auto A = jac * stoich;

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
  solver.compute(A);
  Eigen::VectorXd r = solver.solve(b);

  if (solver.info() != Eigen::Success) {
    throw CanteraError("Condensation::updateROP",
                       "Failed to solve for net rates of progress.");
  }

  /*std::cout << jac << std::endl;
  std::cout << A << std::endl;
  std::cout << A.transpose() * A << std::endl;
  std::cout << "r = " << -r << std::endl;
  std::cout << "S.r = " << -m_stoichMatrix * r << std::endl;*/
  std::cout << "b = " << b << std::endl;

  for (size_t j = 0; j != nReactions(); ++j) {
    m_ropf[j] = std::max(0., -r(j));
    //m_ropf[j] = std::max(0., b(j));
    m_ropr[j] = std::max(0., r(j));
    //m_ropr[j] = std::max(0., -b(j));
    m_ropnet[j] = m_ropf[j] - m_ropr[j];
  }

  m_ROP_ok = true;
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

Eigen::SparseMatrix<double> Condensation::netRatesOfProgress_ddCi()
{
  // set rate constants, m_rfn
  _update_rates_T();

  // set activity concentrations, m_actConc
  _update_rates_C();

  // forward reaction rate coefficients
  Eigen::SparseMatrix<double> jac(nReactions(), nTotalSpecies());

  std::fill(m_rbuf0.begin(), m_rbuf0.end(), 1.0);
  m_reactantStoich.multiply(m_actConc.data(), m_rbuf0.data());

  for (size_t j = 0; j < nReactions(); ++j) {
    // inactive reactions
    if (m_rfn[j] < 0.0) {
      continue;
    }

    double ss = m_rbuf0[j] - m_rfn[j];
    size_t order = m_reactions[j]->reactants.size();
    auto& R = m_reactions[j];

    if (order == 1) {
      size_t ix = kineticsSpeciesIndex(R->reactants.begin()->first);
      size_t iy = kineticsSpeciesIndex(R->products.begin()->first);
      double y = m_actConc[iy];

      if (ss > 0. || (ss < 0. && y > - ss)) {
        jac.coeffRef(j, ix) = 1.;
      } else {
        jac.coeffRef(j, iy) = -1.;
      }
    } else { // order 2
      size_t ix1 = kineticsSpeciesIndex(R->reactants.begin()->first);
      size_t ix2 = kineticsSpeciesIndex(next(R->reactants.begin())->first);
      size_t iy = kineticsSpeciesIndex(R->products.begin()->first);

      double x1 = m_actConc[ix1];
      double x2 = m_actConc[ix2];
      double y = m_actConc[iy];
      double react = x1 + x2;
      double delta = (react - sqrt(react * react - 4 * ss)) / 2.;
      if (ss > 0. || (ss < 0. && y > - delta)) {
        jac.coeffRef(j, ix1) = (1. - (x1 - x2) / sqrt(react * react - 4 * ss)) / 2.;
        jac.coeffRef(j, ix2) = (1. - (x2 - x1) / sqrt(react * react - 4 * ss)) / 2.;
      } else {  // y < -delta
        jac.coeffRef(j, iy) = -1.;
      }
    }
  }

  return jac;
}

}
