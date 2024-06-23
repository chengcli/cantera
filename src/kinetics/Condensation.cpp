// C/C++
#include <numeric>
#include <algorithm>

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Condensation.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/IdealMoistPhase.h"

namespace Cantera
{

// x -> y at constant volume (mole concentration)
inline double satfunc1v(double s, double x, double y)
{
  double rate = x - s;
  if (rate > 0. || (rate < 0. && y > - rate)) {
    return rate;
  }
  return -y;
}

inline void set_jac1v(
    Eigen::SparseMatrix<double> &jac, double s, double const *conc,
    int j, int ix, int iy)
{
  double const& x = conc[ix];
  double const& y = conc[iy];
  double rate = x - s;
  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac.coeffRef(j, ix) = 1.;
  } else {
    jac.coeffRef(j, iy) = -1.;
  }
}

// x1 + x2 -> y at constant volume (mole concentration)
inline double satfunc2v(double s, double x1, double x2, double y)
{
  double delta = (x1 - x2) * (x1 - x2) + 4 * s;
  double rate = (x1 + x2 - sqrt(delta)) / 2.;

  if (rate > 0. || (rate < 0. && y > - rate)) {
    return rate;
  }
  return -y;
}

inline void set_jac2v(
    Eigen::SparseMatrix<double> &jac, double s, double const *conc,
    int j, int ix1, int ix2, int iy)
{
  double const& x1 = conc[ix1];
  double const& x2 = conc[ix2];
  double const& y = conc[iy];
  double delta = (x1 - x2) * (x1 - x2) + 4 * s;
  double rate = (x1 + x2 - sqrt(delta)) / 2.;

  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac.coeffRef(j, ix1) = (1. - (x1 - x2) / sqrt(delta)) / 2.;
    jac.coeffRef(j, ix2) = (1. - (x2 - x1) / sqrt(delta)) / 2.;
  } else {
    jac.coeffRef(j, iy) = -1.;
  }
}

// x -> y at constant pressure (mole fraction)
inline double satfunc1p(double s, double x, double y, double g)
{
  // boil all condensates
  if (s > 1.) {
    return -y;
  }

  std::cout << "s = " << s << ", x = " << x << ", y = " << y << ", g = " << g << std::endl;

  double rate = (x - s * g) / (1. - s);

  if (rate > 0. || (rate < 0. && y > - rate)) {
    return rate;
  }
  return -y;
}

inline void set_jac1p(Eigen::SparseMatrix<double> &jac, 
                      double s, double const* frac, double g,
                      int j, int ix, int iy)
{
  // boil all condensates
  if (s > 1.) {
    jac.coeffRef(j, iy) = -1.;
    return;
  }

  double const& x = frac[ix];
  double const& y = frac[iy];

  double rate = (x - s * g) / (1. - s);

  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac.coeffRef(j, ix) = 1.;
  } else {
    jac.coeffRef(j, iy) = -1.;
  }
}

// x1 + x2 -> y at constant pressure (mole fraction)
inline double satfunc2p(double s, double x1, double x2, double y, double g)
{
  double delta = (x1 - x2) * (x1 - x2) + 4 * s * (g - 2. * x1) * (g - 2. * x2);
  double rate = (x1 + x2 - 4 * g * s - sqrt(delta)) / (2. * (1. - 4. * s));

  if (rate > 0. || (rate < 0. && y > - rate)) {
    return rate;
  }
  return -y;
}

inline void set_jac2p(Eigen::SparseMatrix<double> &jac,
                      double s, double const* frac, double g,
                      int j, int ix1, int ix2, int iy)
{
  double const& x1 = frac[ix1];
  double const& x2 = frac[ix2];
  double const& y = frac[iy];

  double delta = (x1 - x2) * (x1 - x2) + 4 * s * (g - 2. * x1) * (g - 2. * x2);
  double rate = (x1 + x2 - 4 * g * s - sqrt(delta)) / (2. * (1. - 4. * s));

  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac.coeffRef(j, ix1) = (1. - (x1 - x2) / sqrt(delta)) / 2.;
    jac.coeffRef(j, ix2) = (1. - (x2 - x1) / sqrt(delta)) / 2.;
  } else {
    jac.coeffRef(j, iy) = -1.;
  }
}

inline Eigen::VectorXd linear_solve_rop(
    Eigen::SparseMatrix<double> const& jac,
    Eigen::SparseMatrix<double> const& stoich,
    Eigen::VectorXd const& b)
{
  auto A = jac * stoich;

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
  solver.compute(A);
  Eigen::VectorXd r = solver.solve(b);

  if (solver.info() != Eigen::Success) {
    throw CanteraError("Condensation::updateROP",
                       "Failed to solve for net rates of progress.");
  }

  std::cout << jac << std::endl;
  std::cout << A << std::endl;
  std::cout << A.transpose() * A << std::endl;
  std::cout << "b = " << b << std::endl;
  std::cout << "r = " << -r << std::endl;

  return r;
}

void Condensation::resizeReactions()
{
  Kinetics::resizeReactions();
  m_rbuf0.resize(nReactions());
  m_rbuf1.resize(nReactions());

  for (auto& rates : m_interfaceRates) {
    rates->resize(nTotalSpecies(), nReactions(), nPhases());
  }
  m_jac.resize(nReactions(), nTotalSpecies());
}

void Condensation::resizeSpecies()
{
  size_t kOld = m_kk;
  Kinetics::resizeSpecies();
  if (m_kk != kOld && nReactions()) {
      throw CanteraError("InterfaceKinetics::resizeSpecies", "Cannot add"
          " species to InterfaceKinetics after reactions have been added.");
  }

  m_conc.resize(m_kk);
}

void Condensation::getActivityConcentrations(double* const conc)
{
  if (m_use_mole_fraction) {
    _update_rates_X(conc);
  } else {
    _update_rates_C(conc);
  }
}

void Condensation::getFwdRateConstants(double* kfwd)
{
  _update_rates_T(kfwd);
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

  if (rtype == "nucleation") {
    if (r_base->reactants.size() == 1) {
      m_jxy.push_back(i);
    } else if (r_base->reactants.size() == 2) {
      m_jxxy.push_back(i);
    }
  } else if (rtype == "freezing") {
    m_jyy.push_back(i);
  } else {
    throw CanteraError("Condensation::addReaction",
                       "Unknown reaction type '{}'", rtype);
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
  _update_rates_T(m_rfn.data());
  if (m_use_mole_fraction) {
    _update_rates_X(m_conc.data());
  } else {
    _update_rates_C(m_conc.data());
  }

  if (m_ROP_ok) {
    return;
  }

  m_jac.setZero();

  double pres = thermo().pressure();
  double temp = thermo().temperature();
  double dens = pres / (GasConstant * temp);
  double xgas = m_conc[0];

  if (m_use_mole_fraction) {
    size_t nvapor = static_cast<IdealMoistPhase&>(thermo()).nVapor();
    for (size_t i = 1; i <= nvapor; i++) {
      xgas += m_conc[i];
    }
  }
  std::cout << "xgas = " << xgas << std::endl;

  Eigen::VectorXd b(nReactions());
  Eigen::SparseMatrix<double> stoich(m_stoichMatrix);

  // nucleation: x <=> y
  for (auto j : m_jxy) {
    std::cout << "jxy = " << j << std::endl;
    // inactive reactions
    if (m_rfn[j] < 0.0) {
      b(j) = 0.0;
      for (int i = 0; i < nTotalSpecies(); i++) {
        stoich.coeffRef(i,j) = 0.0;
      }
      continue;
    }

    auto& R = m_reactions[j];
    size_t ix = kineticsSpeciesIndex(R->reactants.begin()->first);
    size_t iy = kineticsSpeciesIndex(R->products.begin()->first);

    if (m_use_mole_fraction) {
      b(j) = satfunc1p(m_rfn[j] / dens, m_conc[ix], m_conc[iy], xgas);
      set_jac1p(m_jac, m_rfn[j] / dens, m_conc.data(), xgas, j, ix, iy);
    } else {
      b(j) = satfunc1v(m_rfn[j], m_conc[ix], m_conc[iy]);
      set_jac1v(m_jac, m_rfn[j], m_conc.data(), j, ix, iy);
    }
  }

  // nucleation: x1 + x2 <=> y
  for (auto j : m_jxxy) {
    std::cout << "jxxy = " << j << std::endl;
    // inactive reactions
    if (m_rfn[j] < 0.0) {
      b(j) = 0.0;
      for (int i = 0; i < nTotalSpecies(); i++) {
        stoich.coeffRef(i,j) = 0.0;
      }
      continue;
    }

    auto& R = m_reactions[j];
    size_t ix1 = kineticsSpeciesIndex(R->reactants.begin()->first);
    size_t ix2 = kineticsSpeciesIndex(next(R->reactants.begin())->first);
    size_t iy = kineticsSpeciesIndex(R->products.begin()->first);

    if (m_use_mole_fraction) {
      b(j) = satfunc2p(m_rfn[j] / (dens * dens), m_conc[ix1], m_conc[ix2], m_conc[iy], xgas);
      set_jac2p(m_jac, m_rfn[j] / (dens * dens), m_conc.data(), xgas, j, ix1, ix2, iy);
    } else {
      b(j) = satfunc2v(m_rfn[j], m_conc[ix1], m_conc[ix2], m_conc[iy]);
      set_jac2v(m_jac, m_rfn[j], m_conc.data(), j, ix1, ix2, iy);
    }
  }

  // freezing: y1 <=> y2
  for (auto j : m_jyy) {
    auto& R = m_reactions[j];
    size_t iy1 = kineticsSpeciesIndex(R->reactants.begin()->first);
    size_t iy2 = kineticsSpeciesIndex(R->products.begin()->first);
    
    if (temp > m_rfn[j]) { // higher than freezing temperature
      if (m_conc[iy2] > 0.) {
        b(j) = - m_conc[iy2];
        m_jac.coeffRef(j, iy2) = -1.;
      }
    } else { // lower than freezing temperature
      if (m_conc[iy1] > 0.) {
        b(j) = m_conc[iy1];
        m_jac.coeffRef(j, iy1) = 1.;
      }
    }
  }

  // solve the optimal net rates
  auto r = linear_solve_rop(m_jac, stoich, b);

  for (size_t j = 0; j != nReactions(); ++j) {
    m_ropf[j] = std::max(0., -r(j));
    m_ropr[j] = std::max(0., r(j));
    m_ropnet[j] = m_ropf[j] - m_ropr[j];
  }

  m_ROP_ok = true;
}

void Condensation::_update_rates_T(double *pdata)
{
  // Go find the temperature from the surface
  double T = thermo().temperature();

  if (T != m_temp) {
    m_temp = T;
    m_ROP_ok = false;
  }

  // loop over interface MultiRate evaluators for each reaction type
  for (auto& rates : m_interfaceRates) {
    bool changed = rates->update(thermo(), *this);
    if (changed) {
      rates->getRateConstants(pdata);
      m_ROP_ok = false;
    }
  }
}

void Condensation::_update_rates_C(double *pdata)
{
  thermo().getActivityConcentrations(pdata);
  m_ROP_ok = false;
}

void Condensation::_update_rates_X(double *pdata)
{
  thermo().getMoleFractions(pdata);
  m_ROP_ok = false;
}

Eigen::SparseMatrix<double> Condensation::netRatesOfProgress_ddX()
{
  updateROP();
  return m_jac;
}

Eigen::SparseMatrix<double> Condensation::netRatesOfProgress_ddCi()
{
  updateROP();
  return m_jac;
}

}
