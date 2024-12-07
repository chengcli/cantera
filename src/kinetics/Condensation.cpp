// C/C++
#include <numeric>
#include <algorithm>

#include <Eigen/Dense>

#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Condensation.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/IdealMoistPhase.h"

namespace Cantera
{

// x -> y at constant volume (mole concentration)
inline pair<double, double> satfunc1v(double s, double x, double y,
                                      double logs_ddT = 0.)
{
  double rate = x - s;
  if (rate > 0. || (rate < 0. && y > - rate)) {
    return {rate, - s * logs_ddT};
  }
  return {-y, - s * logs_ddT};
}

inline void set_jac1v(
    Eigen::MatrixXd &jac, double s, double const *conc,
    int j, int ix, int iy)
{
  double const& x = conc[ix];
  double const& y = conc[iy];
  double rate = x - s;
  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac(j, ix) = 1.;
  } else {
    jac(j, iy) = -1.;
  }
}

// x1 + x2 -> y at constant volume (mole concentration)
inline pair<double, double> satfunc2v(double s, double x1, double x2, double y,
                                      double logs_ddT = 0.)
{
  double delta = (x1 - x2) * (x1 - x2) + 4 * s;
  double rate = (x1 + x2 - sqrt(delta)) / 2.;

  if (rate > 0. || (rate < 0. && y > - rate)) {
    return {rate, - s * logs_ddT / sqrt(delta)};
  }
  return {-y, - s * logs_ddT / sqrt(delta)};
}

inline void set_jac2v(
    Eigen::MatrixXd &jac, double s, double const *conc,
    int j, int ix1, int ix2, int iy)
{
  double const& x1 = conc[ix1];
  double const& x2 = conc[ix2];
  double const& y = conc[iy];
  double delta = (x1 - x2) * (x1 - x2) + 4 * s;
  double rate = (x1 + x2 - sqrt(delta)) / 2.;

  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac(j, ix1) = (1. - (x1 - x2) / sqrt(delta)) / 2.;
    jac(j, ix2) = (1. - (x2 - x1) / sqrt(delta)) / 2.;
  } else {
    jac(j, iy) = -1.;
  }
}

inline double handle_singularity(double s)
{
  if (s - 1 < 1.e-8 && s - 1 > -1.e-8) {
    auto sgn = s > 1. ? 1. : -1.;
    return 1. + sgn * 1.e-6;
  }
  return s;
}

// x -> y at constant pressure (mole fraction)
inline double satfunc1p(double s, double x, double y, double g)
{
  // boil all condensates
  if (s > 1.) {
    return -y;
  }
  //s = handle_singularity(s);

  //std::cout << "s = " << s << ", x = " << x << ", y = " << y << ", g = " << g << std::endl;

  double rate = (x - s * g) / (1. - s);
  
  if (rate > 0. || (rate < 0. && y > - rate)) {
    return rate;
  }
  return -y;
}

inline void set_jac1p(Eigen::MatrixXd &jac, 
                      double s, double const* frac, double g,
                      int j, int ix, int iy)
{
  // boil all condensates
  if (s > 1.) {
    jac(j, iy) = -1.;
    return;
  }
  //s = handle_singularity(s);

  double const& x = frac[ix];
  double const& y = frac[iy];

  double rate = (x - s * g) / (1. - s);

  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac(j, ix) = 1.;
  } else {
    jac(j, iy) = -1.;
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

inline void set_jac2p(Eigen::MatrixXd &jac,
                      double s, double const* frac, double g,
                      int j, int ix1, int ix2, int iy)
{
  double const& x1 = frac[ix1];
  double const& x2 = frac[ix2];
  double const& y = frac[iy];

  double delta = (x1 - x2) * (x1 - x2) + 4 * s * (g - 2. * x1) * (g - 2. * x2);
  double rate = (x1 + x2 - 4 * g * s - sqrt(delta)) / (2. * (1. - 4. * s));

  if (rate > 0. || (rate < 0. && y > - rate)) {
    jac(j, ix1) = (1. - (x1 - x2) / sqrt(delta)) / 2.;
    jac(j, ix2) = (1. - (x2 - x1) / sqrt(delta)) / 2.;
  } else {
    jac(j, iy) = -1.;
  }
}

void Condensation::resizeReactions()
{
  Kinetics::resizeReactions();

  for (auto& rates : m_interfaceRates) {
    rates->resize(nTotalSpecies(), nReactions(), nPhases());
  }
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
  m_intEng.resize(m_kk);
  m_cv.resize(m_kk);
}

void Condensation::getActivityConcentrations(double* const pdata)
{
  if (m_use_mole_fraction) {
    _update_rates_X(pdata);
  } else {
    _update_rates_C(pdata);
  }
}

void Condensation::getFwdRateConstants(double* kfwd)
{
  _update_rates_T(kfwd, nullptr);
}

void Condensation::getFwdRateConstants_ddT(double* kfwd)
{
  _update_rates_T(nullptr, kfwd);
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
  } else if (rtype == "evaporation") {
    m_jevap.push_back(i);
  } else {
    m_jcloud.push_back(i);
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
  size_t nfast = m_jxy.size() + m_jxxy.size() + m_jyy.size();

  //! rate jacobian matrix
  Eigen::MatrixXd m_jac(nfast, nTotalSpecies());
  m_jac.setZero();

  //! rate jacobian with respect to temperature
  vector<double> m_rfn_ddT(nReactions());
  m_rfn_ddT.assign(nfast, 0.);

  _update_rates_T(m_rfn.data(), m_rfn_ddT.data());
  if (m_use_mole_fraction) {
    _update_rates_X(m_conc.data());
  } else {
    _update_rates_C(m_conc.data());
  }

  if (m_ROP_ok) {
    return;
  }

  double pres = thermo().pressure();
  double temp = thermo().temperature();
  double dens = pres / (GasConstant * temp);
  double xgas = 0.;

  if (m_use_mole_fraction) {
    size_t ngas = static_cast<IdealMoistPhase&>(thermo()).nGas();
    for (size_t i = 0; i < ngas; i++)
      xgas += m_conc[i];
    //std::cout << "xgas = " << xgas << std::endl;
  }

  Eigen::VectorXd b(nfast);
  Eigen::VectorXd b_ddT(nfast);
  Eigen::MatrixXd stoich(nTotalSpecies(), nfast);
  Eigen::MatrixXd rate_ddT(nfast, nTotalSpecies());

  b.setZero();
  b_ddT.setZero();
  stoich.setZero();
  rate_ddT.setZero();

  // nucleation: x <=> y
  for (auto j : m_jxy) {
    // inactive reactions
    if (m_rfn[j] < 0.0) continue;
      for (int i = 0; i < nTotalSpecies(); i++)
        stoich(i,j) = m_stoichMatrix.coeffRef(i,j);

    auto& R = m_reactions[j];
    size_t ix = kineticsSpeciesIndex(R->reactants.begin()->first);
    size_t iy = kineticsSpeciesIndex(R->products.begin()->first);

    if (m_use_mole_fraction) {
      b(j) = satfunc1p(m_rfn[j] / dens, m_conc[ix], m_conc[iy], xgas);
      set_jac1p(m_jac, m_rfn[j] / dens, m_conc.data(), xgas, j, ix, iy);
    } else {
      auto result = satfunc1v(m_rfn[j], m_conc[ix], m_conc[iy], m_rfn_ddT[j]);
      b(j) = result.first;
      b_ddT(j) = result.second;
      set_jac1v(m_jac, m_rfn[j], m_conc.data(), j, ix, iy);
    }
  }

  // nucleation: x1 + x2 <=> y
  for (auto j : m_jxxy) {
    //std::cout << "jxxy = " << j << std::endl;
    // inactive reactions
    if (m_rfn[j] < 0.0) continue;
      for (int i = 0; i < nTotalSpecies(); i++)
        stoich(i,j) = m_stoichMatrix.coeffRef(i,j);

    auto& R = m_reactions[j];
    size_t ix1 = kineticsSpeciesIndex(R->reactants.begin()->first);
    size_t ix2 = kineticsSpeciesIndex(next(R->reactants.begin())->first);
    size_t iy = kineticsSpeciesIndex(R->products.begin()->first);

    if (m_use_mole_fraction) {
      b(j) = satfunc2p(m_rfn[j] / (dens * dens), m_conc[ix1], m_conc[ix2], m_conc[iy], xgas);
      set_jac2p(m_jac, m_rfn[j] / (dens * dens), m_conc.data(), xgas, j, ix1, ix2, iy);
    } else {
      auto result = satfunc2v(m_rfn[j], m_conc[ix1], m_conc[ix2], m_conc[iy], m_rfn_ddT[j]);
      b(j) = result.first;
      b_ddT(j) = result.second;
      set_jac2v(m_jac, m_rfn[j], m_conc.data(), j, ix1, ix2, iy);
    }
  }

  // freezing: y1 <=> y2
  for (auto j : m_jyy) {
    if (m_rfn[j] < 0.0) continue;
    for (int i = 0; i < nTotalSpecies(); i++)
      stoich(i,j) = m_stoichMatrix.coeffRef(i,j);

    auto& R = m_reactions[j];
    size_t iy1 = kineticsSpeciesIndex(R->reactants.begin()->first);
    size_t iy2 = kineticsSpeciesIndex(R->products.begin()->first);
    
    if (temp > m_rfn[j]) { // higher than freezing temperature
      if (m_conc[iy2] > 0.) {
        b(j) = - m_conc[iy2];
        m_jac(j, iy2) = -1.;
      }
    } else { // lower than freezing temperature
      if (m_conc[iy1] > 0.) {
        b(j) = m_conc[iy1];
        m_jac(j, iy1) = 1.;
      }
    }
  }

  // set up temperature gradient
  if (!m_use_mole_fraction) {
    for (size_t j = 0; j < nfast; ++j) {
      // active reactions
      if (m_rfn[j] > 0. && b_ddT[j] != 0.0)  {
        for (size_t i = 0; i != nTotalSpecies(); ++i)
          rate_ddT(j, i) = b_ddT[j] * m_intEng[i];
      }
    }

    /*std::cout << "u = " << m_intEng << std::endl;
    std::cout << "cc = " << m_conc.dot(m_cv) << std::endl;
    std::cout << "cv = " << m_cv << std::endl;*/
    m_jac -= rate_ddT / m_conc.dot(m_cv);
  }

  // solve the optimal net rates
  Eigen::MatrixXd A = m_jac * stoich;
  Eigen::VectorXd r = A.colPivHouseholderQr().solve(b);

  /*std::cout << m_jac << std::endl;
  std::cout << A << std::endl;
  std::cout << A.transpose() * A << std::endl;
  std::cout << "b = " << b << std::endl;
  std::cout << "r = " << -r << std::endl;*/

  // scale rate down if some species becomes negative
  //Eigen::VectorXd rates = - stoich * r;
  
  for (size_t j = 0; j < nfast; ++j) {
    m_ropf[j] = std::max(0., -r(j));
    m_ropr[j] = std::max(0., r(j));
    m_ropnet[j] = m_ropf[j] - m_ropr[j];
  }

  for (size_t j = nfast; j < nReactions(); ++j) {
    m_ropf[j] = 0.;
    m_ropr[j] = 0.;
    m_ropnet[j] = 0.;
  }

  if (!m_use_mole_fraction) {
    // slow cloud reactions
    for (auto j : m_jcloud) {
      auto& R = m_reactions[j];
      size_t iy1 = kineticsSpeciesIndex(R->reactants.begin()->first);
      m_ropf[j] = m_rfn[j] * m_conc[iy1] * m_dt;
      m_ropr[j] = 0.;
      m_ropnet[j] = m_ropf[j];

      //for (int i = 0; i < nTotalSpecies(); i++)
      //  stoich(i,j) = m_stoichMatrix.coeffRef(i,j);
      //b(j) = (m_conc0[iy1] - m_conc[iy1]) / m_dt - m_rfn[j] * m_conc[iy1];
      //m_jac(j, iy1) = - 1. / m_dt - m_rfn[j];
    }

    // evaporation (only works for y <=> x)
    for (auto j : m_jevap) {
      auto& R = m_reactions[j];
      size_t iy1 = kineticsSpeciesIndex(R->reactants.begin()->first);
      auto vapors = static_cast<IdealMoistPhase&>(thermo()).vaporIndices(iy1);
      auto ivapor = vapors[0] - 1;

      // requires that the reaction indices and vapor indices are aligned
      if (m_boiling && (m_rfn[ivapor] > 1.)) {  // boiling point
        m_ropf[j] = m_conc[iy1];
        m_ropr[j] = 0.;
        m_ropnet[j] = m_ropf[j];
      } else {
        auto [rate, _] = satfunc1v(m_rfn[ivapor], m_conc[iy1], 0.);
        if (rate < 0.) {
          m_ropf[j] = - m_rfn[j] * rate * m_dt;
          m_ropr[j] = 0.;
          m_ropnet[j] = m_ropf[j];
        }
      }

      //for (int i = 0; i < nTotalSpecies(); i++)
      //  stoich(i,j) = m_stoichMatrix.coeffRef(i,j);
      //b(j) = (m_conc0[iy1] - m_conc[iy1]) / m_dt - m_rfn[j] * m_conc[iy1];
      //m_jac(j, iy1) = - 1. / m_dt - m_rfn[j];
    }
  }

  m_ROP_ok = true;
}

void Condensation::_update_rates_T(double *pdata, double *pdata_ddT)
{
  if (nReactions() == 0) {
    m_ROP_ok = true;
    return;
  }

  // Go find the temperature from the surface
  double T = thermo().temperature();

  if (T != m_temp) {
    m_temp = T;
    m_ROP_ok = false;
  }

  if (pdata_ddT != nullptr) {
    std::fill(pdata_ddT, pdata_ddT + nReactions(), 1.);
  }

  // loop over interface MultiRate evaluators for each reaction type
  for (auto& rates : m_interfaceRates) {
    bool changed = rates->update(thermo(), *this);
    if (changed) {
      if (pdata != nullptr) {
        rates->getRateConstants(pdata);
      }
      if (pdata_ddT != nullptr) {
        rates->processRateConstants_ddT(pdata_ddT, nullptr, 0.);
      }
      m_ROP_ok = false;
    }
  }
}

void Condensation::_update_rates_C(double *pdata)
{
  if (nReactions() == 0) {
    m_ROP_ok = true;
    return;
  }

  thermo().getActivityConcentrations(pdata);
  thermo().getIntEnergy_RT(m_intEng.data());
  thermo().getCv_R(m_cv.data());

  for (size_t i = 0; i < m_kk; i++) {
    m_intEng[i] *= thermo().temperature();
  }

  m_ROP_ok = false;
}

void Condensation::_update_rates_X(double *pdata)
{
  if (nReactions() == 0) {
    m_ROP_ok = true;
    return;
  }

  thermo().getMoleFractions(pdata);
  m_ROP_ok = false;
}

/*Eigen::SparseMatrix<double> Condensation::netRatesOfProgress_ddX()
{
  updateROP();
  return m_jac;
}

Eigen::SparseMatrix<double> Condensation::netRatesOfProgress_ddCi()
{
  updateROP();
  return m_jac;
}*/

}
