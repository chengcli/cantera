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
  return {-y, 0.};
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
inline pair<double, double> satfunc2v(double s, double x1, double x2, double y,
                                      double logs_ddT = 0.)
{
  double delta = (x1 - x2) * (x1 - x2) + 4 * s;
  double rate = (x1 + x2 - sqrt(delta)) / 2.;

  if (rate > 0. || (rate < 0. && y > - rate)) {
    return {rate, - s * logs_ddT / sqrt(delta)};
  }
  return {-y, 0.};
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

  //std::cout << "s = " << s << ", x = " << x << ", y = " << y << ", g = " << g << std::endl;

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

bool Condensation::addReaction(shared_ptr<Reaction> r, bool resize)
{
  bool added = BulkKinetics::addReaction(r, resize);
  if (!added) {
    return false;
  }

  string rtype = r->rate()->subType();
  if (rtype == "") {
    rtype = r->rate()->type();
  }

  if (rtype == "nucleation") {
    if (r->reactants.size() == 1) {
      m_jxy.push_back(nReactions());
    } else if (r->reactants.size() == 2) {
      m_jxxy.push_back(nReactions());
    }
  } else if (rtype == "freezing") {
    m_jyy.push_back(nReactions());
  }

  return true;
}

void Condensation::updateROP() {
  BulkKinetics::updateROP();

  /// override rate of progress for thermodynamic reactions

  //! This variable has two interpretations.
  //! If m_use_mole_fraction is true, then it is the vector of mole fractions.
  //! If m_use_mole_fraction is false, then it is the vector of concentrations.
  Eigen::VectorXd m_conc(m_kk);
  Eigen::VectorXd m_intEng(m_kk);
  Eigen::VectorXd m_cv(m_kk);

  size_t nfast = m_jxy.size() + m_jxxy.size() + m_jyy.size();

  //! rate jacobian matrix
  Eigen::SparseMatrix<double> m_jac(nfast, nTotalSpecies());

  //! rate jacobian with respect to temperature
  vector<double> m_rfn_ddT(nfast);

  for (size_t i = 0; i < nfast; ++i) {
    m_bulk_rates[i]->processRateConstants_ddT(m_rfn_ddT.data(), nullptr, 0.);
  }

  if (m_use_mole_fraction) {
    thermo().getMoleFractions(m_conc.data());
  } else {
    thermo().getActivityConcentrations(m_conc.data());
    thermo().getIntEnergy_RT(m_intEng.data());
    thermo().getCv_R(m_cv.data());

    for (size_t i = 0; i < m_kk; i++) {
      m_intEng[i] *= thermo().temperature();
    }
  }

  m_jac.setZero();

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
  Eigen::SparseMatrix<double> stoich(nTotalSpecies(), nfast);
  Eigen::SparseMatrix<double> rate_ddT(nfast, nTotalSpecies());

  b.setZero();
  b_ddT.setZero();
  rate_ddT.setZero();

  // nucleation: x <=> y
  for (auto j : m_jxy) {
    // inactive reactions
    if (m_rfn[j] < 0.0) continue;
    for (size_t i = 0; i < nTotalSpecies(); ++i)
      stoich.coeffRef(i,j) = m_stoichMatrix.coeffRef(i,j);

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
    for (size_t i = 0; i < nTotalSpecies(); ++i)
      stoich.coeffRef(i,j) = m_stoichMatrix.coeffRef(i,j);

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
    for (size_t i = 0; i < nTotalSpecies(); ++i)
      stoich.coeffRef(i,j) = m_stoichMatrix.coeffRef(i,j);

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

  // set up temperature gradient
  if (!m_use_mole_fraction) {
    for (size_t j = 0; j < nfast; ++j) {
      // active reactions
      if (m_rfn[j] > 0. && b_ddT[j] != 0.0)  {
        for (size_t i = 0; i != nTotalSpecies(); ++i)
          rate_ddT.coeffRef(j, i) = b_ddT[j] * m_intEng[i];
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
}
}
