#ifndef CT_CONDENSATION_H
#define CT_CONDENSATION_H

#include "Kinetics.h"

namespace Cantera
{

class Condensation : public Kinetics {
 public:
  Condensation() = default;

  ~Condensation() override {}

  void setQuantityMoleFraction() {
    if (!m_use_mole_fraction) {
      m_ROP_ok = false;
    }
    m_use_mole_fraction = true;
  }

  void setQuantityConcentration() {
    if (m_use_mole_fraction) {
      m_ROP_ok = false;
    }
    m_use_mole_fraction = false;
  }

  void resizeReactions() override;

  string kineticsType() const override {
    return "condensation";
  }

  void getActivityConcentrations(double* const pdata) override;

  bool isReversible(size_t i) override {
    return true;
  }

  void getFwdRateConstants(double* kfwd) override;

  void getFwdRateConstants_ddT(double* dkfwd) override;

  void resizeSpecies() override;

  bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;

  void updateROP() override;

  //Eigen::SparseMatrix<double> netRatesOfProgress_ddCi() override;
  //Eigen::SparseMatrix<double> netRatesOfProgress_ddX() override;

 protected:
  void _update_rates_T(double *pdata, double *pdata_ddT);
  void _update_rates_C(double *pdata);
  void _update_rates_X(double *pdata);

  //! Process derivatives
  //! @param stoich  stoichiometry manager
  //! @param in  rate expression used for the derivative calculation
  //! @param ddX true: w.r.t mole fractions false: w.r.t species concentrations
  //! @return a sparse matrix of derivative contributions for each reaction of
  //! dimensions nTotalReactions by nTotalSpecies
  Eigen::MatrixXd calculateCompositionDerivatives(
      StoichManagerN& stoich, const vector<double>& in, bool ddX=true);

  //! This variable has two interpretations.
  //! If m_use_mole_fraction is true, then it is the vector of mole fractions.
  //! If m_use_mole_fraction is false, then it is the vector of concentrations.
  Eigen::VectorXd m_conc;
  Eigen::VectorXd m_intEng;
  Eigen::VectorXd m_cv;

  //! Number of dimensions of reacting phase (2 for InterfaceKinetics, 1 for
  //! EdgeKinetics)
  size_t m_nDim = 2;

  //! Vector of rate handlers for interface reactions
  vector<unique_ptr<MultiRateBase>> m_interfaceRates;
  map<string, size_t> m_interfaceTypes; //!< Rate handler mapping
  
  //! reaction index for x <=> y
  vector<size_t> m_jxy;

  //! reaction index for x1 + x2 <=> y
  vector<size_t> m_jxxy;

  //! reaction index for freezing reaction x(l) <=> x(s)
  vector<size_t> m_jyy;

  //! rate jacobian matrix
  Eigen::MatrixXd m_jac;

  //! rate jacobian with respect to temperature
  vector<double> m_rfn_ddT;

  bool m_use_mole_fraction = false;
  bool m_ROP_ok = false;

  //! Current temperature of the data
  double m_temp = 0.0;

  //! Buffers for partial rop results with length nReactions()
  vector<double> m_rbuf0;
  vector<double> m_rbuf1;
};

}

#endif
