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
    m_dt = 0.;
  }

  void setQuantityConcentration(double dt = 0.0) {
    if (m_use_mole_fraction) {
      m_ROP_ok = false;
    }

    m_use_mole_fraction = false;
    m_dt = dt;
  }

  void setBoiling(bool boil) {
    m_boiling = boil;
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
  
  //! reaction indices for x <=> y
  vector<size_t> m_jxy;

  //! reaction indices for x1 + x2 <=> y
  vector<size_t> m_jxxy;

  //! reaction indices for freezing reaction x(l) <=> x(s)
  vector<size_t> m_jyy;

  //! reaction indices for slow cloud reactions
  vector<size_t> m_jcloud;

  //! reaction indices for evaporation
  vector<size_t> m_jevap;

  bool m_use_mole_fraction = false;
  bool m_ROP_ok = false;

  //! Current temperature of the data
  double m_temp = 0.0;

  //! Current time step
  double m_dt = 0.0;

  //! Boiling flag
  bool m_boiling = false;
};

}

#endif
