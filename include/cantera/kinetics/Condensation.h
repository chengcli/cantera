#ifndef CT_CONDENSATION_H
#define CT_CONDENSATION_H

#include "Kinetics.h"

namespace Cantera
{

class Condensation : public Kinetics {
 public:
  Condensation() = default;

  ~Condensation() override {}

  void resizeReactions() override;

  string kineticsType() const override {
    return "nucleation";
  }

  void getActivityConcentrations(double* const conc) override;

  bool isReversible(size_t i) override {
    return true;
  }

  void getFwdRateConstants(double* kfwd) override;

  void resizeSpecies() override;

  bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;

  void updateROP() override;

  void _update_rates_T();

  void _update_rates_C();

  Eigen::SparseMatrix<double> netRatesOfProgress_ddCi() override;

 protected:
  //! Process derivatives
  //! @param stoich  stoichiometry manager
  //! @param in  rate expression used for the derivative calculation
  //! @param ddX true: w.r.t mole fractions false: w.r.t species concentrations
  //! @return a sparse matrix of derivative contributions for each reaction of
  //! dimensions nTotalReactions by nTotalSpecies
  Eigen::SparseMatrix<double> calculateCompositionDerivatives(
      StoichManagerN& stoich, const vector<double>& in, bool ddX=true);

  vector<double> m_actConc;
  vector<double> m_satf;

  //! Number of dimensions of reacting phase (2 for InterfaceKinetics, 1 for
  //! EdgeKinetics)
  size_t m_nDim = 2;

  //! Vector of rate handlers for interface reactions
  vector<unique_ptr<MultiRateBase>> m_interfaceRates;
  map<string, size_t> m_interfaceTypes; //!< Rate handler mapping

  bool m_ROP_ok = false;

  //! Current temperature of the data
  double m_temp = 0.0;

  //! Buffers for partial rop results with length nReactions()
  vector<double> m_rbuf0;
  vector<double> m_rbuf1;
};

}

#endif
