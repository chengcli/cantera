#ifndef CT_NUCLEATION_H
#define CT_NUCLEATION_H

#include "Kinetics.h"

namespace Cantera
{

class Nucleation : public Kinetics {
 public:
  Nucleation() = default;

  ~Nucleation() override {}

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

 protected:
  vector<double> m_actConc;

  //! Number of dimensions of reacting phase (2 for InterfaceKinetics, 1 for
  //! EdgeKinetics)
  size_t m_nDim = 2;

  //! Vector of rate handlers for interface reactions
  vector<unique_ptr<MultiRateBase>> m_interfaceRates;
  map<string, size_t> m_interfaceTypes; //!< Rate handler mapping

  bool m_ROP_ok = false;

  //! Current temperature of the data
  double m_temp = 0.0;
};

}

#endif
