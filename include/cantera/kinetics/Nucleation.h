#ifndef CT_NUCLEATION_H
#define CT_NUCLEATION_H

#include <functional>

#include "cantera/base/ct_defs.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/Arrhenius.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyValue;
class AnyMap;

class NucleationRate : public ReactionRate {
 public:
  NucleationRate() = default;
  NucleationRate(const AnyMap& node, const UnitStack& rate_units);

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<NucleationRate, ArrheniusData>>();
  }

  //! Set the rate parameters for this reaction.
  void setRateParameters(const AnyValue& equation,
                         const AnyValue& rate,
                         const AnyMap& node);

  //! return the rate coefficient type
  const string type() const override { return "condensation"; }

  void getParameters(AnyMap& rateNode, const Units& rate_units=Units(0.)) const;
  using ReactionRate::getParameters;

  void validate(const string& equation, const Kinetics& kin) override;

  double evalFromStruct(const ArrheniusData& shared_data) const;

 protected:
  std::function<double(double)> m_svpfunc;

  double m_t3 = 0.0;
  double m_p3 = 0.0;
  double m_beta = 0.0;
  double m_delta = 0.0;
  string m_svp_str = "formula";

  double m_min_temp = 0.0;
  double m_max_temp = 1.0e30;
};

std::function<double(double)> find_svp_function(const Composition& reactants, 
                                                const string& svp_name);

}

#endif
