#ifndef CT_FREEZING_H
#define CT_FREEZING_H

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

class FreezingRate : public ReactionRate {
 public:
  FreezingRate() = default;
  FreezingRate(const AnyMap& node, const UnitStack& rate_units);

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<FreezingRate, ArrheniusData>>();
  }

  //! Set the rate parameters for this reaction.
  void setRateParameters(const AnyValue& equation,
                         const AnyValue& rate,
                         const AnyMap& node);

  //! return the rate coefficient type
  const string type() const override { return "freezing"; }

  void getParameters(AnyMap& rateNode, const Units& rate_units=Units(0.)) const;
  using ReactionRate::getParameters;

  void validate(const string& equation, const Kinetics& kin) override;

  double evalFromStruct(const ArrheniusData& shared_data) const;

  double ddTScaledFromStruct(const ArrheniusData& shared_data) const {
    return 0.;
  }

 protected:
  string m_formula_str = "formula";
  std::function<double(double)> m_meltfunc;

  double m_t3 = 0.0;
};

}

#endif
