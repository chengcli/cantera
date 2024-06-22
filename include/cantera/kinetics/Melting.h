#ifndef CT_MELTING_H
#define CT_MELTING_H

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

class MeltingRate : public ReactionRate {
 public:
  MeltingRate() = default;
  MeltingRate(const AnyMap& node, const UnitStack& rate_units);

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<MeltingRate, ArrheniusData>>();
  }

  //! Set the rate parameters for this reaction.
  void setRateParameters(const AnyValue& equation,
                         const AnyValue& rate,
                         const AnyMap& node);

  //! return the rate coefficient type
  const string type() const override { return "melting"; }

  void getParameters(AnyMap& rateNode, const Units& rate_units=Units(0.)) const;
  using ReactionRate::getParameters;

  void validate(const string& equation, const Kinetics& kin) override;

  double evalFromStruct(const ArrheniusData& shared_data) const;

 protected:
  string m_melt_str = "formula";
  std::function<double(double)> m_meltfunc;

  double m_t3 = 0.0;
};

}

#endif
