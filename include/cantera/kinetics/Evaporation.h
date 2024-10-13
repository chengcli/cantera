#ifndef CT_EVAPORATION_H
#define CT_EVAPORATION_H

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

class EvaporationRate : public ReactionRate {
 public:
  EvaporationRate() = default;
  EvaporationRate(const AnyMap& node, const UnitStack& rate_units);

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<EvaporationRate, ArrheniusData>>();
  }

  //! Set the rate parameters for this reaction.
  void setRateParameters(const AnyValue& equation,
                         const AnyValue& rate,
                         const AnyMap& node);

  //! return the rate coefficient type
  const string type() const override { return "evaporation"; }

  void getParameters(AnyMap& rateNode, const Units& rate_units=Units(0.)) const;
  using ReactionRate::getParameters;

  void validate(const string& equation, const Kinetics& kin) override;

  double evalFromStruct(const ArrheniusData& shared_data) const;

  double ddTScaledFromStruct(const ArrheniusData& shared_data) const;

 protected:
  double m_A = 0.; //!< Pre-exponential factor
  double m_b = 0.; //!< Temperature exponent
  double m_Ea_R = 0.; //!< Activation energy in units of R
};

}

#endif
