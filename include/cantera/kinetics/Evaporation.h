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
  string m_formula_str = "formula";

  //! returns s = svp(T)/RT or svp(T)/(RT)^2
  std::function<double(double)> m_svp;

  //! returns d(log(s))/dT
  std::function<double(double)> m_logsvp_ddT;

  size_t m_order = 1;

  double m_t3 = 0.0;
  double m_p3 = 0.0;
  double m_beta = 0.0;
  double m_delta = 0.0;

  double m_min_temp = 0.0;
  double m_max_temp = 1.0e30;
};

}

#endif
