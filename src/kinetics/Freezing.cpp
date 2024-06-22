#include "cantera/numerics/Func1.h"
#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/Freezing.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

FreezingRate::FreezingRate(const AnyMap& node, const UnitStack& rate_units)
    : FreezingRate()
{
  setParameters(node, rate_units);

  if (!node.hasKey("rate-constant")) {
    throw CanteraError("FreezingRate::FreezingRate",
                       "Missing 'rate-constant' key in rate node.");
  }

  setRateParameters(node["equation"], node["rate-constant"], node);
}

void FreezingRate::setRateParameters(
    const AnyValue& equation,
    const AnyValue& rate,
    const AnyMap& node)
{
  if (rate.empty()) {
    throw InputFileError("FreezingRate::setRateParameters", rate,
                         "Missing rate constant data.");
  }

  if (!rate.is<AnyMap>()) {
    throw InputFileError("FreezingRate::setRateParameters", rate,
                         "Expected a parameter map.");
  }

  auto& rate_map = rate.as<AnyMap>();

  m_t3 = rate_map["T3"].asDouble();
  m_valid = true;
}

void FreezingRate::validate(const string& equation, const Kinetics& kin)
{
  if (!m_valid) {
    throw InputFileError("FreezingRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
  }
}

double FreezingRate::evalFromStruct(const ArrheniusData& shared_data) const
{
  return m_t3;
}

void FreezingRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
  throw NotImplementedError("FreezingRate::getParameters",
                            "Not implemented by '{}' object.", type());
}

}
