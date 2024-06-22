#include "cantera/numerics/Func1.h"
#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/Melting.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

MeltingRate::MeltingRate(const AnyMap& node, const UnitStack& rate_units)
    : MeltingRate()
{
  setParameters(node, rate_units);

  if (!node.hasKey("rate-constant")) {
    throw CanteraError("MeltingRate::MeltingRate",
                       "Missing 'rate-constant' key in rate node.");
  }

  setRateParameters(node["equation"], node["rate-constant"], node);
}

void MeltingRate::setRateParameters(
    const AnyValue& equation,
    const AnyValue& rate,
    const AnyMap& node)
{
  if (rate.empty()) {
    throw InputFileError("MeltingRate::setRateParameters", rate,
                         "Missing rate constant data.");
  }

  if (!rate.is<AnyMap>()) {
    throw InputFileError("MeltingRate::setRateParameters", rate,
                         "Expected a parameter map.");
  }

  auto& rate_map = rate.as<AnyMap>();

  m_t3 = rate_map["T3"].asDouble();
  m_valid = true;
}

void MeltingRate::validate(const string& equation, const Kinetics& kin)
{
  if (!m_valid) {
    throw InputFileError("MeltingRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
  }
}

double MeltingRate::evalFromStruct(const ArrheniusData& shared_data) const
{
  if (shared_data.temperature >= m_t3) {
    return -1;
  }

  return 0;
}

void MeltingRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
  throw NotImplementedError("MeltingRate::getParameters",
                            "Not implemented by '{}' object.", type());
}

}
