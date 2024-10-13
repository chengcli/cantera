#include "cantera/numerics/Func1.h"
#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/Evaporation.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

EvaporationRate::EvaporationRate(const AnyMap& node, const UnitStack& rate_units)
    : EvaporationRate()
{
  setParameters(node, rate_units);

  if (!node.hasKey("rate-constant")) {
    throw CanteraError("EvaporationRate::EvaporationRate",
                       "Missing 'rate-constant' key in rate node.");
  }

  setRateParameters(node["equation"], node["rate-constant"], node);
}

void EvaporationRate::setRateParameters(
    const AnyValue& equation,
    const AnyValue& rate,
    const AnyMap& node)
{
  if (rate.empty()) {
    throw InputFileError("EvaporationRate::setRateParameters", rate,
                         "Missing rate constant data.");
  }

  if (!rate.is<AnyMap>()) {
    throw InputFileError("EvaporationRate::setRateParameters", rate,
                         "Expected a parameter map.");
  }

  auto& rate_map = rate.as<AnyMap>();

  if (rate.hasKey("A")) {
    m_A = rate_map["A"].asDouble();
  };

  if (rate.hasKey("b")) {
    m_b = rate_map["b"].asDouble();
  };

  if (rate.hasKey("Ea")) {
    m_Ea_R = rate_map["Ea"].asDouble();
  };

  m_valid = true;
}

void EvaporationRate::validate(const string& equation, const Kinetics& kin)
{
  if (!m_valid) {
    throw InputFileError("EvaporationRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
  }
}

double EvaporationRate::evalFromStruct(const ArrheniusData& shared_data) const
{
  double T = shared_data.temperature;
  return m_A * exp(m_b * log(T) - m_Ea_R / T);
}

//! Evaluate derivative of reaction rate with respect to temperature
//! divided by reaction rate
double EvaporationRate::ddTScaledFromStruct(const ArrheniusData& shared_data) const
{
  double T = shared_data.temperature;
  return (m_Ea_R / T + m_b) / T * m_A * exp(m_b * log(T) - m_Ea_R / T);
}

void EvaporationRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
  throw NotImplementedError("EvaporationRate::getParameters",
                            "Not implemented by '{}' object.", type());
}

}
