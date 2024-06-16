#include "cantera/numerics/Func1.h"
#include "cantera/kinetics/Condensation.h"

namespace Cantera
{

Condensation::Condensation(const AnyMap& node, const UnitStack& rate_units)
    : Condensation()
{
  setParameters(node, rate_units);

  if (!node.hasKey("rate-constant")) {
    throw CanteraError("Condensation::Condensation",
                       "Missing 'rate-constant' key in rate node.");
  }
  setRateParameters(node["equation"], node["rate-constant"], node.units(), rate_units);
}

void Condensation::setRateParameters(
    const AnyValue& equation,
    const AnyValue& rate,
    const UnitSystem& units, const UnitStack& rate_units)
{
  if (rate.empty()) {
    throw InputFileError("Condensation::setRateParameters", rate,
                         "Missing rate constant data.");
  }

  if (!rate.is<AnyMap>()) {
    throw InputFileError("Condensation::setRateParameters", rate,
                         "Expected a parameter map.");
  }

  auto& rate_map = rate.as<AnyMap>();
  string svp_name = rate_map[m_svp_str].asString();

  if (svp_name == "ideal") {
    m_t3 = rate_map["T3"].asDouble();
    m_p3 = rate_map["P3"].asDouble();
    m_beta = rate_map["beta"].asDouble();
    m_delta = rate_map["delta"].asDouble();
    m_svpfunc = [this](double T) {
      return m_p3 * exp((1. - m_t3 / T) * m_beta - m_delta * log(T / m_t3));
    };
  } else {
    //m_svpfunc = find_svp_function(reactants, svp_name);
  }

  m_valid = true;
}

void Condensation::validate(const string& equation, const Kinetics& kin)
{
  if (!m_valid) {
    throw InputFileError("Condensation::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
  }
}

double Condensation::evalFromStruct(const ArrheniusData& shared_data) const
{
  return m_svpfunc(shared_data.temperature);
}

void Condensation::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
  throw NotImplementedError("Condensation::getParameters",
                            "Not implemented by '{}' object.", type());
}

}
