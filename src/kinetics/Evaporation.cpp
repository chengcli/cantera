#include "cantera/numerics/Func1.h"
#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/Evaporation.h"
#include "cantera/kinetics/svp_funcs.h"
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
  string svp_name = rate_map[m_formula_str].asString();

  if (rate.hasKey("minT")) {
    m_min_temp = rate_map["minT"].asDouble();
  };

  if (rate.hasKey("maxT")) {
    m_max_temp = rate_map["maxT"].asDouble();
  };

  Reaction rtmp;
  parseReactionEquation(rtmp, equation.asString(), node, nullptr);
  m_order = rtmp.reactants.size();

  if (svp_name == "ideal") {
    m_t3 = rate_map["T3"].asDouble();
    m_p3 = rate_map["P3"].asDouble();
    m_beta = rate_map["beta"].asDouble();
    m_delta = rate_map["delta"].asDouble();
    m_svp = [this](double T) {
      return m_p3 * exp((1. - m_t3 / T) * m_beta - m_delta * log(T / m_t3));
    };
    m_logsvp_ddT = [this](double T) {
      return m_beta * m_t3 / (T * T) - m_delta / T;
    };
    if (m_delta > 0.) {
      m_max_temp = std::min(m_max_temp, m_beta * m_t3 / m_delta);
    }
  } else {
    m_svp = find_svp(rtmp.reactants, svp_name);
    m_logsvp_ddT = find_logsvp_ddT(rtmp.reactants, svp_name);
  }

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
  if (T < m_min_temp || T > m_max_temp) {
    return -1;
  }

  double RT = GasConstant * T;
  return m_svp(T) / pow(RT, m_order);
}

double EvaporationRate::ddTScaledFromStruct(const ArrheniusData& shared_data) const
{
  double T = shared_data.temperature;
  if (T < m_min_temp || T > m_max_temp) {
    return 0.;
  }

  return m_logsvp_ddT(T) - m_order / T;
}

void EvaporationRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
  throw NotImplementedError("EvaporationRate::getParameters",
                            "Not implemented by '{}' object.", type());
}

}
