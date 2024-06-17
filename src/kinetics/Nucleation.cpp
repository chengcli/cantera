#include "cantera/numerics/Func1.h"
#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/Nucleation.h"
#include "cantera/kinetics/svp_funcs.h"
#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

string concatenate(const Composition& comp, char sep) {
  if (comp.empty()) return "";
  if (comp.size() == 1) return comp.begin()->first;
  
  std::string result = comp.begin()->first;

  for (auto it = std::next(comp.begin()); it != comp.end(); ++it) {
    result += sep + it->first;
  }

  return result;
}

std::function<double(double)> find_svp_function(const Composition& reactants,
                                                const string& svp_name) 
{
  string name = concatenate(reactants, '-') + '-' + svp_name;

  if (name == "NH3-H2S-lewis" || name == "H2S-NH3-lewis") {
    return svp_nh3_h2s_Lewis;
  } else if (name == "H2O-antoine") {
    return svp_h2o_Antoine;
  } else if (name == "NH3-antoine") {
    return svp_nh3_Antoine;
  }

  throw CanteraError("find_svp_function",
                     "No SVP function found for reaction '{}'.", name);
}

NucleationRate::NucleationRate(const AnyMap& node, const UnitStack& rate_units)
    : NucleationRate()
{
  setParameters(node, rate_units);

  if (!node.hasKey("rate-constant")) {
    throw CanteraError("NucleationRate::NucleationRate",
                       "Missing 'rate-constant' key in rate node.");
  }

  setRateParameters(node["equation"], node["rate-constant"], node);
}

void NucleationRate::setRateParameters(
    const AnyValue& equation,
    const AnyValue& rate,
    const AnyMap& node)
{
  if (rate.empty()) {
    throw InputFileError("NucleationRate::setRateParameters", rate,
                         "Missing rate constant data.");
  }

  if (!rate.is<AnyMap>()) {
    throw InputFileError("NucleationRate::setRateParameters", rate,
                         "Expected a parameter map.");
  }

  auto& rate_map = rate.as<AnyMap>();
  string svp_name = rate_map[m_svp_str].asString();

  if (rate.hasKey("minT")) {
    m_min_temp = rate_map["minT"].asDouble();
  };

  if (rate.hasKey("maxT")) {
    m_max_temp = rate_map["maxT"].asDouble();
  };

  if (svp_name == "ideal") {
    m_t3 = rate_map["T3"].asDouble();
    m_p3 = rate_map["P3"].asDouble();
    m_beta = rate_map["beta"].asDouble();
    m_delta = rate_map["delta"].asDouble();
    m_svpfunc = [this](double T) {
      return m_p3 * exp((1. - m_t3 / T) * m_beta - m_delta * log(T / m_t3));
    };
  } else {
    Reaction rtmp;
    parseReactionEquation(rtmp, equation.asString(), node, nullptr);
    m_svpfunc = find_svp_function(rtmp.reactants, svp_name);
  }

  m_valid = true;
}

void NucleationRate::validate(const string& equation, const Kinetics& kin)
{
  if (!m_valid) {
    throw InputFileError("NucleationRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
  }
}

double NucleationRate::evalFromStruct(const ArrheniusData& shared_data) const
{
  if (shared_data.temperature < m_min_temp || shared_data.temperature > m_max_temp) {
    return -1;
  }

  return m_svpfunc(shared_data.temperature);
}

void NucleationRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{
  throw NotImplementedError("NucleationRate::getParameters",
                            "Not implemented by '{}' object.", type());
}

}
