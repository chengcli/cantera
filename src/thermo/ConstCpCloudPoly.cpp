
#include "cantera/base/AnyMap.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ConstCpCloudPoly.h"

namespace Cantera
{

ConstCpCloudPoly::ConstCpCloudPoly()
  : SpeciesThermoInterpType(0.0, std::numeric_limits<double>::infinity(), 0.0)
{}

ConstCpCloudPoly::ConstCpCloudPoly(
    double tlow, double thigh, double pref, const double *coeffs):
  SpeciesThermoInterpType(tlow, thigh, pref)
{
  setParameters(coeffs[0], coeffs[1], coeffs[2]);
}

void ConstCpCloudPoly::setParameters(double t0, double h0, double cp0)
{
  m_t0 = t0;
  m_h0_R = h0 / GasConstant;
  m_cp0_R = cp0 / GasConstant;
}

void ConstCpCloudPoly::setReaction(std::string name, Kinetics& kin) 
{
  Composition reactants;

  for (int j = 0; j < kin.nReactions(); j++) {
    auto rxn = kin.reaction(j);
    if (rxn->products.size() == 1 && rxn->products.begin()->first == name) 
    {
      m_rate = rxn->rate();
      reactants = rxn->reactants;
      break;
    }
  }

  if (reactants.empty()) {
    throw CanteraError("ConstCpCloudPoly::setReaction",
                       "No reaction found for species " + name);
  }

  auto& thermo = kin.thermo().speciesThermo();

  for (auto& [name, _] : reactants) {
    auto k = kin.kineticsSpeciesIndex(name);
    auto& thermo = kin.thermo().speciesThermo();
    m_vapor_thermo.push_back(thermo.getSpeciesThermo(k));
  }
}

void ConstCpCloudPoly::updatePropertiesTemp(const double temp,
                                   double* cp_R,
                                   double* h_RT,
                                   double* s_R) const
{
  *cp_R = m_cp0_R;
  *h_RT = (m_h0_R + (temp - m_t0) * m_cp0_R) / temp;

  double cpv_R;
  double hv_RT, sum_hv_RT = 0.;
  double sv_R, sum_sv_R = 0.;

  int order = m_vapor_thermo.size();
  double svp = m_rate->eval(temp) * pow(GasConstant * temp, order);
  sum_sv_R -= log(svp / m_Pref);

  for (int i = 0; i < order; i++) {
    m_vapor_thermo[i]->updatePropertiesTemp(temp, &cpv_R, &hv_RT, &sv_R);
    sum_sv_R += sv_R;
    sum_hv_RT += hv_RT;
  }

  *s_R = sum_sv_R - sum_hv_RT + (*h_RT);
}

void ConstCpCloudPoly::reportParameters(
    size_t& n, int& type, double& tlow, double& thigh,
    double& pref, double* const coeffs) const
{
  n = 0;
  type = CONSTANT_CP_CLOUD;
  tlow = m_lowT;
  thigh = m_highT;
  pref = m_Pref;
  coeffs[0] = m_t0;
  coeffs[1] = m_h0_R * GasConstant;
  coeffs[2] = m_cp0_R * GasConstant;
}

void ConstCpCloudPoly::getParameters(AnyMap& thermo) const
{
  thermo["model"] = "constant-cp-cloud";
  SpeciesThermoInterpType::getParameters(thermo);
  thermo["T0"].setQuantity(m_t0, "K");
  thermo["h0"].setQuantity(m_h0_R * GasConstant, "J/kmol");
  thermo["cp0"].setQuantity(m_cp0_R * GasConstant, "J/kmol/K");
}

}
